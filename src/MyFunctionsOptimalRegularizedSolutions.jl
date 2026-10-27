# Numerical functions for the NeurIPS 2026 LQMI model.
# Reference development and validation version: Julia 1.8.5.
# Control convention: u_t = L_t z_t; internal dynamics: A + B*L_t - K_t*H.
# T counts stored samples (1:T); K and L have T-1 gain matrices, and the
# terminal cost is x_T'*Q_T*x_T. The implementation assumes zero initial
# x/z cross-covariance for consistency with the TOD comparison and experiments;
# this is an initialization convention, not a restriction of the LCO framework.
# C and D accept a single matrix or nonempty collections of matrices.
# Only standard libraries are needed here. For the numerical GD benchmark,
# import Optim in the calling script and use the objective functions in (B.2).

# PACKAGES
using LinearAlgebra, Random, Statistics, Base.Threads

#(A): Functions to compute expected cost ----------------------------------------------------------------------
# Use this function to compute the expected cost with moments propagation.
"""
    expected_cost_raw_moments_propagation(T, m, p, k, K, L, A, B, H, C, D,
        x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q, R, Ω_ξ, Ω_ω, Ω_η)

Expected state/control cost from forward propagation of raw second moments.
`K` has size `(m,k,T-1)`, `L` has size `(p,m,T-1)`, `Q` has size `(m,m,T)`,
and `R` has size `(p,p,T)` (its terminal slice is unused). `Ω_ξ`, `Ω_ω`, and
`Ω_η` are covariance matrices. The implementation uses zero initial x/z
cross-covariance for consistency with the TOD comparison and experiments.
LCO itself allows correlated initialization: the general initial cross-moment
is `Cov(x_1,z_1) + x_1_mean*z_1_mean'`. Supporting that case would require an
explicit initial cross-moment and correlated Monte Carlo initial draws;
the current API does not implement it.
`C` and `D` accept single matrices or nonempty collections of matrices.
"""
function expected_cost_raw_moments_propagation(T, dimension_of_state, dimension_of_control, dimension_of_observation, K_matrix, L_matrix, A, B, H, C, D, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η)

    #Initial conditions for non-central moments
    S_t_xx = zeros(dimension_of_state, dimension_of_state, T)
    S_t_zz = zeros(dimension_of_state, dimension_of_state, T)
    S_t_xz = zeros(dimension_of_state, dimension_of_state, T)
    S_t_zx = zeros(dimension_of_state, dimension_of_state, T)
    S_t_xx[:,:,1] .= Σ_1_x[:,:] .+ x_1_mean*transpose(x_1_mean)
    S_t_zz[:,:,1] .= Σ_1_z[:,:] .+ z_1_mean*transpose(z_1_mean)
    S_t_xz[:,:,1] .= x_1_mean*transpose(z_1_mean) #initial values of state and state estimate are assumed to be uncorrelated
    S_t_zx[:,:,1] .= transpose(S_t_xz[:,:,1])

    for t in 2:T
        S_t_xx[:,:,t] .= A*S_t_xx[:,:,t-1]*transpose(A) .+ A*S_t_xz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1]) .+ B*L_matrix[:,:,t-1]*S_t_zx[:,:,t-1]*transpose(A) .+ B*L_matrix[:,:,t-1]*S_t_zz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1]) .+ sumCLMLC(C, L_matrix[:,:,t-1], S_t_zz[:,:,t-1]) .+ Ω_ξ
        S_t_zx[:,:,t] .= K_matrix[:,:,t-1]*H*S_t_xx[:,:,t-1]*transpose(A) .+ K_matrix[:,:,t-1]*H*S_t_xz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1]) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zx[:,:,t-1]*transpose(A) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1])
        S_t_xz[:,:,t] .= transpose(S_t_zx[:,:,t])
        S_t_zz[:,:,t] .= K_matrix[:,:,t-1]*H*S_t_xx[:,:,t-1]*transpose(K_matrix[:,:,t-1]*H) .+ sumKDMDK(D, K_matrix[:,:,t-1], S_t_xx[:,:,t-1]) .+ K_matrix[:,:,t-1]*H*S_t_xz[:,:,t-1]*transpose(A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zx[:,:,t-1]*transpose(K_matrix[:,:,t-1]*H) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zz[:,:,t-1]*transpose(A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H) .+ K_matrix[:,:,t-1]*Ω_ω*transpose(K_matrix[:,:,t-1]) .+ Ω_η
    end

    cost = 0
    for t in 1:T-1
      cost = cost + tr(Q_matrix[:,:,t]*S_t_xx[:,:,t]) + tr(transpose(L_matrix[:,:,t])*R_matrix[:,:,t]*L_matrix[:,:,t]*S_t_zz[:,:,t])
    end
    #final time step
    cost = cost + tr(Q_matrix[:,:,T]*S_t_xx[:,:,T])

    return cost

end

# Compute expected cost using the trace formula proved by induction.
"""
    expected_cost_using_trace_formula_induction(L, K, A, B, H, C, D, T, m,
        x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q, R, Ω_ξ, Ω_ω, Ω_η)

Expected cost from backward multipliers and additive-noise cost accumulation.
Uses the same objective, gain sign, horizon, and initialization assumptions as
`expected_cost_raw_moments_propagation`.
"""
function expected_cost_using_trace_formula_induction(L_matrix, K_matrix, A, B, H, C, D, T, dimension_of_state, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η)

    S_t_xx = zeros(dimension_of_state, dimension_of_state)
    S_t_zz = zeros(dimension_of_state, dimension_of_state)
    S_t_xz = zeros(dimension_of_state, dimension_of_state)
    S_t_zx = zeros(dimension_of_state, dimension_of_state)
    #Initial conditions for non-central moments
    S_t_xx[:,:] .= Σ_1_x[:,:] .+ x_1_mean*transpose(x_1_mean)
    S_t_zz[:,:] .= Σ_1_z[:,:] .+ z_1_mean*transpose(z_1_mean)
    S_t_xz[:,:] .= x_1_mean*transpose(z_1_mean) #initial values of state and state estimate are assumed to be uncorrelated
    S_t_zx[:,:] .= transpose(S_t_xz[:,:])

    # Variables used to compute the expected accumulated cost from time t to the final time T <c_t> = Tr[lambda_multiplier[:,:,t]*S_t^{xx} + omega_multiplier[:,:,t]*S_t^{zz} + nu_multiplier[:,:,t]*S_t^{xz}] + d_multiplier[t]
    # where S_t^{xx} = E[x_t*x_t^T], S_t^{zz} = E[z_t*z_t^T], S_t^{xz} = E[x_t*z_t^T]
    lambda_multiplier = zeros(dimension_of_state, dimension_of_state, T)
    omega_multiplier = zeros(dimension_of_state, dimension_of_state, T)
    nu_multiplier = zeros(dimension_of_state, dimension_of_state, T)
    d_multiplier = zeros(T)

    #multipliers at the last time step [boundary conditions]
    lambda_multiplier[:,:,T] .= Q_matrix[:,:,T]
    omega_multiplier[:,:,T] .= zeros(dimension_of_state, dimension_of_state)
    nu_multiplier[:,:,T] .= zeros(dimension_of_state, dimension_of_state)
    d_multiplier[T] = 0

    #computing the multipliers backwards in time [at fixed K and L]
    for t in (T-1):-1:1
        lambda_multiplier[:,:,t] .= Q_matrix[:,:,t] .+ transpose(A)*lambda_multiplier[:,:,t+1]*A .+ transpose(H)*transpose(K_matrix[:,:,t])*omega_multiplier[:,:,t+1]*K_matrix[:,:,t]*H .+ transpose(H)*transpose(K_matrix[:,:,t])*nu_multiplier[:,:,t+1]*A .+ sumDtGD(D, transpose(K_matrix[:,:,t])*omega_multiplier[:,:,t+1]*K_matrix[:,:,t])
        omega_multiplier[:,:,t] .= transpose(L_matrix[:,:,t])*R_matrix[:,:,t]*L_matrix[:,:,t] .+ transpose(L_matrix[:,:,t])*transpose(B)*lambda_multiplier[:,:,t+1]*B*L_matrix[:,:,t] .+ sumLtCtGCL(C, L_matrix[:,:,t], lambda_multiplier[:,:,t+1]) .+ transpose(A .+ B*L_matrix[:,:,t] .- K_matrix[:,:,t]*H)*omega_multiplier[:,:,t+1]*(A .+ B*L_matrix[:,:,t] .- K_matrix[:,:,t]*H) .+ transpose(A .+ B*L_matrix[:,:,t] .- K_matrix[:,:,t]*H)*nu_multiplier[:,:,t+1]*B*L_matrix[:,:,t]
        nu_multiplier[:,:,t] .= transpose(L_matrix[:,:,t])*transpose(B)*(lambda_multiplier[:,:,t+1] .+ transpose(lambda_multiplier[:,:,t+1]))*A .+ transpose(A .+ B*L_matrix[:,:,t] .- K_matrix[:,:,t]*H)*(omega_multiplier[:,:,t+1] .+ transpose(omega_multiplier[:,:,t+1]))*K_matrix[:,:,t]*H .+ transpose(A .+ B*L_matrix[:,:,t] .- K_matrix[:,:,t]*H)*nu_multiplier[:,:,t+1]*A .+ transpose(L_matrix[:,:,t])*transpose(B)*transpose(nu_multiplier[:,:,t+1])*K_matrix[:,:,t]*H
        d_multiplier[t] = d_multiplier[t+1] + tr(lambda_multiplier[:,:,t+1]*Ω_ξ .+ omega_multiplier[:,:,t+1]*Ω_η .+ omega_multiplier[:,:,t+1]*K_matrix[:,:,t]*Ω_ω*transpose(K_matrix[:,:,t]))
    end

    #Computing the expected cost at the first time step
    cost = tr(lambda_multiplier[:,:,1]*S_t_xx .+ omega_multiplier[:,:,1]*S_t_zz .+ nu_multiplier[:,:,1]*S_t_xz) + d_multiplier[1]

    return cost

end

#(B): Functions for Optimization -- Previous Methods ----------------------------------------------------------------------
#(B.1) (Todorov, 2005) ----------------------------------------------------------------------
function Todorov_optimization_multidimensional_case(A, B, H, C, D, T, dimension_of_state, dimension_of_control, dimension_of_observation, initial_x_state_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η, N_iter = 500)

    Σe_t = zeros(dimension_of_state,dimension_of_state)
    Σxe_t = zeros(dimension_of_state,dimension_of_state)
    Σx_t = zeros(dimension_of_state,dimension_of_state)

    Σe_t_start = zeros(dimension_of_state,dimension_of_state)
    Σxe_t_start = zeros(dimension_of_state,dimension_of_state)
    Σx_t_start = zeros(dimension_of_state,dimension_of_state)

    Σe_t_start[:,:] .= Σ_1_x[:,:] .+ Σ_1_z[:,:]
    Σxe_t_start[:,:] .= -Σ_1_z
    Σx_t_start[:,:] .= Σ_1_z .+ initial_x_state_mean*transpose(initial_x_state_mean)

    #OPTIMIZATION

    L = zeros(dimension_of_control,dimension_of_state,(T-1))
    K = zeros(dimension_of_state,dimension_of_observation,(T-1))

    for i in 1:N_iter

        #OPTIMAL CONTROLLER --> backwards through time

        #Initial (i.e. final) conditions
        Sx_tp1 = Q_matrix[:,:,T]
        Se_tp1 = zeros(dimension_of_state,dimension_of_state)

        for idx in (T-1):-1:1

            L[:,:,idx] = inv(R_matrix[:,:,idx] .+ transpose(B)*Sx_tp1*B .+ sumCtGC(C, (Sx_tp1 .+ Se_tp1)))*transpose(B)*Sx_tp1*A
            cSx_tp1 = Q_matrix[:,:,idx] .+ transpose(A)*Sx_tp1*(A .- B*L[:,:,idx]) .+ sumDtGD(D, transpose(K[:,:,idx])*Se_tp1*K[:,:,idx])
            Se_tp1 = transpose(A)*Sx_tp1*B*L[:,:,idx] .+ transpose(A .- K[:,:,idx]*H)*Se_tp1*(A .- K[:,:,idx]*H)
            Sx_tp1 = cSx_tp1

        end

        #OPTIMAL FILTER --> forward in time

        #Initial conditions --> given by function arguments
        Σe_t[:,:] .= Σe_t_start[:,:]
        Σxe_t[:,:] .= Σxe_t_start[:,:]
        Σx_t[:,:] .= Σx_t_start[:,:]

        for idx in 1:(T-1)

            K[:,:,idx] = A*Σe_t*transpose(H) * inv(H*Σe_t*transpose(H) .+ Ω_ω .+ sumDGDt(D, (Σe_t .+ Σx_t .+ Σxe_t .+ transpose(Σxe_t))))
            cΣe_t = Ω_ξ .+ Ω_η .+ (A .- K[:,:,idx]*H)*Σe_t*transpose(A) .+ sumCLMLC(C, L[:,:,idx], Σx_t)
            cΣx_t = Ω_η .+ K[:,:,idx]*H*Σe_t*transpose(A) .+ (A .- B*L[:,:,idx])*Σx_t*transpose(A .- B*L[:,:,idx]) .+ (A .- B*L[:,:,idx])*Σxe_t*transpose(H)*transpose(K[:,:,idx]) .+ K[:,:,idx]*H*transpose(Σxe_t)*transpose(A .- B*L[:,:,idx])
            cΣxe_t = (A .- B*L[:,:,idx])*Σxe_t*transpose(A .- K[:,:,idx]*H) .- Ω_η

            Σe_t[:,:]  .= cΣe_t
            Σx_t[:,:]  .= cΣx_t
            Σxe_t[:,:] .= cΣxe_t
        end

        # #first step for Σe_t, Σxe_t, Σx_t: the first K is zero, if we assume zero initial uncertainty on the state vector (Σ_1 = 0)
        # idx = 1
        # cΣe_t = Ω_ξ .+ Ω_η .+ (A .- K[:,:,idx]*H)*Σe_t*transpose(A) .+ sumCLMLC(C, L[:,:,idx], Σx_t)
        # Σx_t = Ω_η .+ K[:,:,idx]*H*Σe_t*transpose(A) .+ (A .- B*L[:,:,idx])*Σx_t*transpose((A .- B*L[:,:,idx])) .+ (A .- B*L[:,:,idx])*Σxe_t*transpose(H)*transpose(K[:,:,idx]) .+ K[:,:,idx]*H*transpose(Σxe_t)*transpose(A .- B*L[:,:,idx])
        # Σe_t = cΣe_t
        # Σxe_t = (A .- B*L[:,:,idx])*Σxe_t*transpose(A .- K[:,:,idx]*H) .- Ω_η

        # for idx in 2:(T-1)

        #     K[:,:,idx] = A*Σe_t*transpose(H)*inv(H*Σe_t*transpose(H) .+ Ω_ω .+ sumDGDt(D, (Σe_t .+ Σx_t .+ Σxe_t .+ transpose(Σxe_t))))
        #     cΣe_t = Ω_ξ .+ Ω_η .+ (A .- K[:,:,idx]*H)*Σe_t*transpose(A) .+ sumCLMLC(C, L[:,:,idx], Σx_t)
        #     Σx_t = Ω_η .+ K[:,:,idx]*H*Σe_t*transpose(A) .+ (A .- B*L[:,:,idx])*Σx_t*transpose((A .- B*L[:,:,idx])) .+ (A .- B*L[:,:,idx])*Σxe_t*transpose(H)*transpose(K[:,:,idx]) .+ K[:,:,idx]*H*transpose(Σxe_t)*transpose(A .- B*L[:,:,idx])
        #     Σe_t = cΣe_t
        #     Σxe_t = (A .- B*L[:,:,idx])*Σxe_t*transpose(A .- K[:,:,idx]*H) .- Ω_η

        # end

    end

    return .-L,K

end

#(B.2) (Damiani et al., 2024), Numerical GD  ----------------------------------------------------------------------
#function for GD optimization - OPTIMAL ESTIMATOR
function expected_cost_using_mom_prop_for_GD_ESTIMATION_optimization_raw_moments(x, L_matrix, A, B, H, C, D, T, dimension_of_state, dimension_of_control, dimension_of_observation, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η)

    K_matrix = reshape(x[1:dimension_of_state*dimension_of_observation*(T - 1)], dimension_of_state, dimension_of_observation, T - 1)
    L_matrix_tmp = zeros(dimension_of_control, dimension_of_state, T)
    L_matrix_tmp[:,:,1:T-1] .= L_matrix[:,:,:]

    #Initial conditions for non-central moments
    S_t_xx = zeros(dimension_of_state, dimension_of_state, T)
    S_t_zz = zeros(dimension_of_state, dimension_of_state, T)
    S_t_xz = zeros(dimension_of_state, dimension_of_state, T)
    S_t_zx = zeros(dimension_of_state, dimension_of_state, T)
    S_t_xx[:,:,1] .= Σ_1_x[:,:] .+ x_1_mean*transpose(x_1_mean)
    S_t_zz[:,:,1] .= Σ_1_z[:,:] .+ z_1_mean*transpose(z_1_mean)
    S_t_xz[:,:,1] .= x_1_mean*transpose(z_1_mean) #initial values of state and state estimate are assumed to be uncorrelated
    S_t_zx[:,:,1] .= transpose(S_t_xz[:,:,1])
    cost = 0
    cost = cost + tr(Q_matrix[:,:,1]*S_t_xx[:,:,1]) + tr(transpose(L_matrix[:,:,1])*R_matrix[:,:,1]*L_matrix[:,:,1]*S_t_zz[:,:,1])

    for t in 2:T
        S_t_xx[:,:,t] .= A*S_t_xx[:,:,t-1]*transpose(A) .+ A*S_t_xz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1]) .+ B*L_matrix[:,:,t-1]*S_t_zx[:,:,t-1]*transpose(A) .+ B*L_matrix[:,:,t-1]*S_t_zz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1]) .+ sumCLMLC(C, L_matrix[:,:,t-1], S_t_zz[:,:,t-1]) .+ Ω_ξ
        S_t_zx[:,:,t] .= K_matrix[:,:,t-1]*H*S_t_xx[:,:,t-1]*transpose(A) .+ K_matrix[:,:,t-1]*H*S_t_xz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1]) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zx[:,:,t-1]*transpose(A) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1])
        S_t_xz[:,:,t] .= transpose(S_t_zx[:,:,t])
        S_t_zz[:,:,t] .= K_matrix[:,:,t-1]*H*S_t_xx[:,:,t-1]*transpose(K_matrix[:,:,t-1]*H) .+ sumKDMDK(D, K_matrix[:,:,t-1], S_t_xx[:,:,t-1]) .+ K_matrix[:,:,t-1]*H*S_t_xz[:,:,t-1]*transpose(A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zx[:,:,t-1]*transpose(K_matrix[:,:,t-1]*H) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zz[:,:,t-1]*transpose(A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H) .+ K_matrix[:,:,t-1]*Ω_ω*transpose(K_matrix[:,:,t-1]) .+ Ω_η
        cost = cost + tr(Q_matrix[:,:,t]*S_t_xx[:,:,t]) + tr(transpose(L_matrix_tmp[:,:,t])*R_matrix[:,:,t]*L_matrix_tmp[:,:,t]*S_t_zz[:,:,t])
    end

    return cost

end

#function for GD optimization - OPTIMAL LINEAR CONTROLLER
function expected_cost_using_mom_prop_for_GD_CONTROL_optimization_raw_moments(x, K_matrix, A, B, H, C, D, T, dimension_of_state, dimension_of_control, dimension_of_observation, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η)

    L_matrix = reshape(x[1:dimension_of_control*dimension_of_state*(T - 1)], dimension_of_control, dimension_of_state, T - 1)
    L_matrix_tmp = zeros(dimension_of_control, dimension_of_state, T)
    L_matrix_tmp[:,:,1:T-1] .= L_matrix[:,:,:]

    #Initial conditions for non-central moments
    S_t_xx = zeros(dimension_of_state, dimension_of_state, T)
    S_t_zz = zeros(dimension_of_state, dimension_of_state, T)
    S_t_xz = zeros(dimension_of_state, dimension_of_state, T)
    S_t_zx = zeros(dimension_of_state, dimension_of_state, T)
    S_t_xx[:,:,1] .= Σ_1_x[:,:] .+ x_1_mean*transpose(x_1_mean)
    S_t_zz[:,:,1] .= Σ_1_z[:,:] .+ z_1_mean*transpose(z_1_mean)
    S_t_xz[:,:,1] .= x_1_mean*transpose(z_1_mean) #initial values of state and state estimate are assumed to be uncorrelated
    S_t_zx[:,:,1] .= transpose(S_t_xz[:,:,1])
    cost = 0
    cost = cost + tr(Q_matrix[:,:,1]*S_t_xx[:,:,1]) + tr(transpose(L_matrix[:,:,1])*R_matrix[:,:,1]*L_matrix[:,:,1]*S_t_zz[:,:,1])

    for t in 2:T
        S_t_xx[:,:,t] .= A*S_t_xx[:,:,t-1]*transpose(A) .+ A*S_t_xz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1]) .+ B*L_matrix[:,:,t-1]*S_t_zx[:,:,t-1]*transpose(A) .+ B*L_matrix[:,:,t-1]*S_t_zz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1]) .+ sumCLMLC(C, L_matrix[:,:,t-1], S_t_zz[:,:,t-1]) .+ Ω_ξ
        S_t_zx[:,:,t] .= K_matrix[:,:,t-1]*H*S_t_xx[:,:,t-1]*transpose(A) .+ K_matrix[:,:,t-1]*H*S_t_xz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1]) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zx[:,:,t-1]*transpose(A) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1])
        S_t_xz[:,:,t] .= transpose(S_t_zx[:,:,t])
        S_t_zz[:,:,t] .= K_matrix[:,:,t-1]*H*S_t_xx[:,:,t-1]*transpose(K_matrix[:,:,t-1]*H) .+ sumKDMDK(D, K_matrix[:,:,t-1], S_t_xx[:,:,t-1]) .+ K_matrix[:,:,t-1]*H*S_t_xz[:,:,t-1]*transpose(A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zx[:,:,t-1]*transpose(K_matrix[:,:,t-1]*H) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zz[:,:,t-1]*transpose(A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H) .+ K_matrix[:,:,t-1]*Ω_ω*transpose(K_matrix[:,:,t-1]) .+ Ω_η
        cost = cost + tr(Q_matrix[:,:,t]*S_t_xx[:,:,t]) + tr(transpose(L_matrix_tmp[:,:,t])*R_matrix[:,:,t]*L_matrix_tmp[:,:,t]*S_t_zz[:,:,t])
    end

    return cost

end

#function for GD optimization - OPTIMAL ESTIMATOR AND LINEAR CONTROLLER [joint optimization]
function expected_cost_using_mom_prop_for_GD_WHOLE_optimization_raw_moments(x, A, B, H, C, D, T, dimension_of_state, dimension_of_control, dimension_of_observation, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η)

    L_matrix = reshape(x[1:dimension_of_control*dimension_of_state*(T - 1)], dimension_of_control, dimension_of_state, T - 1)
    K_matrix = reshape(x[dimension_of_control*dimension_of_state*(T - 1)+1:dimension_of_control*dimension_of_state*(T - 1)+dimension_of_state*dimension_of_observation*(T - 1)], dimension_of_state, dimension_of_observation, T - 1)

    L_matrix_tmp = zeros(dimension_of_control, dimension_of_state, T)
    L_matrix_tmp[:,:,1:T-1] .= L_matrix[:,:,:]

    #Initial conditions for non-central moments
    S_t_xx = zeros(dimension_of_state, dimension_of_state, T)
    S_t_zz = zeros(dimension_of_state, dimension_of_state, T)
    S_t_xz = zeros(dimension_of_state, dimension_of_state, T)
    S_t_zx = zeros(dimension_of_state, dimension_of_state, T)
    S_t_xx[:,:,1] .= Σ_1_x[:,:] .+ x_1_mean*transpose(x_1_mean)
    S_t_zz[:,:,1] .= Σ_1_z[:,:] .+ z_1_mean*transpose(z_1_mean)
    S_t_xz[:,:,1] .= x_1_mean*transpose(z_1_mean) #initial values of state and state estimate are assumed to be uncorrelated
    S_t_zx[:,:,1] .= transpose(S_t_xz[:,:,1])
    cost = 0
    cost = cost + tr(Q_matrix[:,:,1]*S_t_xx[:,:,1]) + tr(transpose(L_matrix[:,:,1])*R_matrix[:,:,1]*L_matrix[:,:,1]*S_t_zz[:,:,1])

    for t in 2:T
        S_t_xx[:,:,t] .= A*S_t_xx[:,:,t-1]*transpose(A) .+ A*S_t_xz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1]) .+ B*L_matrix[:,:,t-1]*S_t_zx[:,:,t-1]*transpose(A) .+ B*L_matrix[:,:,t-1]*S_t_zz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1]) .+ sumCLMLC(C, L_matrix[:,:,t-1], S_t_zz[:,:,t-1]) .+ Ω_ξ
        S_t_zx[:,:,t] .= K_matrix[:,:,t-1]*H*S_t_xx[:,:,t-1]*transpose(A) .+ K_matrix[:,:,t-1]*H*S_t_xz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1]) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zx[:,:,t-1]*transpose(A) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zz[:,:,t-1]*transpose(B*L_matrix[:,:,t-1])
        S_t_xz[:,:,t] .= transpose(S_t_zx[:,:,t])
        S_t_zz[:,:,t] .= K_matrix[:,:,t-1]*H*S_t_xx[:,:,t-1]*transpose(K_matrix[:,:,t-1]*H) .+ sumKDMDK(D, K_matrix[:,:,t-1], S_t_xx[:,:,t-1]) .+ K_matrix[:,:,t-1]*H*S_t_xz[:,:,t-1]*transpose(A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zx[:,:,t-1]*transpose(K_matrix[:,:,t-1]*H) .+ (A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H)*S_t_zz[:,:,t-1]*transpose(A .+ B*L_matrix[:,:,t-1] .- K_matrix[:,:,t-1]*H) .+ K_matrix[:,:,t-1]*Ω_ω*transpose(K_matrix[:,:,t-1]) .+ Ω_η
        cost = cost + tr(Q_matrix[:,:,t]*S_t_xx[:,:,t]) + tr(transpose(L_matrix_tmp[:,:,t])*R_matrix[:,:,t]*L_matrix_tmp[:,:,t]*S_t_zz[:,:,t])
    end

    return cost

end

#(C): Functions for Optimization -- Lagrangian Constrained Optimization and Analytical GD ----------------------------------------------------------------------
# Compute Lagrange multipliers backwards in time at fixed K and L.
"""
    Compute_Lagrange_Multipliers(T, K, L, A, B, H, C, D, Q, R)

Return `(λ, ω, ν)`, each of size `(m,m,T)`, for the backward cost-to-go
recursions at fixed gains. Terminal values are `(Q[:,:,T], 0, 0)`.
"""
function Compute_Lagrange_Multipliers(T, K_matrix, L_matrix, A, B, H, C, D, Q_matrix, R_matrix)

    n = size(A, 1)

    λ = zeros(n, n, T)
    ω = zeros(n, n, T)
    ν = zeros(n, n, T)

    M = similar(A)

    Atr = A'
    Btr = B'
    Htr = H'

    @views λ[:, :, T] .= Q_matrix[:, :, T]
    @views ω[:, :, T] .= 0
    @views ν[:, :, T] .= 0

    for t in (T-1):-1:1

        @views K_t = K_matrix[:, :, t]
        @views L_t = L_matrix[:, :, t]
        @views λ_next = λ[:, :, t+1]
        @views ω_next = ω[:, :, t+1]
        @views ν_next = ν[:, :, t+1]
        @views Q_t = Q_matrix[:, :, t]
        @views R_t = R_matrix[:, :, t]

        M .= A + B*L_t - K_t*H

        λ[:, :, t] .= Q_t + Atr*λ_next*A + Htr*K_t'*ω_next*K_t*H + Htr*K_t'*ν_next*A + sumDtGD(D, K_t'*ω_next*K_t)

        ω[:, :, t] .= L_t'*R_t*L_t + L_t'*Btr*λ_next*B*L_t + sumLtCtGCL(C, L_t, λ_next) + M'*ω_next*M + M'*ν_next*B*L_t

        Sλ_next = λ_next + λ_next'
        Sω_next = ω_next + ω_next'

        ν[:, :, t] .= L_t'*Btr*Sλ_next*A + M'*Sω_next*K_t*H + M'*ν_next*A + L_t'*Btr*ν_next'*K_t*H

    end

    return λ, ω, ν
end

#(C.1):Regularised inverses approach ----------------------------------------------------------------------
# Optimal Control and Filter gains at time t through regularised matrix inversions
function Optimal_L_t_Inverse!(t, K_t, L_t, ρ, A, B, H, C, D, Q_matrix, R_matrix, S_t_zz, S_t_xz, λ, ω, ν)

    Btr = B'

    @views λt_p_1 = λ[:, :, t+1]
    @views ωt_p_1 = ω[:, :, t+1]
    @views νt_p_1 = ν[:, :, t+1]
    @views Rt     = R_matrix[:, :, t]

    Sλ = λt_p_1 + λt_p_1'
    Sω = ωt_p_1 + ωt_p_1'
    Sν = νt_p_1 + νt_p_1'

    M_tmp = A - K_t*H

    J = 2*Rt + Btr*Sλ*B + sumCtGC(C, Sλ) + Btr*Sω*B + Btr*Sν*B
    N = Btr*Sλ*A + Btr*Sω*K_t*H + Btr*νt_p_1*A + Btr*νt_p_1'*K_t*H
    P = Btr*Sω*M_tmp + Btr*νt_p_1'*M_tmp

    Rhs = -(N*S_t_xz + P*S_t_zz)

    # Damped versions (for stability if J or S_t_zz is ill-conditioned)
    n = size(J, 1)
    m = size(S_t_zz, 1)
    J̃ = J + ρ * I(n)
    S̃ = S_t_zz + ρ * I(m)

    # X = J̃ \ Rhs #uncomment these lines to use direct solve
    # L_t .= (S̃' \ X')'

    #let's solve with pseudoinverse as backup
    L_t .= pinv(J̃) * Rhs * pinv(S̃) #this one is more stable and robust for small ρ (and works fine also for ρ=0)
    # L_t .= pinv(J) * Rhs * pinv(S_t_zz)

    return nothing
end

function Optimal_K_t_Inverse!(t, K_t, L_t, ρ, Ω_ω, A, B, H, C, D, Q_matrix, R_matrix, S_t_xx, S_t_xz, S_t_zx, S_t_zz, λ, ω, ν)

    @views ωt_p_1 = ω[:, :, t+1]
    @views νt_p_1 = ν[:, :, t+1]

    Sω = ωt_p_1 + ωt_p_1'

    Δxx_xz = S_t_xx - S_t_xz
    Δzz_zx = S_t_zz - S_t_zx
    Δzx_zz = S_t_zx - S_t_zz

    F = Sω
    G = H*Δxx_xz*H' + H*Δzz_zx*H' + sumDGDt(D, S_t_xx) + Ω_ω
    V = Sω*(A + B*L_t)*Δzx_zz*H' + νt_p_1*(A*Δxx_xz + B*L_t*Δzx_zz)*H'

    Rhs = -V

    n = size(F, 1)
    m = size(G, 1)

    F̃ = F + ρ * I(n)
    G̃ = G + ρ * I(m)

    # X = F̃ \ Rhs #uncomment these lines to use direct solve
    # K_t .= (G̃' \ X')'

    #let's solve with pseudoinverse as backup
    K_t .= pinv(F̃) * Rhs * pinv(G̃) #this one is more stable and robust for small ρ (and works fine also for ρ=0)
    # K_t .= pinv(F) * Rhs * pinv(G)

    return nothing
end

# Optimal Control and Filter gains for the whole trial
function Optimal_Control_Gains_with_Regularized_Inverse(K_matrix_old, L_initial, ρ, iters_L, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean)

    # Dimensions inferred from matrices
    n = size(A, 1)

    # Preallocate
    K_matrix = copy(K_matrix_old)
    L_matrix = copy(L_initial)

    # Moment buffers
    S_t_xx = zeros(n, n)
    S_t_zz = zeros(n, n)
    S_t_xz = zeros(n, n)
    S_t_zx = zeros(n, n)
    S_t_xx_old = similar(S_t_xx)
    S_t_zz_old = similar(S_t_zz)
    S_t_xz_old = similar(S_t_xz)
    S_t_zx_old = similar(S_t_zx)

    Atr = A'

    for _ in 1:iters_L

        λ, ω, ν = Compute_Lagrange_Multipliers(
            T, K_matrix, L_matrix, A, B, H, C, D, Q_matrix, R_matrix
        )

        # Initialize non-central moments
        S_t_xx_old .= Σ_1_x .+ x_1_mean * x_1_mean'
        S_t_zz_old .= Σ_1_z .+ z_1_mean * z_1_mean'
        S_t_xz_old .= x_1_mean * z_1_mean'
        S_t_zx_old .= S_t_xz_old'

        for t in 1:(T-1)

            @views L_t = L_matrix[:, :, t]
            @views K_t = K_matrix[:, :, t]

            Optimal_L_t_Inverse!(t, K_t, L_t, ρ, A, B, H, C, D, Q_matrix, R_matrix, S_t_zz_old, S_t_xz_old, λ, ω, ν)

            # Precompute recurrent combinations
            BL = B*L_t
            M  = A + BL - K_t*H
            KH = K_t*H
            BLt = BL'  # = L_t'*B'
            KHt = KH'
            Mt  = M'

            # Update moments
            S_t_xx .= A*S_t_xx_old*A' + A*S_t_xz_old*BLt + BL*S_t_zx_old*Atr + BL*S_t_zz_old*BLt + sumCLMLC(C, L_t, S_t_zz_old) + Ω_ξ
            S_t_zx .= KH*S_t_xx_old*Atr + KH*S_t_xz_old*BLt + M*S_t_zx_old*Atr + M*S_t_zz_old*BLt
            S_t_xz .= S_t_zx'
            S_t_zz .= KH*S_t_xx_old*KHt + sumKDMDK(D, K_t, S_t_xx_old) + KH*S_t_xz_old*Mt + M*S_t_zx_old*KHt + M*S_t_zz_old*Mt + K_t*Ω_ω*K_t' + Ω_η

            # Update old moments for next time step
            S_t_xx_old .= S_t_xx
            S_t_xz_old .= S_t_xz
            S_t_zx_old .= S_t_zx
            S_t_zz_old .= S_t_zz
        end
    end

    return L_matrix

end

function Optimal_Filter_Gains_with_Regularized_Inverse(K_initial, L_matrix_old, ρ, iters_K, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean)

    n = size(A, 1)

    # Preallocate arrays
    K_matrix = copy(K_initial)
    L_matrix = copy(L_matrix_old)

    S_t_xx = zeros(n, n)
    S_t_zz = zeros(n, n)
    S_t_xz = zeros(n, n)
    S_t_zx = zeros(n, n)
    S_t_xx_old = similar(S_t_xx)
    S_t_zz_old = similar(S_t_zz)
    S_t_xz_old = similar(S_t_xz)
    S_t_zx_old = similar(S_t_zx)

    Atr = A'

    for _ in 1:iters_K

        λ, ω, ν = Compute_Lagrange_Multipliers(
            T, K_matrix, L_matrix, A, B, H, C, D, Q_matrix, R_matrix
        )

        # Initialize non-central moments
        S_t_xx_old .= Σ_1_x .+ x_1_mean * x_1_mean'
        S_t_zz_old .= Σ_1_z .+ z_1_mean * z_1_mean'
        S_t_xz_old .= x_1_mean * z_1_mean'
        S_t_zx_old .= S_t_xz_old'

        for t in 1:(T-1)

            @views L_t = L_matrix[:, :, t]
            @views K_t = K_matrix[:, :, t]

            Optimal_K_t_Inverse!(t, K_t, L_t, ρ, Ω_ω, A, B, H, C, D, Q_matrix, R_matrix, S_t_xx_old, S_t_xz_old, S_t_zx_old, S_t_zz_old, λ, ω, ν)

            # Precompute recurrent combinations
            BL  = B*L_t
            M   = A + BL - K_t*H
            KH  = K_t*H
            BLt = BL'
            KHt = KH'
            Mt  = M'

            # Update moments
            S_t_xx .= A*S_t_xx_old*A' + A*S_t_xz_old*BLt + BL*S_t_zx_old*Atr + BL*S_t_zz_old*BLt + sumCLMLC(C, L_t, S_t_zz_old) + Ω_ξ
            S_t_zx .= KH*S_t_xx_old*Atr + KH*S_t_xz_old*BLt + M*S_t_zx_old*Atr + M*S_t_zz_old*BLt
            S_t_xz .= S_t_zx'
            S_t_zz .= KH*S_t_xx_old*KHt + sumKDMDK(D, K_t, S_t_xx_old) + KH*S_t_xz_old*Mt + M*S_t_zx_old*KHt + M*S_t_zz_old*Mt + K_t*Ω_ω*K_t' + Ω_η

            # Update old moments
            S_t_xx_old .= S_t_xx
            S_t_xz_old .= S_t_xz
            S_t_zx_old .= S_t_zx
            S_t_zz_old .= S_t_zz
        end
    end

    return K_matrix

end

# Whole optimization routine with coordinate descent using regularized inverses.
# The complete cost trajectory is retained for diagnostics, while the returned
# gains are the checkpoint with the lowest finite cost encountered.
function Optimal_Control_Estimation_with_Regularized_Inverse(K_initial, L_initial, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean; N_iter_tot = 20, ρ = 1e-10, iters_L=10, iters_K=10)

    n = size(A, 1)
    m = size(B, 2)
    p = size(H, 1)

    K_matrix = copy(K_initial)
    L_matrix = copy(L_initial)

    cost_mom_prop = zeros(2*N_iter_tot + 1)

    cost_mom_prop[1] = expected_cost_raw_moments_propagation(T, n, m, p, K_matrix, L_matrix, A, B, H, C, D, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η)

    best_cost = isfinite(cost_mom_prop[1]) ? cost_mom_prop[1] : Inf
    K_best = copy(K_matrix)
    L_best = copy(L_matrix)

    for i in 1:N_iter_tot

        L_matrix .= Optimal_Control_Gains_with_Regularized_Inverse(K_matrix, L_matrix, ρ, iters_L, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean)

        cost_mom_prop[2*i] = expected_cost_raw_moments_propagation(T, n, m, p, K_matrix, L_matrix, A, B, H, C, D, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η)

        if isfinite(cost_mom_prop[2*i]) && cost_mom_prop[2*i] < best_cost
            best_cost = cost_mom_prop[2*i]
            K_best .= K_matrix
            L_best .= L_matrix
        end

        K_matrix .= Optimal_Filter_Gains_with_Regularized_Inverse(K_matrix, L_matrix, ρ, iters_K, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean)

        cost_mom_prop[2*i+1] = expected_cost_raw_moments_propagation(T, n, m, p, K_matrix, L_matrix, A, B, H, C, D, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η)

        if isfinite(cost_mom_prop[2*i+1]) && cost_mom_prop[2*i+1] < best_cost
            best_cost = cost_mom_prop[2*i+1]
            K_best .= K_matrix
            L_best .= L_matrix
        end
    end

    return K_best, L_best, cost_mom_prop

end

#(C.2): Analytical GD approach ----------------------------------------------------------------------
function Gradient_for_L_t(t, K_t, L_t, A, B, H, C, D, Q_matrix, R_matrix, S_t_zz, S_t_xz, λ, ω, ν)

    Btr = B'

    @views λt_p_1 = λ[:, :, t+1]
    @views ωt_p_1 = ω[:, :, t+1]
    @views νt_p_1 = ν[:, :, t+1]
    @views Rt = R_matrix[:, :, t]

    Sλ = λt_p_1 + λt_p_1'
    Sω = ωt_p_1 + ωt_p_1'
    Sν = νt_p_1 + νt_p_1'

    M_tmp = A - K_t*H

    J = 2*Rt + Btr * Sλ * B + sumCtGC(C, Sλ) + Btr * Sω * B + Btr * Sν * B
    N = Btr * Sλ * A + Btr * Sω * K_t * H + Btr * νt_p_1 * A + Btr * νt_p_1' * K_t * H
    P = Btr * Sω * M_tmp + Btr * νt_p_1' * M_tmp

    return J * L_t * S_t_zz + N * S_t_xz + P * S_t_zz
end

function Gradient_for_K_t(t, K_t, L_t, Ω_ω, A, B, H, C, D, Q_matrix, R_matrix, S_t_xx, S_t_xz, S_t_zx, S_t_zz, λ, ω, ν)

    @views ωt_p_1 = ω[:, :, t+1]
    @views νt_p_1 = ν[:, :, t+1]

    Sω = ωt_p_1 + ωt_p_1'

    Δxx_xz = S_t_xx - S_t_xz
    Δzz_zx = S_t_zz - S_t_zx
    Δzx_zz = S_t_zx - S_t_zz

    F = Sω
    G = H*Δxx_xz*H' + H*Δzz_zx*H' + sumDGDt(D, S_t_xx) + Ω_ω
    V = Sω*(A + B*L_t)*Δzx_zz*H' + νt_p_1*(A*Δxx_xz + B*L_t*Δzx_zz)*H'

    return F*K_t*G + V
end

function gradient_descent_L(t, K_t, L_t, A, B, H, C, D, Q_matrix, R_matrix, S_t_zz, S_t_xz, λ, ω, ν; η=1e-3, maxiter=100, tol=1e-6)

    ∇L = similar(L_t)

    for _ in 1:maxiter
        ∇L .= Gradient_for_L_t(t, K_t, L_t, A, B, H, C, D, Q_matrix, R_matrix, S_t_zz, S_t_xz, λ, ω, ν)
        nrm = norm(∇L)
        nrm < tol && break
        L_t .-= η .* ∇L
    end

    return L_t
end

function gradient_descent_K(t, K_t, L_t, Ω_ω, A, B, H, C, D, Q_matrix, R_matrix, S_t_xx, S_t_xz, S_t_zx, S_t_zz, λ, ω, ν; η=1e-3, maxiter=100, tol=1e-6)

    ∇K = similar(K_t)

    for _ in 1:maxiter
        ∇K .= Gradient_for_K_t(t, K_t, L_t, Ω_ω, A, B, H, C, D, Q_matrix, R_matrix, S_t_xx, S_t_xz, S_t_zx, S_t_zz, λ, ω, ν)
        nrm = norm(∇K)
        nrm < tol && break
        K_t .-= η .* ∇K
    end

    return K_t
end

# ADAM optimizer
function adam_L(t, K_t, L_t, A, B, H, C, D, Q_matrix, R_matrix, S_t_zz, S_t_xz, λ, ω, ν; η=1e-2, β1=0.9, β2=0.999, ϵ=1e-8, maxiter=100, tol=1e-6)

    ∇L = similar(L_t)
    m = zeros(size(L_t))
    v = zeros(size(L_t))
    mhat = similar(L_t)
    vhat = similar(L_t)

    for k in 1:maxiter

        ∇L .= Gradient_for_L_t(t, K_t, L_t, A, B, H, C, D, Q_matrix, R_matrix, S_t_zz, S_t_xz, λ, ω, ν)

        gnorm = norm(∇L)
        gnorm < tol && break

        m .*= β1
        m .+= (1-β1) .* ∇L

        v .*= β2
        v .+= (1-β2) .* (∇L.^2)

        mhat .= m ./ (1 - β1^k)
        vhat .= v ./ (1 - β2^k)

        L_t .-= η .* mhat ./ (sqrt.(vhat) .+ ϵ)
    end

    return L_t

end

function adam_K(t, K_t, L_t, Ω_ω, A, B, H, C, D, Q_matrix, R_matrix, S_t_xx, S_t_xz, S_t_zx, S_t_zz, λ, ω, ν; η=1e-2, β1=0.9, β2=0.999, ϵ=1e-8, maxiter=100, tol=1e-6)

    ∇K = similar(K_t)
    m = zeros(size(K_t))
    v = zeros(size(K_t))
    mhat = similar(K_t)
    vhat = similar(K_t)

    for k in 1:maxiter

        ∇K .= Gradient_for_K_t(t, K_t, L_t, Ω_ω, A, B, H, C, D, Q_matrix, R_matrix, S_t_xx, S_t_xz, S_t_zx, S_t_zz, λ, ω, ν)

        gnorm = norm(∇K)
        gnorm < tol && break

        m .*= β1
        m .+= (1-β1) .* ∇K

        v .*= β2
        v .+= (1-β2) .* (∇K.^2)

        mhat .= m ./ (1 - β1^k)
        vhat .= v ./ (1 - β2^k)

        K_t .-= η .* mhat ./ (sqrt.(vhat) .+ ϵ)
    end

    return K_t

end

function Optimal_Control_Gains_with_Analytical_GD(K_matrix_old, L_initial, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean; η_GD = 1e-3, maxiter_GD = 100, tol_GD = 1e-6, iters_L = 10, β1_Adam = 0.9, β2_Adam = 0.999, ϵ_Adam = 1e-8)

    n = size(A, 1)

    K_matrix = copy(K_matrix_old)
    L_matrix = copy(L_initial)

    S_t_xx = zeros(n, n)
    S_t_zz = zeros(n, n)
    S_t_xz = zeros(n, n)
    S_t_zx = zeros(n, n)
    S_t_xx_old = similar(S_t_xx)
    S_t_zz_old = similar(S_t_zz)
    S_t_xz_old = similar(S_t_xz)
    S_t_zx_old = similar(S_t_zx)

    Atr = A'

    for _ in 1:iters_L

        λ, ω, ν = Compute_Lagrange_Multipliers(T, K_matrix, L_matrix, A, B, H, C, D, Q_matrix, R_matrix)

        # Initialize non-central moments
        S_t_xx_old .= Σ_1_x .+ x_1_mean * x_1_mean'
        S_t_zz_old .= Σ_1_z .+ z_1_mean * z_1_mean'
        S_t_xz_old .= x_1_mean * z_1_mean'
        S_t_zx_old .= S_t_xz_old'

        for t in 1:(T-1)

            @views L_t = L_matrix[:, :, t]
            @views K_t = K_matrix[:, :, t]

            adam_L(t, K_t, L_t, A, B, H, C, D, Q_matrix, R_matrix, S_t_zz_old, S_t_xz_old, λ, ω, ν; η = η_GD, maxiter = maxiter_GD, tol = tol_GD, β1 = β1_Adam, β2 = β2_Adam, ϵ = ϵ_Adam)

            # Precompute recurrent combinations
            BL  = B*L_t
            M   = A + BL - K_t*H
            KH  = K_t*H
            BLt = BL'
            KHt = KH'
            Mt  = M'

            # Update moments
            S_t_xx .= A*S_t_xx_old*A' + A*S_t_xz_old*BLt +
                      BL*S_t_zx_old*Atr   + BL*S_t_zz_old*BLt +
                      sumCLMLC(C, L_t, S_t_zz_old)   + Ω_ξ

            S_t_zx .= KH*S_t_xx_old*Atr + KH*S_t_xz_old*BLt +
                      M*S_t_zx_old*Atr   + M*S_t_zz_old*BLt

            S_t_xz .= S_t_zx'

            S_t_zz .= KH*S_t_xx_old*KHt + sumKDMDK(D, K_t, S_t_xx_old) +
                      KH*S_t_xz_old*Mt   + M*S_t_zx_old*KHt +
                      M*S_t_zz_old*Mt    + K_t*Ω_ω*K_t' + Ω_η

            # Update old moments
            S_t_xx_old .= S_t_xx
            S_t_xz_old .= S_t_xz
            S_t_zx_old .= S_t_zx
            S_t_zz_old .= S_t_zz
        end
    end

    return L_matrix
end

function Optimal_Filter_Gains_with_Analytical_GD(K_initial, L_matrix_old, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean; η_GD = 1e-3, maxiter_GD = 100, tol_GD = 1e-6, iters_K = 10, β1_Adam = 0.9, β2_Adam = 0.999, ϵ_Adam = 1e-8)

    n = size(A, 1)

    K_matrix = copy(K_initial)
    L_matrix = copy(L_matrix_old)

    S_t_xx = zeros(n, n)
    S_t_zz = zeros(n, n)
    S_t_xz = zeros(n, n)
    S_t_zx = zeros(n, n)
    S_t_xx_old = similar(S_t_xx)
    S_t_zz_old = similar(S_t_zz)
    S_t_xz_old = similar(S_t_xz)
    S_t_zx_old = similar(S_t_zx)

    Atr = A'

    for _ in 1:iters_K

        λ, ω, ν = Compute_Lagrange_Multipliers(T, K_matrix, L_matrix, A, B, H, C, D, Q_matrix, R_matrix)

        # Initialize non-central moments
        S_t_xx_old .= Σ_1_x .+ x_1_mean * x_1_mean'
        S_t_zz_old .= Σ_1_z .+ z_1_mean * z_1_mean'
        S_t_xz_old .= x_1_mean * z_1_mean'
        S_t_zx_old .= S_t_xz_old'

        for t in 1:(T-1)

            @views L_t = L_matrix[:, :, t]
            @views K_t = K_matrix[:, :, t]

            adam_K(t, K_t, L_t, Ω_ω, A, B, H, C, D, Q_matrix, R_matrix, S_t_xx_old, S_t_xz_old, S_t_zx_old, S_t_zz_old, λ, ω, ν; η = η_GD, maxiter = maxiter_GD, tol = tol_GD, β1 = β1_Adam, β2 = β2_Adam, ϵ = ϵ_Adam)

            # Precompute recurrent combinations
            BL  = B*L_t
            M   = A + BL - K_t*H
            KH  = K_t*H
            BLt = BL'
            KHt = KH'
            Mt  = M'

            # Update moments
            S_t_xx .= A*S_t_xx_old*A' + A*S_t_xz_old*BLt +
                      BL*S_t_zx_old*Atr   + BL*S_t_zz_old*BLt +
                      sumCLMLC(C, L_t, S_t_zz_old)   + Ω_ξ

            S_t_zx .= KH*S_t_xx_old*Atr + KH*S_t_xz_old*BLt +
                      M*S_t_zx_old*Atr   + M*S_t_zz_old*BLt

            S_t_xz .= S_t_zx'

            S_t_zz .= KH*S_t_xx_old*KHt + sumKDMDK(D, K_t, S_t_xx_old) +
                      KH*S_t_xz_old*Mt   + M*S_t_zx_old*KHt +
                      M*S_t_zz_old*Mt    + K_t*Ω_ω*K_t' + Ω_η

            # Update old moments
            S_t_xx_old .= S_t_xx
            S_t_xz_old .= S_t_xz
            S_t_zx_old .= S_t_zx
            S_t_zz_old .= S_t_zz
        end
    end

    return K_matrix
end

function Optimal_Control_Estimation_with_Analytical_GD(K_initial, L_initial, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean; N_iter_tot = 20, η_GD = 1e-3, maxiter_GD = 100, tol_GD = 1e-6, iters_K = 10, iters_L = 10, β1_Adam = 0.9, β2_Adam = 0.999, ϵ_Adam = 1e-8)

    n = size(A, 1)
    m = size(B, 2)
    p = size(H, 1)

    K_matrix = copy(K_initial)
    L_matrix = copy(L_initial)

    cost_mom_prop = zeros(N_iter_tot + 1)

    cost_mom_prop[1] = expected_cost_raw_moments_propagation(T, n, m, p, K_matrix, L_matrix, A, B, H, C, D, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η)

    for i in 1:N_iter_tot

        L_matrix .= Optimal_Control_Gains_with_Analytical_GD(K_matrix, L_matrix, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean; η_GD = η_GD, maxiter_GD = maxiter_GD, tol_GD = tol_GD, iters_L = iters_L, β1_Adam = β1_Adam, β2_Adam = β2_Adam, ϵ_Adam = ϵ_Adam)

        K_matrix .= Optimal_Filter_Gains_with_Analytical_GD(K_matrix, L_matrix, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean; η_GD = η_GD, maxiter_GD = maxiter_GD, tol_GD = tol_GD, iters_K = iters_K, β1_Adam = β1_Adam, β2_Adam = β2_Adam, ϵ_Adam = ϵ_Adam)

        cost_mom_prop[i+1] = expected_cost_raw_moments_propagation(T, n, m, p, K_matrix, L_matrix, A, B, H, C, D, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η)
    end

    return K_matrix, L_matrix, cost_mom_prop
end

# Standard GD
function Optimal_Control_Gains_with_Standard_GD(K_matrix_old, L_initial, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean; η_GD = 1e-3, maxiter_GD = 100, tol_GD = 1e-6, iters_L = 10)

    n = size(A, 1)

    K_matrix = copy(K_matrix_old)
    L_matrix = copy(L_initial)

    S_t_xx = zeros(n, n)
    S_t_zz = zeros(n, n)
    S_t_xz = zeros(n, n)
    S_t_zx = zeros(n, n)
    S_t_xx_old = similar(S_t_xx)
    S_t_zz_old = similar(S_t_zz)
    S_t_xz_old = similar(S_t_xz)
    S_t_zx_old = similar(S_t_zx)

    Atr = A'

    for _ in 1:iters_L

        λ, ω, ν = Compute_Lagrange_Multipliers(T, K_matrix, L_matrix, A, B, H, C, D, Q_matrix, R_matrix)

        # Initialize non-central moments
        S_t_xx_old .= Σ_1_x .+ x_1_mean * x_1_mean'
        S_t_zz_old .= Σ_1_z .+ z_1_mean * z_1_mean'
        S_t_xz_old .= x_1_mean * z_1_mean'
        S_t_zx_old .= S_t_xz_old'

        for t in 1:(T-1)

            @views L_t = L_matrix[:, :, t]
            @views K_t = K_matrix[:, :, t]

            # Standard GD on L_t
            gradient_descent_L(t, K_t, L_t, A, B, H, C, D, Q_matrix, R_matrix, S_t_zz_old, S_t_xz_old, λ, ω, ν; η = η_GD, maxiter = maxiter_GD, tol = tol_GD)

            # Precompute recurrent combinations
            BL  = B*L_t
            M   = A + BL - K_t*H
            KH  = K_t*H
            BLt = BL'
            KHt = KH'
            Mt  = M'

            # Update moments
            S_t_xx .= A*S_t_xx_old*A' + A*S_t_xz_old*BLt +
                      BL*S_t_zx_old*Atr   + BL*S_t_zz_old*BLt +
                      sumCLMLC(C, L_t, S_t_zz_old)  + Ω_ξ

            S_t_zx .= KH*S_t_xx_old*Atr + KH*S_t_xz_old*BLt +
                      M*S_t_zx_old*Atr   + M*S_t_zz_old*BLt

            S_t_xz .= S_t_zx'

            S_t_zz .= KH*S_t_xx_old*KHt + sumKDMDK(D, K_t, S_t_xx_old) +
                      KH*S_t_xz_old*Mt   + M*S_t_zx_old*KHt +
                      M*S_t_zz_old*Mt    + K_t*Ω_ω*K_t' + Ω_η

            # Update old moments
            S_t_xx_old .= S_t_xx
            S_t_xz_old .= S_t_xz
            S_t_zx_old .= S_t_zx
            S_t_zz_old .= S_t_zz
        end
    end

    return L_matrix
end

function Optimal_Filter_Gains_with_Standard_GD(K_initial, L_matrix_old, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean; η_GD = 1e-3, maxiter_GD = 100, tol_GD = 1e-6, iters_K = 10)

    n = size(A, 1)

    K_matrix = copy(K_initial)
    L_matrix = copy(L_matrix_old)

    S_t_xx = zeros(n, n)
    S_t_zz = zeros(n, n)
    S_t_xz = zeros(n, n)
    S_t_zx = zeros(n, n)
    S_t_xx_old = similar(S_t_xx)
    S_t_zz_old = similar(S_t_zz)
    S_t_xz_old = similar(S_t_xz)
    S_t_zx_old = similar(S_t_zx)

    Atr = A'

    for _ in 1:iters_K

        λ, ω, ν = Compute_Lagrange_Multipliers(T, K_matrix, L_matrix, A, B, H, C, D, Q_matrix, R_matrix)

        # Initialize non-central moments
        S_t_xx_old .= Σ_1_x .+ x_1_mean * x_1_mean'
        S_t_zz_old .= Σ_1_z .+ z_1_mean * z_1_mean'
        S_t_xz_old .= x_1_mean * z_1_mean'
        S_t_zx_old .= S_t_xz_old'

        for t in 1:(T-1)

            @views L_t = L_matrix[:, :, t]
            @views K_t = K_matrix[:, :, t]

            # Standard GD on K_t
            gradient_descent_K(t, K_t, L_t, Ω_ω, A, B, H, C, D, Q_matrix, R_matrix, S_t_xx_old, S_t_xz_old, S_t_zx_old, S_t_zz_old, λ, ω, ν; η = η_GD, maxiter = maxiter_GD, tol = tol_GD)

            # Precompute recurrent combinations
            BL  = B*L_t
            M   = A + BL - K_t*H
            KH  = K_t*H
            BLt = BL'
            KHt = KH'
            Mt  = M'

            # Update moments
            S_t_xx .= A*S_t_xx_old*A' + A*S_t_xz_old*BLt +
                      BL*S_t_zx_old*Atr   + BL*S_t_zz_old*BLt +
                      sumCLMLC(C, L_t, S_t_zz_old)   + Ω_ξ

            S_t_zx .= KH*S_t_xx_old*Atr + KH*S_t_xz_old*BLt +
                      M*S_t_zx_old*Atr   + M*S_t_zz_old*BLt

            S_t_xz .= S_t_zx'

            S_t_zz .= KH*S_t_xx_old*KHt + sumKDMDK(D, K_t, S_t_xx_old) +
                      KH*S_t_xz_old*Mt   + M*S_t_zx_old*KHt +
                      M*S_t_zz_old*Mt    + K_t*Ω_ω*K_t' + Ω_η

            # Update old moments
            S_t_xx_old .= S_t_xx
            S_t_xz_old .= S_t_xz
            S_t_zx_old .= S_t_zx
            S_t_zz_old .= S_t_zz
        end
    end

    return K_matrix
end

function Optimal_Control_Estimation_with_Standard_GD(K_initial, L_initial, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean; N_iter_tot = 20, η_GD = 1e-3, maxiter_GD = 100, tol_GD = 1e-6, iters_K = 10, iters_L = 10)

    n = size(A, 1)
    m = size(B, 2)
    p = size(H, 1)

    K_matrix = copy(K_initial)
    L_matrix = copy(L_initial)

    cost_mom_prop = zeros(N_iter_tot + 1)

    cost_mom_prop[1] = expected_cost_raw_moments_propagation(T, n, m, p, K_matrix, L_matrix, A, B, H, C, D, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η)

    for i in 1:N_iter_tot

        L_matrix .= Optimal_Control_Gains_with_Standard_GD(K_matrix, L_matrix, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean; η_GD = η_GD, maxiter_GD = maxiter_GD, tol_GD = tol_GD, iters_L = iters_L)

        K_matrix .= Optimal_Filter_Gains_with_Standard_GD(K_matrix, L_matrix, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean; η_GD = η_GD, maxiter_GD = maxiter_GD, tol_GD = tol_GD, iters_K = iters_K)

        cost_mom_prop[i+1] = expected_cost_raw_moments_propagation(T, n, m, p, K_matrix, L_matrix, A, B, H, C, D, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η)
    end

    return K_matrix, L_matrix, cost_mom_prop

end

#Helper for choosing the method
function run_analytical_GD(K0, L0, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean; N_iter_tot = 20, η_GD = 1e-3, maxiter_GD = 100, tol_GD = 1e-6, iters_K = 10, iters_L = 10, gd_method = :adam, β1_Adam = 0.9, β2_Adam = 0.999, ϵ_Adam = 1e-8)

    if gd_method == :adam
        return Optimal_Control_Estimation_with_Analytical_GD(K0, L0, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean; N_iter_tot = N_iter_tot, η_GD = η_GD, maxiter_GD = maxiter_GD, tol_GD = tol_GD, iters_K = iters_K, iters_L = iters_L, β1_Adam = β1_Adam, β2_Adam = β2_Adam, ϵ_Adam = ϵ_Adam)
    elseif gd_method == :gd
        return Optimal_Control_Estimation_with_Standard_GD(K0, L0, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean; N_iter_tot = N_iter_tot, η_GD = η_GD, maxiter_GD = maxiter_GD, tol_GD = tol_GD, iters_K = iters_K, iters_L = iters_L)
    else
        error("Unknown gd_method = $(gd_method). Use :adam or :gd.")
    end

end

#(C.3): Full solver ----------------------------------------------------------------------
"""
High-level routine for solving the joint optimal **control and estimation** problem
using a robust, multi-stage procedure that combines *regularized inverse* and
*analytical gradient descent (GD or Adam)* methods.

### Overview of the procedure
1. **Initialization:**
   Compute the baseline cost of the initial (K, L) pair.

2. **ρ-scan with regularized inverses:**
   Try all values in `ρ_list` (and, if necessary, `ρ_list_large`) in **parallel**,
   each acting as a Tikhonov-like regularization term to stabilize matrix inversions.
   The smallest ρ producing an admissible result (finite, nonnegative cost no greater than
   the initialization up to numerical tolerance) is selected automatically.

3. **Analytical GD/Adam refinement (optional):**
   Refine the best regularized solution via gradient-based optimization of K and L.

4. **GD/Adam from initial condition (optional):**
   Optionally reruns optimization from the *same initial (K, L)* as Step 1,
   to ensure the regularized path didn’t miss a better basin.

5. **Failure behavior:**
   If no explicitly enabled route produces a valid solution, throw an error;
   no implicit GD/Adam fallback is run.

When a valid solution is found, the procedure returns the selected gains,
their cost, and diagnostic information. There is no convergence-status output
or outer stopping tolerance: each route uses the configured iteration counts.
The ideal exact-coordinate LCO guarantee does not establish monotonicity for
positive regularization or GD/Adam. Inspect the complete recorded trajectory.

Candidate and final acceptance use `cost_init + ε`, where
`ε = 1e-10*max(abs(cost_init),1.0)`. Optional GD/Adam costs and histories must
also be finite, nonnegative, and within the existing trajectory cost ceiling.
Admissible GD results replace an existing solution only when their cost is lower.

---

### **Arguments**

- `K_init`, `L_init` — Initial estimation and control gain matrices, respectively.
- `A, B, H, C, D` — Problem definition.
- `T` — Number of time samples in the trial.
- `Ω_ξ, Ω_ω, Ω_η` — Covariance matrices of process, observation, and internal noise.
- `Q_matrix`, `R_matrix` — State and control cost matrices.
- `Σ_1_x`, `Σ_1_z` — Initial state and estimate covariances.
- `x_1_mean`, `z_1_mean` — Initial state and estimate means.

The initial cross-covariance convention follows the forward cost routine and
the TOD comparison; it is not a mathematical restriction of LCO. The manuscript's
Theorem 1 allows a general initial joint distribution. Theorem 2 separately
assumes initial orthogonality, `S_0^{xz}=S_0^{zz}` in the paper's zero-based
notation; this is not equivalent to zero initial x/z cross-covariance.

---

### **Keyword parameters**

#### ρ-scan parameters
- `ρ_list = [0.0, 1e-16, 1e-12]`
  Small regularization values tried first for inverse stabilization.
  Automatically parallelized across all entries.
- `ρ_list_large = [1e-10, 1e-8, 1e-6]`
  Larger fallback values used only if all entries in `ρ_list` fail.

#### Iteration controls
- `N_iter_tot = 1000`
  Number of outer alternations between control and estimation updates.
- `iters_K = 1`, `iters_L = 1`
  Inner iterations for refining K and L within each outer loop.
- `cost_ceiling_multiple = 1e4`
  Reject a ρ-run if any cost in its trajectory exceeds this multiple of the
  initialization-cost magnitude. NaN/Inf and materially negative trajectory
  values are also rejected.
- `N_iter_GD_refine = 100`
    Number of outer iterations for GD/Adam refinement of the best ρ-solution. [alternation between lagrange multipliers computation and GD optimization]
- `N_iter_GD_from_IC = 1000`
    Number of outer iterations for GD/Adam run from the initial (K, L). [alternation between lagrange multipliers computation and GD optimization]

#### Gradient-based optimization
- `gd_method = :adam`
  Choose between `:adam` (default, adaptive momentum) or `:gd`
  (plain gradient descent with fixed step size).
  Note that if using ':gd' the learning rate needs to be adjusted accordingly: if too high the method might diverge and produce NaNs.
- `η_GD = 1e-3`
  Learning rate for gradient updates.
- `maxiter_GD = 50`
  Maximum inner iterations for each gradient update step (for L_t and K_t).
- `tol_GD = 1e-6`
  Gradient norm threshold for early stopping.
- `β1_Adam = 0.9`, `β2_Adam = 0.999`, `ϵ_Adam = 1e-8`
  Standard Adam hyperparameters (momentum decay and numerical stability).

#### Behavioral flags
- `use_GD_refinement = false`
  After the best ρ-solution is found, refine it via GD/Adam.
  Disabled by default; set to `true` to enable.
- `use_GD_from_IC = false`
  Also run GD/Adam from the *initial (K, L)* even if a good ρ-solution exists.
  Useful to compare regularized vs. unregularized basins.
- `verbose = true`
  Print progress messages, warnings, and costs.

---

### **Returns**
A named tuple with:
- `K_best`, `L_best` — Optimal gain matrices found.
- `best_cost` — Expected cost of the returned gains (the lowest retained checkpoint for regularized-inverse runs).
- `best_method` — Which explicitly enabled route produced it (e.g., `"RegularizedInverse"`, `"GD_from_IC"`).
- `best_ρ` — The ρ value used (if any; `NaN` otherwise).
- `cost_init` — Baseline cost of the initial gains.
- `best_cost_list` — Complete cost trajectory of the selected run; it may continue beyond the returned best checkpoint.

---

### **Default behavior summary**
By default, the function:
1. Tries small ρ values in parallel and picks the smallest ρ with an admissible result.
2. Does not run GD/Adam unless its corresponding flag is explicitly enabled.
3. Throws an error if no enabled route produces an admissible solution within the initialization-relative tolerance.
4. Returns the gains from the smallest admissible ρ, the corresponding optimization cost trajectory (best_cost_list), and diagnostic details (best_method, best_ρ, cost_init).

    [It seems like that running the refinement is not strictly needed, and it might be enough to just run more iterations of the regularized inverse method, so the default is set to false for now.
    It is indeed more efficient this way. However, we should always check that the solutions come from a small enough ρ, otherwise the regularization might have biased the solution too much.
    Therefore, it's important to always read the Info printed during the optimization and in case run some refinement steps or increase the number of iterations of the regularized inverse method,
    to be sure we are close enough to the true optimum of the unregularized problem.]

"""
# NOTE:
# If no GD/Adam refinement is enabled (use_GD_refinement = false and use_GD_from_IC = false),
# the variable `best_cost_list` contains the full cost trajectory of the selected optimization [Note that the cost is computed after each K or L optimisation -- so the list is 2*N_iter_tot +1 long].
# The returned best K/L checkpoint may occur before the end of this trajectory.
# When refinement is enabled, `best_cost_list` instead stores only the cost history
# of the final optimization stage that actually improved the solution.

function Optimal_Control_Estimation_FullProcedure(K_init, L_init, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean;
    ρ_list = [0.0, 1e-16, 1e-12],
    ρ_list_large = [1e-10, 1e-8, 1e-6],
    N_iter_tot = 1000,
    iters_K = 1, iters_L = 1,
    cost_ceiling_multiple = 1e4,
    N_iter_GD_refine = 100,
    N_iter_GD_from_IC = 1000,
    # GD / Adam hyperparameters
    gd_method = :adam,                # :adam or :gd
    η_GD = 1e-3,
    maxiter_GD = 50,
    tol_GD = 1e-6,
    β1_Adam = 0.9,
    β2_Adam = 0.999,
    ϵ_Adam = 1e-8,
    # Flags
    use_GD_refinement = false,         # refine best ρ-solution with GD/Adam
    use_GD_from_IC = false,           # run GD/Adam from initial IC even if ρ works
    verbose = true)

    # ----------------------------------------------------------
    # Step 1 — Initialization
    # ----------------------------------------------------------

    n = size(A,1); m = size(B,2); p = size(H,1)

    cost_init = expected_cost_raw_moments_propagation(
        T, n, m, p, K_init, L_init, A, B, H, C, D,
        x_1_mean, z_1_mean, Σ_1_x, Σ_1_z,
        Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η
    )

    isfinite(cost_ceiling_multiple) && cost_ceiling_multiple > 1 ||
        error("cost_ceiling_multiple must be finite and greater than one.")
    initialization_tolerance = 1e-10 * max(abs(cost_init), 1.0)
    if !isfinite(cost_init) || cost_init < -initialization_tolerance
        error("The initialization has an inadmissible expected cost: $cost_init")
    end
    trajectory_cost_floor = -initialization_tolerance
    trajectory_cost_ceiling = cost_ceiling_multiple * max(abs(cost_init), 1.0)

    verbose && @info "Initial cost = $cost_init"

    best_cost   = Inf
    best_K      = deepcopy(K_init)
    best_L      = deepcopy(L_init)
    best_method = "none"
    best_ρ      = NaN
    best_cost_list = [cost_init]

    # ----------------------------------------------------------
    # Step 2.1 — Try all ρ values with regularized inverses
    # ----------------------------------------------------------

    function try_rho_list(ρ_list_local)
        ResultType = NamedTuple{(:ρ, :cost, :K, :L, :cost_list),
                                Tuple{Float64, Float64, typeof(K_init), typeof(L_init), Vector{Float64}}}

        nρ = length(ρ_list_local)
        results = Vector{ResultType}(undef, nρ)

        Threads.@threads for i in 1:nρ
            ρ = ρ_list_local[i]
            local K_try, L_try, cost_list, cost_final
            verbose && @info "Trying regularized inverse with ρ = $ρ"
            try
                K_try, L_try, cost_list = Optimal_Control_Estimation_with_Regularized_Inverse(
                    K_init, L_init, A, B, H, C, D, T,
                    Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix,
                    Σ_1_x, Σ_1_z, x_1_mean, z_1_mean;
                    ρ = ρ, N_iter_tot = N_iter_tot,
                    iters_L = iters_L, iters_K = iters_K,
                )
                # The solver returns the lowest-cost K/L checkpoint, whereas
                # cost_list remains the complete (possibly non-monotone) path.
                cost_final = expected_cost_raw_moments_propagation(
                    T, n, m, p, K_try, L_try, A, B, H, C, D,
                    x_1_mean, z_1_mean, Σ_1_x, Σ_1_z,
                    Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η
                )
                trajectory_is_admissible = all(
                    cost -> isfinite(cost) &&
                            cost >= trajectory_cost_floor &&
                            cost <= trajectory_cost_ceiling,
                    cost_list,
                )
                if !trajectory_is_admissible
                    verbose && @warn "ρ=$ρ rejected because its cost trajectory left the admissible interval [$trajectory_cost_floor, $trajectory_cost_ceiling]."
                    cost_final = Inf
                elseif !isfinite(cost_final) ||
                       cost_final < trajectory_cost_floor ||
                       cost_final > trajectory_cost_ceiling
                    cost_final = Inf
                end
                verbose && @info "ρ=$ρ finished with cost=$cost_final"
            catch e
                verbose && @warn "Regularized inverse failed for ρ=$ρ — $(e)"
                K_try = K_init
                L_try = L_init
                cost_list = [Inf]
                cost_final = Inf
            end
            results[i] = (ρ = ρ, cost = cost_final, K = K_try, L = L_try, cost_list = cost_list)
        end

        return results
    end

    results = try_rho_list(ρ_list)

    # ----------------------------------------------------------
    # Step 2.2 — Select the smallest ρ among admissible results
    # ----------------------------------------------------------

    cost_tolerance = 1e-10 * max(abs(cost_init), 1.0)
    admissible_result(r) = isfinite(r.cost) &&
                           r.cost >= -cost_tolerance &&
                           r.cost <= cost_init + cost_tolerance
    admissible_GD_result(cost, cost_list) = isfinite(cost) &&
                                          cost >= 0 &&
                                          cost <= cost_init + cost_tolerance &&
                                          all(c -> isfinite(c) && 0 <= c <= trajectory_cost_ceiling,
                                              cost_list)
    valid_results = filter(admissible_result, results)

    if isempty(valid_results)
        verbose && @warn "No valid solutions in first ρ_list; trying larger ρ values..."
        results = try_rho_list(ρ_list_large)
        valid_results = filter(admissible_result, results)
    end

    if !isempty(valid_results)
        selected_idx = findmin(r -> r.ρ, valid_results)[2]
        best = valid_results[selected_idx]
        best_cost   = best.cost
        best_K      = best.K
        best_L      = best.L
        best_ρ      = best.ρ
        best_cost_list = best.cost_list
        best_method = "RegularizedInverse"
        verbose && @info "Selected smallest admissible ρ = $(best_ρ) with cost = $(best_cost)"
    else
        if use_GD_from_IC
            verbose && @warn "All regularized inverse runs failed or increased cost; the explicitly enabled GD/Adam-from-IC route will still be tried."
        else
            verbose && @warn "All regularized inverse runs failed or increased cost; no implicit GD/Adam fallback will be run."
        end
    end

    # ----------------------------------------------------------
    # Step 2.4 — Optional analytical GD refinement
    # ----------------------------------------------------------

    if use_GD_refinement && isfinite(best_cost)
        verbose && @info "Running analytical $(gd_method) refinement from best ρ-solution..."
        K_GD, L_GD, cost_list_GD = run_analytical_GD(
            best_K, best_L, A, B, H, C, D, T,
            Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix,
            Σ_1_x, Σ_1_z, x_1_mean, z_1_mean;
            N_iter_tot = N_iter_GD_refine,
            η_GD = η_GD, maxiter_GD = maxiter_GD, tol_GD = tol_GD,
            iters_K = iters_K, iters_L = iters_L,
            gd_method = gd_method,
            β1_Adam = β1_Adam, β2_Adam = β2_Adam, ϵ_Adam = ϵ_Adam,
        )
        cost_GD_final = isempty(cost_list_GD) ? Inf : cost_list_GD[end]
        if admissible_GD_result(cost_GD_final, cost_list_GD) && cost_GD_final < best_cost
            verbose && @info "GD/Adam refinement improved cost from $best_cost to $cost_GD_final"
            best_K          = K_GD
            best_L          = L_GD
            best_cost       = cost_GD_final
            best_method     = "RegularizedInverse+GD"
            best_cost_list  = cost_list_GD
        else
            verbose && @info "GD/Adam refinement was inadmissible or did not improve cost (candidate=$cost_GD_final, incumbent=$best_cost)"
        end
    end

    # ----------------------------------------------------------
    # Step 2.5 — Analytical GD from initial condition
    # ----------------------------------------------------------

    if use_GD_from_IC
        verbose && @info "Running analytical $(gd_method) from initial condition (IC)..."
        K_GD_IC, L_GD_IC, cost_list_GD_IC = run_analytical_GD(
            K_init, L_init, A, B, H, C, D, T,
            Ω_ξ, Ω_ω, Ω_η, Q_matrix, R_matrix,
            Σ_1_x, Σ_1_z, x_1_mean, z_1_mean;
            N_iter_tot = N_iter_GD_from_IC,
            η_GD = η_GD, maxiter_GD = maxiter_GD, tol_GD = tol_GD,
            iters_K = iters_K, iters_L = iters_L,
            gd_method = gd_method,
            β1_Adam = β1_Adam, β2_Adam = β2_Adam, ϵ_Adam = ϵ_Adam,
        )
        cost_GD_IC_final = isempty(cost_list_GD_IC) ? Inf : cost_list_GD_IC[end]
        # Initialization acceptance uses ε above; compare quality with the actual incumbent.
        incumbent_cost = best_cost
        if admissible_GD_result(cost_GD_IC_final, cost_list_GD_IC) && cost_GD_IC_final < incumbent_cost
            verbose && @info "GD/Adam from IC improved the incumbent cost from $incumbent_cost to $cost_GD_IC_final"
            best_K          = K_GD_IC
            best_L          = L_GD_IC
            best_cost       = cost_GD_IC_final
            best_method     = "GD_from_IC"
            best_ρ          = NaN
            best_cost_list  = cost_list_GD_IC
        else
            verbose && @info "GD/Adam from IC was inadmissible or did not beat the incumbent cost (candidate=$cost_GD_IC_final, incumbent=$incumbent_cost)"
        end
    end

    # ----------------------------------------------------------
    # Step 3 — Require a valid result from an explicitly enabled route
    # ----------------------------------------------------------

    if !isfinite(best_cost) || best_cost > cost_init + cost_tolerance
        error(
            "No valid solution was produced by the enabled optimization routes. " *
            "No regularized-inverse or explicitly enabled GD/Adam run produced " *
            "an admissible cost within the initialization-relative tolerance. " *
            "No implicit GD/Adam fallback was executed."
        )
    end

    # ----------------------------------------------------------
    # Step 4 — Return final solution
    # ----------------------------------------------------------

    return (
        K_best    = best_K,
        L_best    = best_L,
        best_cost = best_cost,
        best_method = best_method,
        best_ρ    = best_ρ,
        cost_init = cost_init,
        best_cost_list = best_cost_list
    )
end

#(D): Monte Carlo simulation ----------------------------------------------------------------------
# Run N_trials Monte Carlo rollouts of the closed-loop LQMI system (Eqs. 5-7 of the paper),
# given fixed control gains L_matrix (u_t = L_t z_t convention, i.e. the SAME sign used by the
# moment-propagation and LCO functions) and filter gains K_matrix.
#
# Returns per-trial traces for the true state x, estimate z, control u, observation y, together
# with the empirical (Monte Carlo) cost and its standard error. Multiplicative noise supports
# multiple C_i and D_j: pass C and D either as a single Matrix (legacy c=d=1) or as a Vector of
# matrices [C1,C2,...]/[D1,D2,...]; Σε = I, Σρ = I (independent unit-variance channels), matching
# the analytic functions.
#
# Trace array conventions (to match the [:,:,idx] gain indexing used elsewhere):
#   x_traces[:, t, trial]  for t = 1..T      (states x_1 .. x_T)
#   z_traces[:, t, trial]  for t = 1..T      (estimates z_1 .. z_T)
#   u_traces[:, t, trial]  for t = 1..T-1    (controls  u_1 .. u_{T-1}, u_t = L_t z_t)
#   y_traces[:, t, trial]  for t = 1..T-1    (observations y_1 .. y_{T-1}, generated from x_t)
#
# The empirical cost uses the same objective as the analytic routines:
#   C = Σ_{t=1}^{T-1} (x_t'Q_t*x_t + u_t'R_t*u_t)
#       + x_T'Q_T*x_T.
#
# NOTE: this MC estimate should converge (as N_trials -> ∞) to the analytic cost returned by
# expected_cost_raw_moments_propagation / expected_cost_using_trace_formula_induction, providing
# an independent numerical check of the analytic formulas.

"""
    monte_carlo_simulation(L, K, A, B, H, C, D, T, m, p, k,
        x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q, R, Ω_ξ, Ω_ω, Ω_η, N_trials;
        rng=Random.GLOBAL_RNG)

Simulate the LQMI dynamics with independent Gaussian draws of all five noise
sources. Return `x_traces`, `z_traces`, `u_traces`, `y_traces`, `cost_per_trial`,
`empirical_cost`, and `empirical_cost_SE`. Traces use `(dimension,time,trial)`.
Initial x/z draws are independent, following the zero-cross-covariance convention
used for the TOD comparison. This is not a general restriction of LCO; correlated
initialization would require joint initial-state sampling and the corresponding
initial cross-moment in the analytical routines. C/D channels have independent unit variance.
The analytical costs also apply to non-Gaussian noise with the same moments.
"""
function monte_carlo_simulation(L_matrix, K_matrix, A, B, H, C, D, T,
                                dimension_of_state, dimension_of_control, dimension_of_observation,
                                x_1_mean, z_1_mean, Σ_1_x, Σ_1_z,
                                Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η,
                                N_trials; rng = Random.GLOBAL_RNG)

    # multiplicative-noise matrices as lists (single Matrix -> 1-tuple, Vector -> itself)
    Clist = _aslist(C)
    Dlist = _aslist(D)
    c = length(Clist)   # number of control-noise channels   (Σε = I_c)
    d = length(Dlist)   # number of observation-noise channels (Σρ = I_d)

    # Cholesky factors of the additive-noise covariances, so that
    # L_fac * randn gives a sample with the desired covariance. We symmetrise;
    # singular covariances use an eigen-based factor with negative eigenvalues clipped.
    chol_factor(Σ) = begin
        Σs = Symmetric(0.5 .* (Σ .+ transpose(Σ)))
        try
            return Matrix(cholesky(Σs).L)
        catch
            # fall back to eigen-based square root for PSD (possibly singular) covariances
            vals, vecs = eigen(Σs)
            vals_clipped = map(v -> v > 0 ? v : 0.0, vals)
            return vecs * Diagonal(sqrt.(vals_clipped))
        end
    end

    Lξ  = chol_factor(Ω_ξ)                 # process noise ξ_t  ∈ R^m
    Lω  = chol_factor(Ω_ω)                 # obs noise ω_t      ∈ R^k
    Lη  = chol_factor(Ω_η)                 # internal noise η_t ∈ R^m
    Lx0 = chol_factor(Σ_1_x)              # initial true-state spread
    Lz0 = chol_factor(Σ_1_z)              # initial estimate spread

    # storage
    x_traces = zeros(dimension_of_state,       T,     N_trials)
    z_traces = zeros(dimension_of_state,       T,     N_trials)
    u_traces = zeros(dimension_of_control,     T - 1, N_trials)
    y_traces = zeros(dimension_of_observation, T - 1, N_trials)
    cost_per_trial = zeros(N_trials)

    for trial in 1:N_trials

        # ----- initial conditions (t = 1) -----
        x = x_1_mean .+ Lx0 * randn(rng, dimension_of_state)
        z = z_1_mean .+ Lz0 * randn(rng, dimension_of_state)
        x_traces[:, 1, trial] .= x
        z_traces[:, 1, trial] .= z

        trial_cost = 0.0

        for t in 1:(T - 1)

            # control  u_t = L_t z_t
            u = L_matrix[:, :, t] * z
            u_traces[:, t, trial] .= u

            # Instantaneous state and control cost at t.
            trial_cost += dot(x, Q_matrix[:, :, t] * x) + dot(u, R_matrix[:, :, t] * u)

            # observation of the CURRENT state x_t :  y_t = H x_t + ω_t + Σ_j ρ^j_t D_j x_t
            y = H * x .+ Lω * randn(rng, dimension_of_observation)
            for j in 1:d
                y .+= randn(rng) .* (Dlist[j] * x)
            end
            y_traces[:, t, trial] .= y

            # state update : x_{t+1} = A x_t + B u_t + ξ_t + Σ_i ε^i_t C_i u_t
            x_next = A * x .+ B * u .+ Lξ * randn(rng, dimension_of_state)
            for i in 1:c
                x_next .+= randn(rng) .* (Clist[i] * u)
            end

            # estimate update : z_{t+1} = A z_t + B u_t + K_t (y_t - H z_t) + η_t
            z_next = A * z .+ B * u .+ K_matrix[:, :, t] * (y .- H * z) .+ Lη * randn(rng, dimension_of_state)

            x = x_next
            z = z_next
            x_traces[:, t + 1, trial] .= x
            z_traces[:, t + 1, trial] .= z
        end

        # terminal cost : x_T' Q_T x_T
        trial_cost += dot(x, Q_matrix[:, :, T] * x)
        cost_per_trial[trial] = trial_cost
    end

    empirical_cost      = mean(cost_per_trial)
    empirical_cost_SE   = std(cost_per_trial) / sqrt(N_trials)

    return (
        x_traces          = x_traces,
        z_traces          = z_traces,
        u_traces          = u_traces,
        y_traces          = y_traces,
        cost_per_trial    = cost_per_trial,
        empirical_cost    = empirical_cost,
        empirical_cost_SE = empirical_cost_SE,
    )
end

#(E): Numerical utilities ----------------------------------------------------------------------

# ============================================================================
# Helpers for multiplicative-noise sums over multiple C_i and D_j matrices.
# C and D may be passed EITHER as a single Matrix (legacy, c=d=1) OR as a
# Vector of Matrices [C1, C2, ...] / [D1, D2, ...]. The `_aslist` wrapper makes
# both forms work, so single-element usage reproduces the original results
# exactly (Σε = I, Σρ = I : independent unit-variance noise per channel).
# ============================================================================
_aslist(M::AbstractMatrix) = (M,)          # legacy single matrix -> 1-tuple
_aslist(M) = M                              # already a vector/tuple of matrices

# Σ_i (Ci*L) * Mid * (Ci*L)'      (control-noise term in the x-x moment / Σxx)
sumCLMLC(C, L, Mid) = sum((Ci*L)*Mid*(Ci*L)' for Ci in _aslist(C))
# Σ_i (Ci*L)' * G * (Ci*L)        (control-noise term in backward recursions)
sumLtCtGCL(C, L, G) = sum((Ci*L)'*G*(Ci*L) for Ci in _aslist(C))
# Σ_i (K*Di) * Mid * (K*Di)'      (observation-noise term in the z-z moment / Σzz)
sumKDMDK(D, K, Mid) = sum((K*Di)*Mid*(K*Di)' for Di in _aslist(D))
# Σ_i Ci' * G * Ci                (e.g. Todorov controller term, LCO J term)
sumCtGC(C, G) = sum(Ci'*G*Ci for Ci in _aslist(C))
# Σ_i Di' * G * Di                (e.g. LCO λ-multiplier term)
sumDtGD(D, G) = sum(Di'*G*Di for Di in _aslist(D))
# Σ_i Di * G * Di'                (e.g. Kalman innovation / observation covariance)
sumDGDt(D, G) = sum(Di*G*Di' for Di in _aslist(D))
# ============================================================================
