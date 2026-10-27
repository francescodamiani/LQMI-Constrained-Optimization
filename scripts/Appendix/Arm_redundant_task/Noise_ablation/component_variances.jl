# Exact deterministic moments for fixed-policy noise ablation.
# This is an analysis helper; stochastic trajectories use the validated core simulator.
using LinearAlgebra, Statistics

function compute_analytical_component_variances(
    K_matrix, L_matrix, A, B, H, C, D, T, temporal_window,
    x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Ω_ξ, Ω_ω, Ω_η
)

    dimension_of_state = size(A,1)

    S_t_xx = zeros(dimension_of_state, dimension_of_state, T)
    S_t_zz = zeros(dimension_of_state, dimension_of_state, T)
    S_t_xz = zeros(dimension_of_state, dimension_of_state, T)
    S_t_zx = zeros(dimension_of_state, dimension_of_state, T)

    x_mean = zeros(dimension_of_state, T)
    z_mean = zeros(dimension_of_state, T)

    x_mean[:,1] .= x_1_mean
    z_mean[:,1] .= z_1_mean

    S_t_xx[:,:,1] .=
        Σ_1_x .+ x_1_mean*transpose(x_1_mean)

    S_t_zz[:,:,1] .=
        Σ_1_z .+ z_1_mean*transpose(z_1_mean)

    S_t_xz[:,:,1] .=
        x_1_mean*transpose(z_1_mean)

    S_t_zx[:,:,1] .=
        transpose(S_t_xz[:,:,1])

    for t in 2:T

        L_t = L_matrix[:,:,t-1]
        K_t = K_matrix[:,:,t-1]
        M_t = A .+ B*L_t .- K_t*H

        x_mean[:,t] .=
            A*x_mean[:,t-1] .+
            B*L_t*z_mean[:,t-1]

        z_mean[:,t] .=
            K_t*H*x_mean[:,t-1] .+
            M_t*z_mean[:,t-1]

        S_t_xx[:,:,t] .=
            A*S_t_xx[:,:,t-1]*transpose(A) .+
            A*S_t_xz[:,:,t-1]*transpose(B*L_t) .+
            B*L_t*S_t_zx[:,:,t-1]*transpose(A) .+
            B*L_t*S_t_zz[:,:,t-1]*transpose(B*L_t) .+
            C*L_t*S_t_zz[:,:,t-1]*transpose(C*L_t) .+
            Ω_ξ

        S_t_zx[:,:,t] .=
            K_t*H*S_t_xx[:,:,t-1]*transpose(A) .+
            K_t*H*S_t_xz[:,:,t-1]*transpose(B*L_t) .+
            M_t*S_t_zx[:,:,t-1]*transpose(A) .+
            M_t*S_t_zz[:,:,t-1]*transpose(B*L_t)

        S_t_xz[:,:,t] .=
            transpose(S_t_zx[:,:,t])

        S_t_zz[:,:,t] .=
            K_t*H*S_t_xx[:,:,t-1]*transpose(K_t*H) .+
            K_t*D*S_t_xx[:,:,t-1]*transpose(K_t*D) .+
            K_t*H*S_t_xz[:,:,t-1]*transpose(M_t) .+
            M_t*S_t_zx[:,:,t-1]*transpose(K_t*H) .+
            M_t*S_t_zz[:,:,t-1]*transpose(M_t) .+
            K_t*Ω_ω*transpose(K_t) .+
            Ω_η

    end

    variance_x_over_time =
        zeros(dimension_of_state, T)

    variance_z_over_time =
        zeros(dimension_of_state, T)

    for t in 1:T
        variance_x_over_time[:,t] .=
            diag(S_t_xx[:,:,t]) .- x_mean[:,t].^2

        variance_z_over_time[:,t] .=
            diag(S_t_zz[:,:,t]) .- z_mean[:,t].^2
    end

    variance_x =
        vec(mean(
            variance_x_over_time[:,temporal_window],
            dims=2
        ))

    variance_z =
        vec(mean(
            variance_z_over_time[:,temporal_window],
            dims=2
        ))

    return variance_x, variance_z

end

