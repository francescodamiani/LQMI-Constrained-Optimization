# Calibrated elbow-reaching model for Figure 2.
using LinearAlgebra, Random, Statistics

function reaching_problem(; Δt=0.01, T_move=0.40, σ_ξ=0.01)
N_transitions = round(Int, T_move/Δt)
T = N_transitions + 1

dimension_of_state = 6
dimension_of_control = 9
dimension_of_observation = 3

x_1_mean = [0.0, -pi/9, 0.0, 0.0, 0.0, 0.0]
z_1_mean = copy(x_1_mean)
x_target = zeros(dimension_of_state)

Σ_1_x = zeros(dimension_of_state, dimension_of_state)
Σ_1_z = zeros(dimension_of_state, dimension_of_state)

Q_run = Diagonal([
    1000.0, 1e-5, 1000.0,
    100.0, 1e-5, 100.0
])

Q_terminal = Diagonal([
    1000.0, 5000.0, 1000.0,
    100.0, 100.0, 100.0
])

r = 5e-3

σ_ω_add_sensory_noise = deg2rad(0.9)    #Proprioceptive precision: van Beers et al. (1998).
σ_ϵ_mult_control_noise = 0.10           #Signal-dependent motor noise: Liu & Todorov (2007).
σ_ρ_mult_sensory_noise = 0.01           #Multiplicative sensory-noise amplitude.

#DEFINITION OF THE ARM MODEL -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Same arm dynamics and actuator map as in the stabilization task
M_joint = Diagonal([1.2, 0.8, 0.5])
D_joint = Diagonal(2.0 .* ones(3))

S = [
    1.2  -1.0   0.0    0.8  -0.6   0.0    0.5   0.0   0.0;
    0.0   0.0   1.0   -0.4   0.6  -0.5    0.0   0.5   0.0;
    0.0   0.0   0.0    0.0   0.0   1.0    0.0  -0.3   0.6
]

A = zeros(dimension_of_state, dimension_of_state)
B = zeros(dimension_of_state, dimension_of_control)
H = zeros(dimension_of_observation, dimension_of_state)
D = zeros(dimension_of_observation, dimension_of_state)

M_inv = inv(M_joint)

A[1:3,1:3] .= I(3)
A[1:3,4:6] .= Δt .* I(3)
A[4:6,1:3] .= 0
A[4:6,4:6] .= I(3) .- Δt .* (M_inv*D_joint)

B[4:6,:] .= Δt .* (M_inv*S)
C = σ_ϵ_mult_control_noise .* B

H[:,1:3] .= I(3)
D .= σ_ρ_mult_sensory_noise .* H

Ω_ξ = Diagonal(fill(σ_ξ^2, dimension_of_state))
Ω_ω = Diagonal(fill(σ_ω_add_sensory_noise^2, dimension_of_observation))

#COST MATRICES -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Q_matrix = zeros(dimension_of_state, dimension_of_state, T)
R_matrix = zeros(dimension_of_control, dimension_of_control, T)

for t in 1:T-1
    Q_matrix[:,:,t] .= Q_run
    R_matrix[:,:,t] .= Diagonal(fill(r, dimension_of_control))
end

Q_matrix[:,:,T] .= Q_terminal
R_matrix[:,:,T] .= Diagonal(fill(0, dimension_of_control))


    return (;Δt, T_move, N_transitions, T, dimension_of_state, dimension_of_control,
        dimension_of_observation, x_1_mean, z_1_mean, x_target, Σ_1_x, Σ_1_z,
        Q_run, Q_terminal, r, σ_ξ, σ_ω_add_sensory_noise,
        σ_ϵ_mult_control_noise, σ_ρ_mult_sensory_noise, A, B, H, C, D, S,
        Ω_ξ, Ω_ω, Q_matrix, R_matrix)
end

function reaching_cost(p, L, K, σ_η)
    return expected_cost_raw_moments_propagation(p.T, p.dimension_of_state,
        p.dimension_of_control, p.dimension_of_observation, K, L, p.A, p.B,
        p.H, p.C, p.D, p.x_1_mean, p.z_1_mean, p.Σ_1_x, p.Σ_1_z,
        p.Q_matrix, p.R_matrix, p.Ω_ξ, p.Ω_ω,
        Matrix(Diagonal(fill(σ_η^2,p.dimension_of_state))))
end

function solve_reaching(p, σ_η; N_iter_Todorov=500, N_iter_tot=500)
    Ω_η = Matrix(Diagonal(fill(σ_η^2,p.dimension_of_state)))
    L_TOD, K_TOD = Todorov_optimization_multidimensional_case(
        p.A, p.B, p.H, p.C, p.D, p.T, p.dimension_of_state,
        p.dimension_of_control, p.dimension_of_observation, p.x_1_mean,
        p.Σ_1_x, p.Σ_1_z, p.Q_matrix, p.R_matrix, p.Ω_ξ, p.Ω_ω, Ω_η,
        N_iter_Todorov)
    result = Optimal_Control_Estimation_FullProcedure(
        copy(K_TOD), copy(L_TOD), p.A, p.B, p.H, p.C, p.D, p.T,
        p.Ω_ξ, p.Ω_ω, Ω_η, p.Q_matrix, p.R_matrix,
        p.Σ_1_x, p.Σ_1_z, p.x_1_mean, p.z_1_mean;
        ρ_list=[0.0], ρ_list_large=[1e-10,1e-8,1e-6],
        N_iter_tot=N_iter_tot, iters_K=1, iters_L=1,
        N_iter_GD_refine=0, N_iter_GD_from_IC=0,
        gd_method=:adam, η_GD=1e-3, maxiter_GD=50, tol_GD=1e-6,
        β1_Adam=0.9, β2_Adam=0.999, ϵ_Adam=1e-8,
        use_GD_refinement=false, use_GD_from_IC=false, verbose=false)
    return (;L_TOD, K_TOD, L_LCO=result.L_best, K_LCO=result.K_best,
        cost_history_LCO=collect(result.best_cost_list), method_LCO=result.best_method,
        ρ_LCO=result.best_ρ, cost_init_LCO=result.cost_init)
end

function simulate_reaching(p, L, K, σ_η, N_trials_MC; random_seed=1)
    return monte_carlo_simulation(L,K,p.A,p.B,p.H,p.C,p.D,p.T,
        p.dimension_of_state,p.dimension_of_control,p.dimension_of_observation,
        p.x_1_mean,p.z_1_mean,p.Σ_1_x,p.Σ_1_z,p.Q_matrix,p.R_matrix,
        Matrix(p.Ω_ξ),Matrix(p.Ω_ω),Matrix(Diagonal(fill(σ_η^2,p.dimension_of_state))),
        N_trials_MC; rng=MersenneTwister(random_seed))
end
