# Historical 3-DOF stabilization task; not the calibrated camera-ready reaching model.
using LinearAlgebra, Random, Statistics, Distributions

function arm_problem(; T=300)
dimension_of_state = 6
dimension_of_control = 9
dimension_of_observation = 3

Δt = 0.010


# Dynamics: inertia and damping
M_joint = Diagonal([1.2, 0.8, 0.5])


D_joint = Diagonal(2.0 .* ones(3))

# Muscle-to-joint map S
S = [
    1.2  -1.0   0.0    0.8  -0.6   0.0    0.5   0.0   0.0;
    0.0   0.0   1.0   -0.4   0.6  -0.5    0.0   0.5   0.0;
    0.0   0.0   0.0    0.0   0.0   1.0    0.0  -0.3   0.6
]

# Control cost weight
r = 1e-2

# Noise parameters
σ_ξ = 0.1
σ_ω_add_sensory_noise = 0.1
σ_ρ_mult_sensory_noise = 0.1
σ_ϵ_mult_control_noise = 0.1


# Initial conditions
x_state_mean_0 = zeros(dimension_of_state)
x_state_estimate_mean_0 = zeros(dimension_of_state)

Σ_1_x = zeros(dimension_of_state, dimension_of_state)
Σ_1_z = copy(Σ_1_x)

x_1_mean = deepcopy(x_state_mean_0)
z_1_mean = deepcopy(x_state_estimate_mean_0)

A = zeros(dimension_of_state, dimension_of_state)
B = zeros(dimension_of_state, dimension_of_control)
H = zeros(dimension_of_observation, dimension_of_state)
D = zeros(dimension_of_observation, dimension_of_state)

M_inv = inv(M_joint)

# Euler discretization
A[1:3,1:3] .= I(3)
A[1:3,4:6] .= Δt .* I(3)
A[4:6,1:3] .= 0
A[4:6,4:6] .= I(3) .- Δt .* (M_inv * D_joint)

B[4:6,:] .= Δt .* (M_inv * S)
C = σ_ϵ_mult_control_noise .* B

H[:,1:3] .= I(3)
D .= σ_ρ_mult_sensory_noise .* H

Ω_ξ = Diagonal(fill(σ_ξ^2, dimension_of_state))
Ω_ω = Diagonal(fill(σ_ω_add_sensory_noise^2, dimension_of_observation))

Q_matrix = zeros(dimension_of_state, dimension_of_state, T)
R_matrix = zeros(dimension_of_control, dimension_of_control, T)

for i in 1:T
    R_matrix[:,:,i] .= Diagonal(fill(r, dimension_of_control))
end

R_matrix[:,:,T] .= Diagonal(fill(0, dimension_of_control))

q_angles_intermediate = 1.0
q_vels_intermediate   = 1e-3

Q_intermediate = Diagonal([
    q_angles_intermediate, q_angles_intermediate, q_angles_intermediate,
    q_vels_intermediate,   q_vels_intermediate,   q_vels_intermediate
])

for i in 1:T-1
    Q_matrix[:,:,i] .= Q_intermediate
end

q_angles_final = 1.0
q_vels_final   = 1e-3
Q_matrix[:,:,T] .= Diagonal([
    q_angles_final, q_angles_final, q_angles_final,
    q_vels_final,   q_vels_final,   q_vels_final
])

    return (; dimension_of_state, dimension_of_control, dimension_of_observation, Δt, T, M_joint, D_joint, S, r, σ_ξ, σ_ω_add_sensory_noise, σ_ρ_mult_sensory_noise, σ_ϵ_mult_control_noise, x_state_mean_0, x_state_estimate_mean_0, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, A, B, H, C, D, Ω_ξ, Ω_ω, Q_matrix, R_matrix)
end

function solve_arm(p, σ_η; N_iter_Todorov=500, N_iter_tot=500, verbose=true)
    (; A, B, H, C, D, T, dimension_of_state, dimension_of_control,
       dimension_of_observation, x_state_mean_0, x_1_mean, z_1_mean,
       Σ_1_x, Σ_1_z, Q_matrix, R_matrix, Ω_ξ, Ω_ω) = p
    Ω_η = Diagonal(fill(σ_η^2, dimension_of_state))
    t_TOD = @elapsed L_TOD, K_TOD = Todorov_optimization_multidimensional_case(
        A, B, H, C, D, T, dimension_of_state, dimension_of_control,
        dimension_of_observation, x_state_mean_0, Σ_1_x, Σ_1_z,
        Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η, N_iter_Todorov)
    K_initial, L_initial = copy(K_TOD), copy(L_TOD)
    t_LCO = @elapsed res = Optimal_Control_Estimation_FullProcedure(
        K_initial, L_initial, A, B, H, C, D, T, Ω_ξ, Ω_ω, Ω_η,
        Q_matrix, R_matrix, Σ_1_x, Σ_1_z, x_1_mean, z_1_mean;
        ρ_list=[0.0], ρ_list_large=[1e-10,1e-8,1e-6],
        N_iter_tot=N_iter_tot, iters_K=1, iters_L=1,
        N_iter_GD_refine=0, N_iter_GD_from_IC=0,
        gd_method=:adam, η_GD=1e-3, maxiter_GD=50, tol_GD=1e-6,
        β1_Adam=0.9, β2_Adam=0.999, ϵ_Adam=1e-8,
        use_GD_refinement=false, use_GD_from_IC=false, verbose=verbose)
    return (; L_TOD, K_TOD, L_Full_Solver=res.L_best, K_Full_Solver=res.K_best,
        cost_Full_Solver=res.best_cost, best_cost_list=collect(res.best_cost_list),
        best_method=res.best_method, best_ρ=res.best_ρ, cost_init=res.cost_init,
        t_TOD, t_LCO)
end

function arm_analytical_cost(p, L, K, σ_η)
    return expected_cost_raw_moments_propagation(p.T, p.dimension_of_state,
        p.dimension_of_control, p.dimension_of_observation, K, L,
        p.A, p.B, p.H, p.C, p.D, p.x_1_mean, p.z_1_mean,
        p.Σ_1_x, p.Σ_1_z, p.Q_matrix, p.R_matrix, p.Ω_ξ, p.Ω_ω,
        Diagonal(fill(σ_η^2,p.dimension_of_state)))
end

function simulate_arm(p, L, K, σ_η, N_trials_MC; rng=MersenneTwister(1234))
    return monte_carlo_simulation(L, K, p.A, p.B, p.H, p.C, p.D, p.T,
        p.dimension_of_state, p.dimension_of_control, p.dimension_of_observation,
        p.x_1_mean, p.z_1_mean, p.Σ_1_x, p.Σ_1_z, p.Q_matrix, p.R_matrix,
        Matrix(p.Ω_ξ), Matrix(p.Ω_ω), Matrix(Diagonal(fill(σ_η^2,p.dimension_of_state))),
        N_trials_MC; rng=rng)
end

# Distribution-specific draws pass through the definitive core simulator.
# All five independent noise channels have mean zero and variance one.
# Initial x/z covariances are exactly zero in this experiment.
struct ArmNoiseRNG{R<:AbstractRNG} <: AbstractRNG
    rng::R
    noise_type::Symbol
    function ArmNoiseRNG(rng::R, noise_type::Symbol) where R<:AbstractRNG
        noise_type in (:gaussian, :studentt, :beta) || error("Unknown noise type: $noise_type")
        new{R}(rng,noise_type)
    end
end
function Random.randn(rng::ArmNoiseRNG)
    rng.noise_type === :gaussian && return randn(rng.rng)
    rng.noise_type === :studentt && return sqrt(3/5)*rand(rng.rng,TDist(5))
    α, β = 2, 5
    μ = α/(α+β)
    v = α*β/((α+β)^2*(α+β+1))
    return (rand(rng.rng,Beta(α,β))-μ)/sqrt(v)
end
Random.randn(rng::ArmNoiseRNG, n::Integer) = [randn(rng) for _ in 1:n]
