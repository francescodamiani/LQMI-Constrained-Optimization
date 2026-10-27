# Figure 1d: adapted from Monotonic_Convergence_LCO.jl.
using Optim, LinearAlgebra, Random, Statistics, Printf, DataFrames, JLD2

function simulate_convergence(data_dir; smoke=false)
random_seed_master = 1234
Random.seed!(random_seed_master)

# fixed problem instance
n_fixed    = smoke ? 2 : 5
seed_fixed = 1

# horizon
T = smoke ? 4 : 10
ρ_max = 0.9

# fixed noise parameters except σ_η
σ_ξ, σ_ω, σ_ρ, σ_ϵ = 0.2, 0.2, 0.2, 0.2

# internal-noise sweep
sigma_eta_list = smoke ? [0.0, 0.2, 0.7] : [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7]

# TOD settings
N_iter_Todorov = smoke ? 3 : 500

# LCO settings
ρ_list       = [0.0]
ρ_list_large = [1e-10, 1e-8, 1e-6]
N_iter_tot   = smoke ? 3 : 500
iters_K      = 1
iters_L      = 1

gd_method  = :adam
η_GD       = 1e-3
maxiter_GD = 50
tol_GD     = 1e-6
β1_Adam    = 0.9
β2_Adam    = 0.999
ϵ_Adam     = 1e-8
verbose = false

# IMPORTANT for this figure:
use_GD_refinement = false
use_GD_from_IC    = false
N_iter_GD_refine  = 0
N_iter_GD_from_IC = 0

# BUILD ONE FIXED PROBLEM INSTANCE -------------------------------------------------------------------------------------
rng_problem = MersenneTwister(random_seed_master + 10_000*n_fixed + seed_fixed)

(dimension_of_state, dimension_of_control, dimension_of_observation, Tloc,
 A, B, H, C, D, Ω_ξ, Ω_ω, Ω_η_dummy,
 x0, z0, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix) =
    build_random_square_task(n_fixed; T=T, ρ_max=ρ_max,
                             σ_ξ=σ_ξ, σ_ω=σ_ω, σ_ρ=σ_ρ, σ_ϵ=σ_ϵ, σ_η=0.0,
                             rng=rng_problem)

Nσ = length(sigma_eta_list)

# CONTAINERS -----------------------------------------------------------------------------------------------------------
df_summary = DataFrame(
    sigma_eta = Float64[],
    n_iter_recorded = Int[],
    cost_init = Float64[],
    cost_final = Float64[],
    solve_time_sec = Float64[],
    best_method = String[],
    best_rho = Float64[]
)

cost_trajectories = Dict{String,Any}()
cost_trajectories["sigma_eta_list"] = sigma_eta_list
cost_trajectories["traj"] = Vector{Vector{Float64}}(undef, Nσ)

problem_data = Dict{String,Any}()
problem_data["A"] = A
problem_data["B"] = B
problem_data["H"] = H
problem_data["C"] = C
problem_data["D"] = D
problem_data["Q_matrix"] = Q_matrix
problem_data["R_matrix"] = R_matrix
problem_data["Ω_ξ"] = Ω_ξ
problem_data["Ω_ω"] = Ω_ω
problem_data["x0"] = x0
problem_data["z0"] = z0
problem_data["x_1_mean"] = x_1_mean
problem_data["z_1_mean"] = z_1_mean
problem_data["Σ_1_x"] = Σ_1_x
problem_data["Σ_1_z"] = Σ_1_z
problem_data["n_fixed"] = n_fixed
problem_data["seed_fixed"] = seed_fixed
problem_data["Tloc"] = Tloc

# RUN ---------------------------------------------------------------------------------------------------------------
println("\n====================== LCO COST ACROSS ITERATIONS ======================\n")
println("Fixed problem: n=$n_fixed, seed=$seed_fixed")
println("σ_η list = ", sigma_eta_list, "\n")

for (iσ, σ_η) in enumerate(sigma_eta_list)

    println("\n------------------- σ_η = $σ_η -------------------\n")

    Ω_η = Diagonal(fill(σ_η^2, dimension_of_state))

    L_TOD_tmp = K_TOD_tmp = res_LCO = nothing
    # ---- TOD initializer ----
    t_TOD = @elapsed begin
        L_TOD, K_TOD = Todorov_optimization_multidimensional_case(
            A, B, H, C, D, Tloc,
            dimension_of_state, dimension_of_control, dimension_of_observation,
            x0, Σ_1_x, Σ_1_z,
            Q_matrix, R_matrix,
            Ω_ξ, Ω_ω, Ω_η,
            N_iter_Todorov
        )
        L_TOD_tmp = L_TOD
        K_TOD_tmp = K_TOD
    end

    # ---- LCO from TOD ----
    K_initial = copy(K_TOD_tmp)
    L_initial = copy(L_TOD_tmp)

    t_LCO = @elapsed begin
        res = Optimal_Control_Estimation_FullProcedure(
            K_initial, L_initial,
            A, B, H, C, D, Tloc,
            Ω_ξ, Ω_ω, Ω_η,
            Q_matrix, R_matrix,
            Σ_1_x, Σ_1_z, x_1_mean, z_1_mean;
            ρ_list = ρ_list,
            ρ_list_large = ρ_list_large,
            N_iter_tot = N_iter_tot,
            iters_K = iters_K,
            iters_L = iters_L,
            N_iter_GD_refine = N_iter_GD_refine,
            N_iter_GD_from_IC = N_iter_GD_from_IC,
            gd_method = gd_method,
            η_GD = η_GD,
            maxiter_GD = maxiter_GD,
            tol_GD = tol_GD,
            β1_Adam = β1_Adam,
            β2_Adam = β2_Adam,
            ϵ_Adam = ϵ_Adam,
            use_GD_refinement = use_GD_refinement,
            use_GD_from_IC = use_GD_from_IC,
            verbose = verbose
        )
        res_LCO = res
    end

    cost_traj = collect(res_LCO.best_cost_list)

    cost_trajectories["traj"][iσ] = cost_traj

    push!(df_summary, (
        σ_η,
        length(cost_traj),
        cost_traj[1],
        cost_traj[end],
        t_LCO,
        res_LCO.best_method,
        res_LCO.best_ρ
    ))

    println(@sprintf(
        "σ_η=%.3f | TOD time=%.3g s | LCO time=%.3g s | n_iter=%d | cost_init=%.10g | cost_final=%.10g | method=%s | ρ=%.3g",
        σ_η, t_TOD, t_LCO, length(cost_traj), cost_traj[1], cost_traj[end],
        res_LCO.best_method, res_LCO.best_ρ
    ))
end


    settings = Dict("smoke"=>smoke, "random_seed_master"=>random_seed_master,
        "N_iter_Todorov"=>N_iter_Todorov, "N_iter_tot"=>N_iter_tot)
    JLD2.jldsave(joinpath(data_dir, "Fig1d_LCO_cost_trajectories.jld2");
        cost_trajectories, df_summary, problem_data, settings)
    return nothing
end
