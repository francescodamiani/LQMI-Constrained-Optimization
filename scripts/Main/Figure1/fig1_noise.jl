# Figure 1e–g: adapted from RandomTask_sweep_internal_noise.jl.
using Optim, LinearAlgebra, Random, Statistics, Printf, DataFrames, JLD2

function simulate_internal_noise(data_dir; smoke=false)
random_seed_master = 1234
Random.seed!(random_seed_master)

# Choose ONE dimension and ONE seed (fixed problem instance)
n_fixed    = smoke ? 2 : 5
seed_fixed = 1

# Horizon
T = smoke ? 4 : 10
ρ_max = 0.9

# Noise (keep all fixed except σ_η)
σ_ξ, σ_ω, σ_ρ, σ_ϵ = 0.2, 0.2, 0.2, 0.2

# Sweep over σ_η
# The saved seven-level archive is the publication run; the supplied source
# had subsequently been edited to a single-noise exploratory run.
sigma_eta_list = smoke ? [0.0, 0.2, 0.7] : [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7]

# Monte Carlo settings
N_trials_MC = smoke ? 100 : 1000
#N_trials_MC = 10000

# TOD / LCO settings
N_iter_Todorov = smoke ? 3 : 500

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

use_GD_refinement = false
use_GD_from_IC    = false
N_iter_GD_refine  = 0
N_iter_GD_from_IC = 0

# GD baseline
algorithm_GD  = Optim.GradientDescent()
# 10,000 iterations identifies the archived publication run (filename and Table 1).
# The supplied exploratory source left this assignment unfinished.
iterations_GD = smoke ? 3 : 10000

options_GD = Optim.Options(
    iterations   = iterations_GD,
    store_trace  = false,
    show_trace   = false
)

# BUILD ONE FIXED PROBLEM INSTANCE -------------------------------------------------------------------------------------
rng_problem = MersenneTwister(random_seed_master + 10_000*n_fixed + seed_fixed)

(dimension_of_state, dimension_of_control, dimension_of_observation, Tloc,
 A, B, H, C, D, Ω_ξ, Ω_ω, Ω_η_dummy,
 x0, z0, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix) =
    build_random_square_task(n_fixed; T=T, ρ_max=ρ_max,
                             σ_ξ=σ_ξ, σ_ω=σ_ω, σ_ρ=σ_ρ, σ_ϵ=σ_ϵ, σ_η=0.0,
                             rng=rng_problem)

Tm1 = Tloc - 1
Nσ  = length(sigma_eta_list)

# DATAFRAMES / CONTAINERS ----------------------------------------------------------------------------------------------
# Per-method summary table
df_sweep = DataFrame(
    n = Int[],
    seed = Int[],
    sigma_eta = Float64[],
    method = String[],
    solve_time_sec = Float64[],
    mc_cost_mean = Float64[],
    mc_cost_sem  = Float64[],
    raw_objective = Float64[],
    N_trials_MC = Int[]
)

# Save all solutions for every σ_η
# Shapes:
#   L_* : (dim_u, dim_x, T-1, Nσ)
#   K_* : (dim_x, dim_y, T-1, Nσ)
#   M_* : (dim_x, dim_x, T-1, Nσ)
solutions_all = Dict{String, Any}()

solutions_all["sigma_eta_list"] = sigma_eta_list

solutions_all["L_TOD"] = zeros(dimension_of_control, dimension_of_state, Tm1, Nσ)
solutions_all["K_TOD"] = zeros(dimension_of_state, dimension_of_observation, Tm1, Nσ)
solutions_all["M_TOD"] = zeros(dimension_of_state, dimension_of_state, Tm1, Nσ)

solutions_all["L_LCO"] = zeros(dimension_of_control, dimension_of_state, Tm1, Nσ)
solutions_all["K_LCO"] = zeros(dimension_of_state, dimension_of_observation, Tm1, Nσ)
solutions_all["M_LCO"] = zeros(dimension_of_state, dimension_of_state, Tm1, Nσ)

solutions_all["L_GD"]  = zeros(dimension_of_control, dimension_of_state, Tm1, Nσ)
solutions_all["K_GD"]  = zeros(dimension_of_state, dimension_of_observation, Tm1, Nσ)
solutions_all["M_GD"]  = zeros(dimension_of_state, dimension_of_state, Tm1, Nσ)

# Optional packed summaries useful for debugging later
solutions_all["raw_objective_TOD"] = zeros(Nσ)
solutions_all["raw_objective_LCO"] = zeros(Nσ)
solutions_all["raw_objective_GD"]  = zeros(Nσ)

solutions_all["mc_cost_mean_TOD"] = zeros(Nσ)
solutions_all["mc_cost_mean_LCO"] = zeros(Nσ)
solutions_all["mc_cost_mean_GD"]  = zeros(Nσ)

solutions_all["mc_cost_sem_TOD"] = zeros(Nσ)
solutions_all["mc_cost_sem_LCO"] = zeros(Nσ)
solutions_all["mc_cost_sem_GD"]  = zeros(Nσ)

solutions_all["solve_time_sec_TOD"] = zeros(Nσ)
solutions_all["solve_time_sec_LCO"] = zeros(Nσ)
solutions_all["solve_time_sec_GD"]  = zeros(Nσ)

# Save the fixed problem instance too, so you can recompute any diagnostic later
problem_data = Dict{String, Any}()
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
problem_data["ρ_max"] = ρ_max
problem_data["σ_ξ"] = σ_ξ
problem_data["σ_ω"] = σ_ω
problem_data["σ_ρ"] = σ_ρ
problem_data["σ_ϵ"] = σ_ϵ
problem_data["random_seed_master"] = random_seed_master

# Full new-run diagnostics (historical files remain untouched).
cost_per_trial = Dict(method => zeros(N_trials_MC, Nσ) for method in ["TOD", "LCO", "GD"])
lco_diagnostics = Dict{String,Any}()
gd_diagnostics = Dict{String,Any}()
mc_seeds = [random_seed_master + 400_000 + round(Int, 1e6*σ) for σ in sigma_eta_list]

# RUN SWEEP ------------------------------------------------------------------------------------------------------------
println("\n====================== INTERNAL NOISE SWEEP (fixed problem) ======================\n")
println("Fixed problem: n=$n_fixed, seed=$seed_fixed")
println("σ_η list = ", sigma_eta_list)
println("N_trials_MC = $N_trials_MC\n")

for (iσ, σ_η) in enumerate(sigma_eta_list)
    println("\n------------------- σ_η = $σ_η -------------------\n")

    # Update ONLY internal noise covariance
    Ω_η = Diagonal(fill(σ_η^2, dimension_of_state))

    # Deterministic seed for MC evaluation
    random_seed_MC = random_seed_master + 400_000 + round(Int, 1e6*σ_η)

    # Raw analytical objective for this σ_η
    cost_fun(xvec) = expected_cost_using_mom_prop_for_GD_WHOLE_optimization_raw_moments(
        xvec, A, B, H, C, D, Tloc,
        dimension_of_state, dimension_of_control, dimension_of_observation,
        x_1_mean, z_1_mean, Σ_1_x, Σ_1_z,
        Q_matrix, R_matrix,
        Ω_ξ, Ω_ω, Ω_η
    )

    L_TOD_tmp = K_TOD_tmp = L_LCO_tmp = K_LCO_tmp = x_opt = nothing
    # ---- (1) TOD ----
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

    M_TOD = build_M(A, B, H, L_TOD_tmp, K_TOD_tmp, Tloc)

    raw_obj_TOD = cost_fun(pack_LK(L_TOD_tmp, K_TOD_tmp))

    # Dense covariances preserve values and support the singular σ_η=0 case on Julia 1.8.5.
    mc = monte_carlo_simulation(
        L_TOD_tmp, K_TOD_tmp, A, B, H, C, D, Tloc,
        dimension_of_state, dimension_of_control, dimension_of_observation,
        x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix,
        Matrix(Ω_ξ), Matrix(Ω_ω), Matrix(Ω_η), N_trials_MC; rng=MersenneTwister(random_seed_MC))
    mean_total_cost_TOD = mc.empirical_cost
    sem_total_cost_TOD = mc.empirical_cost_SE
    cost_per_trial["TOD"][:, iσ] = mc.cost_per_trial

    push!(df_sweep, (
        n_fixed, seed_fixed, σ_η, "TOD",
        t_TOD, mean_total_cost_TOD, sem_total_cost_TOD, raw_obj_TOD, N_trials_MC
    ))

    # Save TOD full solution
    solutions_all["L_TOD"][:, :, :, iσ] .= L_TOD_tmp
    solutions_all["K_TOD"][:, :, :, iσ] .= K_TOD_tmp
    solutions_all["M_TOD"][:, :, :, iσ] .= M_TOD
    solutions_all["raw_objective_TOD"][iσ] = raw_obj_TOD
    solutions_all["mc_cost_mean_TOD"][iσ] = mean_total_cost_TOD
    solutions_all["mc_cost_sem_TOD"][iσ] = sem_total_cost_TOD
    solutions_all["solve_time_sec_TOD"][iσ] = t_TOD

    # ---- (2) LCO (initialized from TOD) ----
    K_initial = copy(K_TOD_tmp)
    L_initial = copy(L_TOD_tmp)

    t_LCO = @elapsed begin
        K_LCO, L_LCO, cost_LCO_internal, method_LCO,
        ρ_LCO, cost_init_LCO, best_cost_list_LCO =
            Optimal_Control_Estimation_FullProcedure(
                K_initial, L_initial,
                A, B, H, C, D, Tloc,
                Ω_ξ, Ω_ω, Ω_η,
                Q_matrix, R_matrix,
                Σ_1_x, Σ_1_z, x_1_mean, z_1_mean;
                ρ_list = ρ_list, ρ_list_large = ρ_list_large,
                N_iter_tot = N_iter_tot, iters_K = iters_K, iters_L = iters_L,
                N_iter_GD_refine = N_iter_GD_refine, N_iter_GD_from_IC = N_iter_GD_from_IC,
                gd_method = gd_method, η_GD = η_GD,
                maxiter_GD = maxiter_GD, tol_GD = tol_GD,
                β1_Adam = β1_Adam, β2_Adam = β2_Adam, ϵ_Adam = ϵ_Adam,
                use_GD_refinement = use_GD_refinement, use_GD_from_IC = use_GD_from_IC,
                verbose = verbose
            )
        K_LCO_tmp = K_LCO
        L_LCO_tmp = L_LCO
    end

    lco_diagnostics[string(σ_η)] = Dict("best_cost_list"=>best_cost_list_LCO,
        "best_cost"=>cost_LCO_internal, "cost_init"=>cost_init_LCO,
        "best_method"=>method_LCO, "best_rho"=>ρ_LCO)
    M_LCO = build_M(A, B, H, L_LCO_tmp, K_LCO_tmp, Tloc)

    raw_obj_LCO = cost_fun(pack_LK(L_LCO_tmp, K_LCO_tmp))

    mc = monte_carlo_simulation(
        L_LCO_tmp, K_LCO_tmp, A, B, H, C, D, Tloc,
        dimension_of_state, dimension_of_control, dimension_of_observation,
        x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix,
        Matrix(Ω_ξ), Matrix(Ω_ω), Matrix(Ω_η), N_trials_MC; rng=MersenneTwister(random_seed_MC))
    mean_total_cost_LCO = mc.empirical_cost
    sem_total_cost_LCO = mc.empirical_cost_SE
    cost_per_trial["LCO"][:, iσ] = mc.cost_per_trial

    push!(df_sweep, (
        n_fixed, seed_fixed, σ_η, "LCO",
        t_LCO, mean_total_cost_LCO, sem_total_cost_LCO, raw_obj_LCO, N_trials_MC
    ))

    # Save LCO full solution
    solutions_all["L_LCO"][:, :, :, iσ] .= L_LCO_tmp
    solutions_all["K_LCO"][:, :, :, iσ] .= K_LCO_tmp
    solutions_all["M_LCO"][:, :, :, iσ] .= M_LCO
    solutions_all["raw_objective_LCO"][iσ] = raw_obj_LCO
    solutions_all["mc_cost_mean_LCO"][iσ] = mean_total_cost_LCO
    solutions_all["mc_cost_sem_LCO"][iσ] = sem_total_cost_LCO
    solutions_all["solve_time_sec_LCO"][iσ] = t_LCO

    # ---- (3) GD (initialized from LCO) ----
    x0_gd = pack_LK(L_LCO_tmp, K_LCO_tmp) #INITIAL CONDITONS!!! 
    #x0_gd = pack_LK(L_TOD_tmp, K_TOD_tmp)

    t_GD = @elapsed begin
        res = Optim.optimize(cost_fun, x0_gd, algorithm_GD, options_GD)
        x_opt = res.minimizer
    end

    gd_diagnostics[string(σ_η)] = Dict("iterations"=>Optim.iterations(res),
        "converged"=>Optim.converged(res), "minimum"=>Optim.minimum(res))
    L_GD, K_GD = unpack_LK(x_opt, dimension_of_control, dimension_of_state, dimension_of_observation, Tloc)
    M_GD = build_M(A, B, H, L_GD, K_GD, Tloc)

    raw_obj_GD = cost_fun(pack_LK(L_GD, K_GD))

    mc = monte_carlo_simulation(
        L_GD, K_GD, A, B, H, C, D, Tloc,
        dimension_of_state, dimension_of_control, dimension_of_observation,
        x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix,
        Matrix(Ω_ξ), Matrix(Ω_ω), Matrix(Ω_η), N_trials_MC; rng=MersenneTwister(random_seed_MC))
    mean_total_cost_GD = mc.empirical_cost
    sem_total_cost_GD = mc.empirical_cost_SE
    cost_per_trial["GD"][:, iσ] = mc.cost_per_trial

    push!(df_sweep, (
        n_fixed, seed_fixed, σ_η, "GD",
        t_GD, mean_total_cost_GD, sem_total_cost_GD, raw_obj_GD, N_trials_MC
    ))

    # Save GD full solution
    solutions_all["L_GD"][:, :, :, iσ] .= L_GD
    solutions_all["K_GD"][:, :, :, iσ] .= K_GD
    solutions_all["M_GD"][:, :, :, iσ] .= M_GD
    solutions_all["raw_objective_GD"][iσ] = raw_obj_GD
    solutions_all["mc_cost_mean_GD"][iσ] = mean_total_cost_GD
    solutions_all["mc_cost_sem_GD"][iσ] = sem_total_cost_GD
    solutions_all["solve_time_sec_GD"][iσ] = t_GD

    println(@sprintf(
        "σ_η=%.3g | solve_time(s): TOD=%.3g, LCO=%.3g, GD=%.3g | raw obj: TOD=%.8g, LCO=%.8g, GD=%.8g | cost mean±SEM: TOD=%.6g±%.3g, LCO=%.6g±%.3g, GD=%.6g±%.3g",
        σ_η, t_TOD, t_LCO, t_GD,
        raw_obj_TOD, raw_obj_LCO, raw_obj_GD,
        mean_total_cost_TOD, sem_total_cost_TOD,
        mean_total_cost_LCO, sem_total_cost_LCO,
        mean_total_cost_GD,  sem_total_cost_GD
    ))
end


    settings = Dict("smoke"=>smoke, "N_iter_Todorov"=>N_iter_Todorov,
        "N_iter_tot"=>N_iter_tot, "iterations_GD"=>iterations_GD,
        "MC_indexing"=>"dimension,time,trial", "GD_initialization"=>"LCO")
    JLD2.jldsave(joinpath(data_dir, "Fig_internal_noise_sweep_full_solutions_IC_at_LCO_solutions_10000_GD.jld2");
        df_sweep, solutions_all, problem_data, sigma_eta_list, n_fixed, seed_fixed, N_trials_MC,
        cost_per_trial, lco_diagnostics, gd_diagnostics, mc_seeds, settings)
    return nothing
end
