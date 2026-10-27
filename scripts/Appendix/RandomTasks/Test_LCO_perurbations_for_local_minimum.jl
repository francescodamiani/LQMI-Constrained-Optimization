# LCO local-optimality experiment: relative perturbations of control and estimation gains.
# The loop order, relative perturbations, RNG streams and solver settings are unchanged.
using LinearAlgebra, Random, Printf, Statistics, DataFrames, JLD2

function simulate_perturbations(data_dir; n_fixed=5, T=10,
    N_iter_Todorov=500, N_iter_tot=500,
    N_dir_L=100000, N_dir_K=100000, N_dir_joint=100000,
    epsilon_list=[1e-4, 1e-3, 1e-2, 1e-1])

# PARAMETERS ------------------------------------------------------------------------------------------------------------
random_seed_master = 1234
Random.seed!(random_seed_master)

seed_fixed = 1

ρ_max = 0.9

σ_ξ, σ_ω, σ_ρ, σ_ϵ = 0.2, 0.2, 0.2, 0.2
σ_η = 0.2


ρ_list       = [0.0]
ρ_list_large = [1e-10, 1e-8, 1e-6]
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

# perturbation scan settings
rng_scan = MersenneTwister(20260328)

# good first-pass values; if fast, raise to 500 or 1000

# relative perturbation sizes: from 0.01% to 10% of baseline norm

# BUILD ONE FIXED PROBLEM INSTANCE -------------------------------------------------------------------------------------
rng_problem = MersenneTwister(random_seed_master + 10_000*n_fixed + seed_fixed)

(dimension_of_state, dimension_of_control, dimension_of_observation, Tloc,
 A, B, H, C, D, Ω_ξ, Ω_ω, Ω_η_dummy,
 x0, z0, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q_matrix, R_matrix) =
    build_random_square_task(n_fixed; T=T, ρ_max=ρ_max,
                             σ_ξ=σ_ξ, σ_ω=σ_ω, σ_ρ=σ_ρ, σ_ϵ=σ_ϵ, σ_η=0.0,
                             rng=rng_problem)

Ω_η = Diagonal(fill(σ_η^2, dimension_of_state))

println("\n====================== LOCAL PERTURBATION SCAN AROUND LCO ======================\n")
println("Fixed problem: n=$n_fixed, seed=$seed_fixed")
println("σ_η = $σ_η")
println("N_dir_L = $N_dir_L, N_dir_K = $N_dir_K, N_dir_joint = $N_dir_joint")
println("epsilon_list = ", epsilon_list, "\n")

# COMPUTE TOD -----------------------------------------------------------------------------------------------------------
println("Computing TOD solution...")

t_TOD = @elapsed begin
    L_TOD_tmp, K_TOD_tmp = Todorov_optimization_multidimensional_case(
        A, B, H, C, D, Tloc,
        dimension_of_state, dimension_of_control, dimension_of_observation,
        x0, Σ_1_x, Σ_1_z,
        Q_matrix, R_matrix,
        Ω_ξ, Ω_ω, Ω_η,
        N_iter_Todorov
    )
end

println(@sprintf("TOD done in %.3g s", t_TOD))

# COMPUTE LCO -----------------------------------------------------------------------------------------------------------
println("Computing LCO solution...")

K_initial = copy(K_TOD_tmp)
L_initial = copy(L_TOD_tmp)

t_LCO = @elapsed begin
    K_LCO_tmp, L_LCO_tmp, cost_LCO_internal, method_LCO,
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
end

println(@sprintf("LCO done in %.3g s", t_LCO))

# BASELINE COST ---------------------------------------------------------------------------------------------------------
cost_LCO = safe_cost_eval(
    L_LCO_tmp, K_LCO_tmp,
    A, B, H, C, D, Tloc, dimension_of_state,
    x_1_mean, z_1_mean, Σ_1_x, Σ_1_z,
    Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η
)

println(@sprintf("Analytical cost at LCO = %.12g\n", cost_LCO))

# PREP ------------------------------------------------------------------------------------------------------------------
norm_L0 = max(norm(L_LCO_tmp), 1e-12)
norm_K0 = max(norm(K_LCO_tmp), 1e-12)
x_LCO = pack_LK(L_LCO_tmp, K_LCO_tmp)
norm_x0 = max(norm(x_LCO), 1e-12)

df_scan = DataFrame(
    perturbation_type = String[],
    direction_id = Int[],
    sign = Int[],
    epsilon = Float64[],
    epsilon_label = String[],
    cost_perturbed = Float64[],
    delta_cost = Float64[],
    relative_delta_cost = Float64[]
)

# RUN SCAN -------------------------------------------------------------------------------------------------------------
println("Running perturbation scan...\n")

# L only
for dir_id in 1:N_dir_L
    ΔL = randn(rng_scan, size(L_LCO_tmp)...)
    normalize_fro!(ΔL)

    for ε in epsilon_list
        ε_label = epsilon_to_label(ε)
        for sgn in (-1, 1)
            Lp = L_LCO_tmp .+ sgn * ε * norm_L0 * ΔL
            # Lp = L_LCO_tmp .+ sgn * ε * ΔL #absolute scan
            Kp = K_LCO_tmp

            cpert = safe_cost_eval(
                Lp, Kp,
                A, B, H, C, D, Tloc, dimension_of_state,
                x_1_mean, z_1_mean, Σ_1_x, Σ_1_z,
                Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η
            )

            Δc = cpert - cost_LCO
            relΔc = Δc / max(abs(cost_LCO), 1e-12)

            push!(df_scan, ("L_only", dir_id, sgn, ε, ε_label, cpert, Δc, relΔc))
        end
    end
end

# K only
for dir_id in 1:N_dir_K
    ΔK = randn(rng_scan, size(K_LCO_tmp)...)
    normalize_fro!(ΔK)

    for ε in epsilon_list
        ε_label = epsilon_to_label(ε)
        for sgn in (-1, 1)
            Lp = L_LCO_tmp
            Kp = K_LCO_tmp .+ sgn * ε * norm_K0 * ΔK 
            # Kp = K_LCO_tmp .+ sgn * ε * ΔK #absolute scan

            cpert = safe_cost_eval(
                Lp, Kp,
                A, B, H, C, D, Tloc, dimension_of_state,
                x_1_mean, z_1_mean, Σ_1_x, Σ_1_z,
                Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η
            )

            Δc = cpert - cost_LCO
            relΔc = Δc / max(abs(cost_LCO), 1e-12)

            push!(df_scan, ("K_only", dir_id, sgn, ε, ε_label, cpert, Δc, relΔc))
        end
    end
end

# joint
for dir_id in 1:N_dir_joint
    Δx = randn(rng_scan, length(x_LCO))
    normalize_vec!(Δx)

    for ε in epsilon_list
        ε_label = epsilon_to_label(ε)
        for sgn in (-1, 1)
            xp = x_LCO .+ sgn * ε * norm_x0 * Δx
            # xp = x_LCO .+ sgn * ε * Δx #absolute scan
            Lp, Kp = unpack_LK(xp, dimension_of_control, dimension_of_state, dimension_of_observation, Tloc)

            cpert = safe_cost_eval(
                Lp, Kp,
                A, B, H, C, D, Tloc, dimension_of_state,
                x_1_mean, z_1_mean, Σ_1_x, Σ_1_z,
                Q_matrix, R_matrix, Ω_ξ, Ω_ω, Ω_η
            )

            Δc = cpert - cost_LCO
            relΔc = Δc / max(abs(cost_LCO), 1e-12)

            push!(df_scan, ("joint", dir_id, sgn, ε, ε_label, cpert, Δc, relΔc))
        end
    end
end

println("Perturbation scan finished.")
println(@sprintf("Total evaluations: %d\n", nrow(df_scan)))

# SUMMARY --------------------------------------------------------------------------------------------------------------
df_summary = combine(
    groupby(df_scan, [:perturbation_type, :epsilon, :epsilon_label]),
    :delta_cost => safe_min => :min_delta_cost,
    :delta_cost => (x -> safe_quantile(x, 0.05)) => :p05_delta_cost,
    :delta_cost => safe_median => :median_delta_cost,
    :delta_cost => safe_mean => :mean_delta_cost,
    :delta_cost => safe_max => :max_delta_cost,
    :delta_cost => safe_frac_negative => :frac_negative
)

sort!(df_summary, [:perturbation_type, :epsilon])

println("Summary:\n")
show(df_summary, allrows=true, allcols=true)
println("\n")

# SAVE -----------------------------------------------------------------------------------------------------------------
# Simulation results include settings and the reference gains.
validation_only = !(n_fixed == 5 && T == 10 && N_iter_Todorov == 500 &&
    N_iter_tot == 500 && N_dir_L == 100000 && N_dir_K == 100000 &&
    N_dir_joint == 100000 && epsilon_list == [1e-4, 1e-3, 1e-2, 1e-1])
settings = (;random_seed_master, scan_seed=20260328, T, ρ_max,
    σ_ξ, σ_ω, σ_ρ, σ_ϵ, N_iter_Todorov, N_iter_tot, ρ_list, ρ_list_large,
    iters_K, iters_L, gd_method, η_GD, maxiter_GD, tol_GD, β1_Adam, β2_Adam,
    ϵ_Adam, use_GD_refinement, use_GD_from_IC, N_iter_GD_refine, N_iter_GD_from_IC,
    validation_only, perturbation_scaling="relative")
problem_data = (;A, B, H, C, D, Ω_ξ, Ω_ω, Ω_η, x_1_mean, z_1_mean,
    Σ_1_x, Σ_1_z, Q_matrix, R_matrix, L_LCO_tmp, K_LCO_tmp,
    cost_LCO_internal, method_LCO, ρ_LCO, cost_init_LCO, best_cost_list_LCO)
results = Dict{String,Any}("df_scan"=>df_scan, "df_summary"=>df_summary,
    "cost_LCO"=>cost_LCO, "epsilon_list"=>epsilon_list, "σ_η"=>σ_η,
    "n_fixed"=>n_fixed, "seed_fixed"=>seed_fixed, "N_dir_L"=>N_dir_L,
    "N_dir_K"=>N_dir_K, "N_dir_joint"=>N_dir_joint,
    "settings"=>settings, "problem_data"=>problem_data)
isdir(data_dir) || mkpath(data_dir)
tag_sigma = replace(@sprintf("%.3f", σ_η), "." => "p")
jld_out = joinpath(data_dir, "LCO_local_perturbation_scan_sigma_eta_" * tag_sigma * ".jld2")
JLD2.jldopen(jld_out, "w"; compress=true) do file
    for (key, value) in results
        file[key] = value
    end
end
println("Saved data to: ", jld_out)
return results
end
