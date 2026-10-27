# Figure 1b,c: dimension benchmark with optional reduced validation settings.
using Optim, DataFrames, JLD2, Printf

function simulate_dimensions(data_dir; smoke=false)
random_seed_master = 1234
Random.seed!(random_seed_master)

# Horizon
T = smoke ? 4 : 10

# Noise parameters
ρ_max = 0.9
σ_ξ, σ_ω, σ_ρ, σ_ϵ, σ_η = 0.2, 0.2, 0.2, 0.2, 0.5

# Two n-grids
n_list_fast = smoke ? [1, 2] : [1, 2, 5, 10, 15, 20, 40, 60, 80, 100]  # TOD + LCO
n_list_gd   = smoke ? [1, 2] : [1, 2, 5, 10]                            # GD only

# Seeds / realizations
N_seeds_fast = smoke ? 2 : 10
N_seeds_gd   = smoke ? 2 : 10

# Todorov iterations (initializer)
N_iter_Todorov = smoke ? 3 : 500

# LCO FullProcedure parameters (match arm-task script)
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
iterations_GD = smoke ? 3 : 10000

options_GD_default = Optim.Options(
    iterations = iterations_GD,
    store_trace = false,
    show_trace = false
)

# RESULTS ------------------------------------------------------------------------------------------------------------

df = DataFrame(
    n = Int[],
    seed = Int[],
    method = String[],
    time_sec = Float64[],
    cost = Float64[]
)

solutions = Dict{String,Any}()

# RUN ---------------------------------------------------------------------------------------------------------------

println("\n====================== FIGURE 1 BENCHMARK (b,c) ======================\n")
println("n_list_fast = ", n_list_fast)
println("n_list_gd   = ", n_list_gd)
println("N_seeds_fast = ", N_seeds_fast, " | N_seeds_gd = ", N_seeds_gd, "\n")

for n in n_list_fast
    println("\n------------------- dimension n = $n -------------------\n")

    for s in 1:N_seeds_fast

        rng = MersenneTwister(random_seed_master + 10_000*n + s)

        (dim_x, dim_u, dim_y, Tloc,
         A, B, H, C, D, Ω_ξ, Ω_ω, Ω_η,
         x0, z0, x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Q, R) =
            build_random_square_task(n; T=T, ρ_max=ρ_max,
                                     σ_ξ=σ_ξ, σ_ω=σ_ω, σ_ρ=σ_ρ, σ_ϵ=σ_ϵ, σ_η=σ_η,
                                     rng=rng)

        L_TOD_tmp = K_TOD_tmp = nothing
        L_LCO_tmp = K_LCO_tmp = nothing
        x_opt = nothing
        # ---- (1) TOD ----
        t_TOD = @elapsed begin
            L_TOD, K_TOD = Todorov_optimization_multidimensional_case(
                A, B, H, C, D, Tloc,
                dim_x, dim_u, dim_y,
                x0, Σ_1_x, Σ_1_z,
                Q, R,
                Ω_ξ, Ω_ω, Ω_η,
                N_iter_Todorov
            )
            L_TOD_tmp = L_TOD
            K_TOD_tmp = K_TOD
        end

        cost_TOD = expected_cost_using_trace_formula_induction(
            L_TOD_tmp, K_TOD_tmp, A, B, H, C, D, Tloc, dim_x, x0, x0,
            Σ_1_x, Σ_1_z, Q, R, Ω_ξ, Ω_ω, Ω_η
        )

        push!(df, (n=n, seed=s, method="TOD", time_sec=t_TOD, cost=cost_TOD))

        # ---- (2) LCO ----
        K_initial = copy(K_TOD_tmp)
        L_initial = copy(L_TOD_tmp)

        t_LCO = @elapsed begin
            K_LCO, L_LCO, cost_LCO_internal, method_LCO,
            ρ_LCO, cost_init_LCO, best_cost_list_LCO =
                Optimal_Control_Estimation_FullProcedure(
                    K_initial, L_initial,
                    A, B, H, C, D, Tloc,
                    Ω_ξ, Ω_ω, Ω_η,
                    Q, R,
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

        cost_LCO = expected_cost_using_trace_formula_induction(
            L_LCO_tmp, K_LCO_tmp, A, B, H, C, D, Tloc, dim_x, x0, x0,
            Σ_1_x, Σ_1_z, Q, R, Ω_ξ, Ω_ω, Ω_η
        )

        push!(df, (n=n, seed=s, method="LCO", time_sec=t_LCO, cost=cost_LCO))

        # ---- (3) GD (subset of n) ----
        did_run_gd = (n in n_list_gd) && (s <= N_seeds_gd)
        if did_run_gd

            x0_gd = pack_LK(L_TOD_tmp, K_TOD_tmp) #Initial conditions of GD!

            cost_fun(xvec) = expected_cost_using_mom_prop_for_GD_WHOLE_optimization_raw_moments(
                xvec, A, B, H, C, D, Tloc,
                dim_x, dim_u, dim_y,
                x_1_mean, z_1_mean, Σ_1_x, Σ_1_z,
                Q, R,
                Ω_ξ, Ω_ω, Ω_η
            )

            t_GD = @elapsed begin
                res = Optim.optimize(cost_fun, x0_gd, algorithm_GD, options_GD_default)
                x_opt = res.minimizer
            end

            L_GD, K_GD = unpack_LK(x_opt, dim_u, dim_x, dim_y, Tloc)

            cost_GD = expected_cost_using_trace_formula_induction(
                L_GD, K_GD, A, B, H, C, D, Tloc, dim_x, x0, x0,
                Σ_1_x, Σ_1_z, Q, R, Ω_ξ, Ω_ω, Ω_η
            )

            push!(df, (n=n, seed=s, method="GD", time_sec=t_GD, cost=cost_GD))
        end

        record = Dict{String,Any}("L_TOD"=>L_TOD_tmp, "K_TOD"=>K_TOD_tmp,
            "L_LCO"=>L_LCO_tmp, "K_LCO"=>K_LCO_tmp,
            "best_cost_list"=>best_cost_list_LCO, "best_cost"=>cost_LCO_internal,
            "best_method"=>method_LCO, "best_rho"=>ρ_LCO)
        if did_run_gd
            record["L_GD"] = L_GD
            record["K_GD"] = K_GD
            record["GD_iterations"] = Optim.iterations(res)
            record["GD_converged"] = Optim.converged(res)
        end
        solutions["n$(n)_seed$(s)"] = record
        println(@sprintf("n=%d seed=%d | time(s): TOD=%.3g, LCO=%.3g%s | cost: TOD=%.6g, LCO=%.6g%s",
            n, s, t_TOD, t_LCO,
            did_run_gd ? @sprintf(", GD=%.3g", t_GD) : " (GD skipped)",
            cost_TOD, cost_LCO,
            did_run_gd ? @sprintf(", GD=%.6g", cost_GD) : ""
        ))
    end
end


    df_min = select(df, [:n, :seed, :method, :time_sec, :cost])
    settings = Dict("smoke"=>smoke, "T"=>T, "n_list_fast"=>n_list_fast,
        "n_list_gd"=>n_list_gd, "N_seeds_fast"=>N_seeds_fast, "N_seeds_gd"=>N_seeds_gd,
        "random_seed_master"=>random_seed_master, "N_iter_Todorov"=>N_iter_Todorov,
        "N_iter_tot"=>N_iter_tot, "iterations_GD"=>iterations_GD)
    JLD2.jldsave(joinpath(data_dir, "Fig1_ab_data_10000_GD.jld2"); df_min, solutions, settings)
    return nothing
end
