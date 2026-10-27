# Under-observed random-task benchmark: TOD, LCO and direct GD.
using Optim, LinearAlgebra, Random, Printf, DataFrames, JLD2
include(joinpath(@__DIR__,"..","..","..","Main","Figure1","fig1_tasks.jl"))

function build_random_underobserved_task(
    m::Int;
    T::Int=10,
    ρ_max=0.9,
    σ_ξ=0.2,
    σ_ω=0.2,
    σ_ρ=0.2,
    σ_ϵ=0.2,
    σ_η=0.5,
    rng=Random.default_rng()
)

    dim_x = m
    dim_u = m
    dim_y = 1

    A = generate_stable_matrix(m; ρ_max, rng=rng)
    B = generate_stable_matrix(m; ρ_max, rng=rng)

    H = zeros(dim_y, dim_x)
    H[1,1] = 1.0

    D_tmp_full =
        abs.(generate_stable_matrix(
            m;
            ρ_max,
            rng=rng
        ))

    D = σ_ρ .* D_tmp_full[1:1,:]
    C = σ_ϵ .* B

    Ω_ω = Diagonal(fill(σ_ω^2, dim_y))
    Ω_ξ = Diagonal(fill(σ_ξ^2, dim_x))
    Ω_η = Diagonal(fill(σ_η^2, dim_x))

    x0 = ones(dim_x)
    z0 = copy(x0)
    x_1_mean = copy(x0)
    z_1_mean = copy(z0)

    Σ_1_x = zeros(dim_x, dim_x)
    Σ_1_z = zeros(dim_x, dim_x)

    Q = zeros(dim_x, dim_x, T)
    R = zeros(dim_u, dim_u, T)

    for t in 1:T
        Q[:,:,t] .= I(dim_x)
        R[:,:,t] .= I(dim_u)
    end

    R[:,:,T] .= zeros(dim_u, dim_u)
    Q[:,:,1] .= zeros(dim_x, dim_x)

    return (
        dim_x, dim_u, dim_y, T,
        A, B, H, C, D,
        Ω_ξ, Ω_ω, Ω_η,
        x0, z0,
        x_1_mean, z_1_mean,
        Σ_1_x, Σ_1_z,
        Q, R
    )

end

function simulate_underobserved_random(output_directory; validation_only=false)
random_seed_master = 1234
Random.seed!(random_seed_master)

dimension_list = validation_only ? [2] : [6, 8, 10]
T = validation_only ? 4 : 10

ρ_max = 0.9
σ_ξ = 0.2
σ_ω = 0.2
σ_ρ = 0.2
σ_ϵ = 0.2
σ_η = 0.5

N_iter_Todorov = validation_only ? 3 : 500

ρ_list = [0.0]
ρ_list_large = [1e-10, 1e-8, 1e-6]
N_iter_tot = validation_only ? 3 : 500
iters_K = 1
iters_L = 1

gd_method = :adam
η_GD = 1e-3
maxiter_GD = 50
tol_GD = 1e-6
β1_Adam = 0.9
β2_Adam = 0.999
ϵ_Adam = 1e-8

verbose = false

use_GD_refinement = false
use_GD_from_IC = false
N_iter_GD_refine = 0
N_iter_GD_from_IC = 0

algorithm_GD = GradientDescent()
iterations_GD = validation_only ? 2 : 10000

options_GD = Optim.Options(
    iterations=iterations_GD,
    store_trace=false,
    show_trace=false
)

#The LCO routine has a fixed pass count rather than an early-stopping test.
#This tolerance is used only to report a final relative convergence residual.
lco_residual_tolerance = 1e-6
lco_monotonicity_relative_tolerance = 1e-10

# OUTPUT AND RESUMABLE STORAGE ------------------------------------------------------------------------------------------------------------

full_folder_save_results =
    joinpath(output_directory, "")

mkpath(full_folder_save_results)

results_file = joinpath(output_directory, "random_task_severely_underobserved_results.jld2")

df_results = DataFrame(
    m = Int[],
    p = Int[],
    k = Int[],
    T = Int[],
    seed = Int[],
    completed = Bool[],
    failure_message = String[],
    observability_rank = Int[],
    observability_condition_number = Float64[],
    observability_minimum_singular_value = Float64[],
    observability_rank_tolerance = Float64[],
    cost_TOD = Float64[],
    cost_LCO = Float64[],
    cost_GD = Float64[],
    relative_LCO_GD_cost_gap_percent = Float64[],
    relative_LCO_improvement_over_TOD_percent = Float64[],
    runtime_TOD_s = Float64[],
    runtime_LCO_s = Float64[],
    runtime_GD_s = Float64[],
    LCO_updates = Int[],
    LCO_passes = Int[],
    LCO_convergence_residual = Float64[],
    LCO_converged = Bool[],
    LCO_maximum_positive_cost_increment = Float64[],
    LCO_monotonic = Bool[],
    LCO_selected_method = String[],
    LCO_selected_rho = Float64[],
    GD_iterations = Int[],
    GD_converged = Bool[],
    GD_hit_iteration_limit = Bool[],
    all_costs_and_runtimes_finite = Bool[]
)

lco_cost_histories =
    Dict{Int,Vector{Float64}}()

solutions = Dict{Int,Any}()

if isfile(results_file)
    saved_results =
        JLD2.load(results_file)
    get(saved_results,"validation_only",false) == validation_only || error("Cannot mix validation and publication runs")
    solutions = get(saved_results,"solutions",Dict{Int,Any}())

    df_results =
        saved_results["df_results"]

    lco_cost_histories =
        saved_results["lco_cost_histories"]

    println(
        "Loaded ",
        nrow(df_results),
        " saved dimension rows."
    )
end

function save_progress()

    JLD2.jldsave(
        results_file;
        generated_by_current_core=true,
        validation_only=validation_only,
        df_results=df_results,
        lco_cost_histories=lco_cost_histories,
        solutions=solutions,
        dimension_list=dimension_list,
        random_seed_master=random_seed_master,
        T=T,
        ρ_max=ρ_max,
        σ_ξ=σ_ξ,
        σ_ω=σ_ω,
        σ_ρ=σ_ρ,
        σ_ϵ=σ_ϵ,
        σ_η=σ_η,
        iterations_GD=iterations_GD,
        N_iter_Todorov=N_iter_Todorov,
        N_iter_tot=N_iter_tot,
        lco_residual_tolerance=
            lco_residual_tolerance,
        lco_monotonicity_relative_tolerance=
            lco_monotonicity_relative_tolerance
    )



end

# RUN ---------------------------------------------------------------------------------------------------------------

println(
    "\n================ SEVERELY UNDER-OBSERVED RANDOM BENCHMARK ================\n"
)

println("Dimensions: ", dimension_list)
println("p = m, k = 1, T = ", T)
println("H = [1, 0, ..., 0]")
println("GD iterations: ", iterations_GD)

for m_task in dimension_list

    completed_rows =
        findall(
            (df_results.m .== m_task) .&
            df_results.completed
        )

    if !isempty(completed_rows)
        println(
            "\nDimension m = ",
            m_task,
            " already completed."
        )
        continue
    end

    if nrow(df_results) > 0
        filter!(
            row -> row.m != m_task,
            df_results
        )
    end

    seed =
        random_seed_master +
        10_000*m_task +
        1

    println(
        "\n------------------- dimension m = ",
        m_task,
        " | seed = ",
        seed,
        " -------------------\n"
    )

    try

        rng = MersenneTwister(seed)

        (
            dim_x, dim_u, dim_y, Tloc,
            A, B, H, C, D,
            Ω_ξ, Ω_ω, Ω_η,
            x0, z0,
            x_1_mean, z_1_mean,
            Σ_1_x, Σ_1_z,
            Q, R
        ) =
            build_random_underobserved_task(
                m_task;
                T=T,
                ρ_max=ρ_max,
                σ_ξ=σ_ξ,
                σ_ω=σ_ω,
                σ_ρ=σ_ρ,
                σ_ϵ=σ_ϵ,
                σ_η=σ_η,
                rng=rng
            )

        observability =
            observability_diagnostics(
                A, H, Tloc
            )

        println(
            "Observability rank: ",
            observability.numerical_rank,
            "/",
            m_task,
            " | condition number: ",
            observability.condition_number
        )

        local L_TOD
        local K_TOD

        runtime_TOD =
            @elapsed begin
                L_TOD, K_TOD =
                    Todorov_optimization_multidimensional_case(
                        A, B, H, C, D, Tloc,
                        dim_x, dim_u, dim_y,
                        x0, Σ_1_x, Σ_1_z,
                        Q, R,
                        Ω_ξ, Ω_ω, Ω_η,
                        N_iter_Todorov
                    )
            end

        cost_TOD =
            expected_cost_using_trace_formula_induction(
                L_TOD, K_TOD,
                A, B, H, C, D,
                Tloc, dim_x,
                x0, x0,
                Σ_1_x, Σ_1_z,
                Q, R,
                Ω_ξ, Ω_ω, Ω_η
            )

        K_initial_LCO = copy(K_TOD)
        L_initial_LCO = copy(L_TOD)

        K_initial_GD = copy(K_TOD)
        L_initial_GD = copy(L_TOD)

        local result_LCO

        runtime_LCO =
            @elapsed begin
                result_LCO =
                    Optimal_Control_Estimation_FullProcedure(
                        K_initial_LCO,
                        L_initial_LCO,
                        A, B, H, C, D, Tloc,
                        Ω_ξ, Ω_ω, Ω_η,
                        Q, R,
                        Σ_1_x, Σ_1_z,
                        x_1_mean, z_1_mean;
                        ρ_list=ρ_list,
                        ρ_list_large=ρ_list_large,
                        N_iter_tot=N_iter_tot,
                        iters_K=iters_K,
                        iters_L=iters_L,
                        N_iter_GD_refine=
                            N_iter_GD_refine,
                        N_iter_GD_from_IC=
                            N_iter_GD_from_IC,
                        gd_method=gd_method,
                        η_GD=η_GD,
                        maxiter_GD=maxiter_GD,
                        tol_GD=tol_GD,
                        β1_Adam=β1_Adam,
                        β2_Adam=β2_Adam,
                        ϵ_Adam=ϵ_Adam,
                        use_GD_refinement=
                            use_GD_refinement,
                        use_GD_from_IC=
                            use_GD_from_IC,
                        verbose=verbose
                    )
            end

        K_LCO = result_LCO.K_best
        L_LCO = result_LCO.L_best

        cost_LCO =
            expected_cost_using_trace_formula_induction(
                L_LCO, K_LCO,
                A, B, H, C, D,
                Tloc, dim_x,
                x0, x0,
                Σ_1_x, Σ_1_z,
                Q, R,
                Ω_ξ, Ω_ω, Ω_η
            )

        lco_cost_history =
            collect(result_LCO.best_cost_list)

        lco_cost_histories[m_task] =
            lco_cost_history

        LCO_updates =
            length(lco_cost_history)-1

        LCO_passes =
            fld(LCO_updates, 2)

        LCO_convergence_residual =
            LCO_updates > 0 ?
            abs(
                lco_cost_history[end] -
                lco_cost_history[end-1]
            )/
            max(
                abs(lco_cost_history[end-1]),
                eps(Float64)
            ) :
            Inf

        cost_increments =
            diff(lco_cost_history)

        maximum_positive_increment =
            isempty(cost_increments) ?
            0.0 :
            max(0.0, maximum(cost_increments))

        monotonicity_tolerance =
            lco_monotonicity_relative_tolerance*
            max(
                abs(lco_cost_history[1]),
                1.0
            )

        LCO_monotonic =
            maximum_positive_increment <=
            monotonicity_tolerance

        LCO_converged =
            isfinite(LCO_convergence_residual) &&
            LCO_convergence_residual <=
            lco_residual_tolerance

        x0_GD =
            pack_LK(
                L_initial_GD,
                K_initial_GD
            )

        cost_function_GD(x_vector) =
            expected_cost_using_mom_prop_for_GD_WHOLE_optimization_raw_moments(
                x_vector,
                A, B, H, C, D, Tloc,
                dim_x, dim_u, dim_y,
                x_1_mean, z_1_mean,
                Σ_1_x, Σ_1_z,
                Q, R,
                Ω_ξ, Ω_ω, Ω_η
            )

        local result_GD

        runtime_GD =
            @elapsed begin
                result_GD =
                    optimize(
                        cost_function_GD,
                        x0_GD,
                        algorithm_GD,
                        options_GD
                    )
            end

        L_GD, K_GD =
            unpack_LK(
                result_GD.minimizer,
                dim_u,
                dim_x,
                dim_y,
                Tloc
            )

        cost_GD =
            expected_cost_using_trace_formula_induction(
                L_GD, K_GD,
                A, B, H, C, D,
                Tloc, dim_x,
                x0, x0,
                Σ_1_x, Σ_1_z,
                Q, R,
                Ω_ξ, Ω_ω, Ω_η
            )

        solutions[m_task] = (;L_TOD,K_TOD,L_LCO,K_LCO,L_GD,K_GD)

        GD_iterations =
            Optim.iterations(result_GD)

        GD_converged =
            Optim.converged(result_GD)

        GD_hit_iteration_limit =
            !GD_converged &&
            GD_iterations >= iterations_GD

        relative_LCO_GD_gap =
            100*(cost_LCO-cost_GD)/
            max(abs(cost_GD), eps(Float64))

        relative_LCO_improvement =
            100*(cost_TOD-cost_LCO)/
            max(abs(cost_TOD), eps(Float64))

        all_finite =
            all(isfinite, [
                cost_TOD,
                cost_LCO,
                cost_GD,
                runtime_TOD,
                runtime_LCO,
                runtime_GD,
                LCO_convergence_residual,
                maximum_positive_increment
            ])

        push!(
            df_results,
            (
                m_task,
                dim_u,
                dim_y,
                Tloc,
                seed,
                true,
                "",
                observability.numerical_rank,
                observability.condition_number,
                observability.minimum_singular_value,
                observability.rank_tolerance,
                cost_TOD,
                cost_LCO,
                cost_GD,
                relative_LCO_GD_gap,
                relative_LCO_improvement,
                runtime_TOD,
                runtime_LCO,
                runtime_GD,
                LCO_updates,
                LCO_passes,
                LCO_convergence_residual,
                LCO_converged,
                maximum_positive_increment,
                LCO_monotonic,
                result_LCO.best_method,
                result_LCO.best_ρ,
                GD_iterations,
                GD_converged,
                GD_hit_iteration_limit,
                all_finite
            )
        )

        sort!(df_results, :m)
        save_progress()

        flags = String[]

        observability.numerical_rank < m_task &&
            push!(flags, "RANK-DEFICIENT O")

        !all_finite &&
            push!(flags, "NON-FINITE")

        !LCO_converged &&
            push!(flags, "LCO RESIDUAL")

        !LCO_monotonic &&
            push!(flags, "NON-MONOTONIC LCO")

        !GD_converged &&
            push!(flags, "GD INCOMPLETE")

        flag_text =
            isempty(flags) ?
            "OK" :
            join(flags, "; ")

        println(@sprintf(
            "m=%d | rank=%d/%d cond=%.3e | cost TOD=%.8g LCO=%.8g GD=%.8g | gap=%.4g%% | runtime(s) TOD=%.3g LCO=%.3g GD=%.3g | %s",
            m_task,
            observability.numerical_rank,
            m_task,
            observability.condition_number,
            cost_TOD,
            cost_LCO,
            cost_GD,
            relative_LCO_GD_gap,
            runtime_TOD,
            runtime_LCO,
            runtime_GD,
            flag_text
        ))

    catch error_object

        if error_object isa InterruptException
            save_progress()
            println(
                "\nInterrupted during m = ",
                m_task,
                ". Previously completed dimensions remain saved."
            )
            rethrow(error_object)
        end

        error_message =
            sprint(showerror, error_object)

        @warn(
            "Dimension m = $m_task failed: $error_message"
        )

        push!(
            df_results,
            (
                m_task,
                m_task,
                1,
                T,
                seed,
                false,
                error_message,
                0,
                Inf,
                NaN,
                NaN,
                NaN,
                NaN,
                NaN,
                NaN,
                NaN,
                NaN,
                NaN,
                NaN,
                0,
                0,
                NaN,
                false,
                NaN,
                false,
                "FAILED",
                NaN,
                0,
                false,
                false,
                false
            )
        )

        sort!(df_results, :m)
        save_progress()

    end

end


save_progress()
return JLD2.load(results_file)
end
