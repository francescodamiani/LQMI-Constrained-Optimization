# Under-observed arm: retain the existing stabilization plant and only its elbow sensor.
using LinearAlgebra, Random, Printf, DataFrames, JLD2
include(joinpath(@__DIR__,"..","arm_task.jl"))

function underobserved_arm_problem(;T=300)
    p = arm_problem(;T)
    dimension_of_observation = 1
    H = zeros(dimension_of_observation,p.dimension_of_state)
    H[1,2] = 1.0
    D = p.σ_ρ_mult_sensory_noise .* H
    Ω_ω = Diagonal(fill(p.σ_ω_add_sensory_noise^2,dimension_of_observation))
    return merge(p,(;dimension_of_observation,H,D,Ω_ω))
end

function simulate_underobserved_arm(output_directory;validation_only=false)
    p = underobserved_arm_problem(T=validation_only ? 12 : 300)
    (;dimension_of_state,dimension_of_control,dimension_of_observation,T,
      A,B,H,C,D,Ω_ξ,Ω_ω,Q_matrix,R_matrix,Σ_1_x,Σ_1_z,x_1_mean,z_1_mean,
      x_state_mean_0) = p
# OPTIMIZATION PARAMETERS ----------------------------------------------------------------------------------------------

# These parameters match the redundant-arm internal-noise sweep.
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

use_GD_refinement = false
use_GD_from_IC = false
N_iter_GD_refine = 0
N_iter_GD_from_IC = 0
verbose = false

# These tolerances are diagnostics only; the LCO run uses 500 fixed passes.
lco_residual_tolerance = 1e-6
lco_monotonicity_relative_tolerance = 1e-10

sigma_eta_list = [
    0.05,
    0.10,
    0.15,
    0.20,
    0.30,
    0.40
]

validation_only && (sigma_eta_list = [0.1])

# OBSERVABILITY --------------------------------------------------------------------------------------------------------

observability =
    observability_diagnostics(
        A,
        H,
        T
    )

println(
    "Observability rank: ",
    observability.numerical_rank,
    "/",
    dimension_of_state
)

println(
    "Observability singular values: ",
    observability.singular_values
)

println(
    "Observability condition number: ",
    observability.condition_number
)

if observability.numerical_rank < dimension_of_state
    @warn(
        "The finite-horizon plant is rank-deficient from the elbow-only observation. "*
        "The plant will be retained without resampling or modification."
    )
end

# OUTPUT AND RESUMABLE STORAGE ------------------------------------------------------------------------------------------

full_folder_save_results =
    joinpath(output_directory, "")

mkpath(full_folder_save_results)

results_file = joinpath(output_directory, "arm_redundant_task_severely_underobserved_results.jld2")

df_results = DataFrame(
    sigma_eta = Float64[],
    completed = Bool[],
    failure_message = String[],
    observability_rank = Int[],
    observability_condition_number = Float64[],
    observability_nonzero_condition_number = Float64[],
    observability_rank_tolerance = Float64[],
    observability_singular_value_1 = Float64[],
    observability_singular_value_2 = Float64[],
    observability_singular_value_3 = Float64[],
    observability_singular_value_4 = Float64[],
    observability_singular_value_5 = Float64[],
    observability_singular_value_6 = Float64[],
    cost_TOD = Float64[],
    cost_LCO = Float64[],
    relative_LCO_improvement_over_TOD_percent = Float64[],
    runtime_TOD_s = Float64[],
    runtime_LCO_s = Float64[],
    LCO_updates = Int[],
    LCO_passes = Int[],
    LCO_convergence_residual = Float64[],
    LCO_converged = Bool[],
    LCO_maximum_positive_cost_increment = Float64[],
    LCO_monotonic = Bool[],
    LCO_selected_method = String[],
    LCO_selected_rho = Float64[],
    numerical_warnings = String[],
    all_costs_and_runtimes_finite = Bool[]
)

lco_cost_histories =
    Dict{Float64,Vector{Float64}}()

K_TOD_all =
    Dict{Float64,Array{Float64,3}}()

L_TOD_all =
    Dict{Float64,Array{Float64,3}}()

K_LCO_all =
    Dict{Float64,Array{Float64,3}}()

L_LCO_all =
    Dict{Float64,Array{Float64,3}}()

if isfile(results_file)

    saved_results =
        JLD2.load(results_file)
    get(saved_results,"validation_only",false) == validation_only || error("Cannot mix validation and publication runs")

    df_results =
        saved_results["df_results"]

    lco_cost_histories =
        saved_results["lco_cost_histories"]

    K_TOD_all =
        saved_results["K_TOD_all"]

    L_TOD_all =
        saved_results["L_TOD_all"]

    K_LCO_all =
        saved_results["K_LCO_all"]

    L_LCO_all =
        saved_results["L_LCO_all"]

    println(
        "Loaded ",
        nrow(df_results),
        " saved noise-level rows."
    )

    # Retain only the noise levels used in the current nonzero sweep.
    filter!(
        row ->
            any(
                σ_η ->
                    isapprox(
                        row.sigma_eta,
                        σ_η;
                        atol=1e-12,
                        rtol=0.0
                    ),
                sigma_eta_list
            ),
        df_results
    )

    for solution_dictionary in (
        lco_cost_histories,
        K_TOD_all,
        L_TOD_all,
        K_LCO_all,
        L_LCO_all
    )
        filter!(
            noise_solution ->
                any(
                    σ_η ->
                        isapprox(
                            first(noise_solution),
                            σ_η;
                            atol=1e-12,
                            rtol=0.0
                        ),
                    sigma_eta_list
                ),
            solution_dictionary
        )
    end

end

function save_progress()

    JLD2.jldsave(
        results_file;
        generated_by_current_core=true,
        validation_only=validation_only,
        df_results=df_results,
        lco_cost_histories=lco_cost_histories,
        K_TOD_all=K_TOD_all,
        L_TOD_all=L_TOD_all,
        K_LCO_all=K_LCO_all,
        L_LCO_all=L_LCO_all,
        sigma_eta_list=sigma_eta_list,
        A=A,
        B=B,
        H=H,
        C=C,
        D=D,
        Ω_ξ=Ω_ξ,
        Ω_ω=Ω_ω,
        Q_matrix=Q_matrix,
        R_matrix=R_matrix,
        Σ_1_x=Σ_1_x,
        Σ_1_z=Σ_1_z,
        x_1_mean=x_1_mean,
        z_1_mean=z_1_mean,
        T=T,
        N_iter_Todorov=N_iter_Todorov,
        N_iter_tot=N_iter_tot,
        iters_K=iters_K,
        iters_L=iters_L,
        ρ_list=ρ_list,
        ρ_list_large=ρ_list_large,
        lco_residual_tolerance=
            lco_residual_tolerance,
        lco_monotonicity_relative_tolerance=
            lco_monotonicity_relative_tolerance,
        observability=observability
    )



end

# RUN ------------------------------------------------------------------------------------------------------------------

println(
    "\n================ SEVERELY UNDER-OBSERVED ARM BENCHMARK ================\n"
)

println("sigma_eta values: ", sigma_eta_list)
println("m = 6, p = 9, k = 1, T = 300")
println("H = [0, 1, 0, 0, 0, 0]")

for σ_η_internal_noise in sigma_eta_list

    completed_rows =
        findall(
            isapprox.(
                df_results.sigma_eta,
                σ_η_internal_noise;
                atol=1e-12,
                rtol=0.0
            ) .&
            df_results.completed
        )

    if !isempty(completed_rows)
        println(
            "\nsigma_eta = ",
            σ_η_internal_noise,
            " already completed."
        )
        continue
    end

    if nrow(df_results) > 0
        filter!(
            row ->
                !isapprox(
                    row.sigma_eta,
                    σ_η_internal_noise;
                    atol=1e-12,
                    rtol=0.0
                ),
            df_results
        )
    end

    println(
        "\n------------------- sigma_eta = ",
        σ_η_internal_noise,
        " -------------------\n"
    )

    Ω_η =
        Diagonal(
            fill(
                σ_η_internal_noise^2,
                dimension_of_state
            )
        )

    try

        local L_TOD
        local K_TOD

        runtime_TOD =
            @elapsed begin
                L_TOD, K_TOD =
                    Todorov_optimization_multidimensional_case(
                        A, B, H, C, D, T,
                        dimension_of_state,
                        dimension_of_control,
                        dimension_of_observation,
                        x_state_mean_0,
                        Σ_1_x,
                        Σ_1_z,
                        Q_matrix,
                        R_matrix,
                        Ω_ξ,
                        Ω_ω,
                        Ω_η,
                        N_iter_Todorov
                    )
            end

        cost_TOD =
            expected_cost_using_trace_formula_induction(
                L_TOD,
                K_TOD,
                A, B, H, C, D,
                T,
                dimension_of_state,
                x_1_mean,
                z_1_mean,
                Σ_1_x,
                Σ_1_z,
                Q_matrix,
                R_matrix,
                Ω_ξ,
                Ω_ω,
                Ω_η
            )

        K_initial = copy(K_TOD)
        L_initial = copy(L_TOD)

        local result_LCO

        runtime_LCO =
            @elapsed begin
                result_LCO =
                    Optimal_Control_Estimation_FullProcedure(
                        K_initial,
                        L_initial,
                        A, B, H, C, D, T,
                        Ω_ξ,
                        Ω_ω,
                        Ω_η,
                        Q_matrix,
                        R_matrix,
                        Σ_1_x,
                        Σ_1_z,
                        x_1_mean,
                        z_1_mean;
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
                L_LCO,
                K_LCO,
                A, B, H, C, D,
                T,
                dimension_of_state,
                x_1_mean,
                z_1_mean,
                Σ_1_x,
                Σ_1_z,
                Q_matrix,
                R_matrix,
                Ω_ξ,
                Ω_ω,
                Ω_η
            )

        lco_cost_history =
            collect(
                result_LCO.best_cost_list
            )

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
            max(
                0.0,
                maximum(cost_increments)
            )

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
            isfinite(
                LCO_convergence_residual
            ) &&
            LCO_convergence_residual <=
            lco_residual_tolerance

        relative_LCO_improvement =
            100*(cost_TOD-cost_LCO)/
            max(
                abs(cost_TOD),
                eps(Float64)
            )

        all_finite =
            all(
                isfinite,
                [
                    cost_TOD,
                    cost_LCO,
                    runtime_TOD,
                    runtime_LCO,
                    LCO_convergence_residual,
                    maximum_positive_increment
                ]
            )

        numerical_warning_list =
            String[]

        observability.numerical_rank <
        dimension_of_state &&
            push!(
                numerical_warning_list,
                "rank-deficient observability matrix"
            )

        !all_finite &&
            push!(
                numerical_warning_list,
                "non-finite cost or diagnostic"
            )

        !LCO_converged &&
            push!(
                numerical_warning_list,
                "LCO final residual exceeds tolerance"
            )

        !LCO_monotonic &&
            push!(
                numerical_warning_list,
                "positive LCO cost increment exceeds numerical tolerance"
            )

        occursin(
            "GD",
            result_LCO.best_method
        ) &&
            push!(
                numerical_warning_list,
                "non-regularized LCO fallback selected"
            )

        numerical_warnings =
            isempty(numerical_warning_list) ?
            "" :
            join(
                numerical_warning_list,
                "; "
            )

        lco_cost_histories[
            σ_η_internal_noise
        ] =
            lco_cost_history

        K_TOD_all[
            σ_η_internal_noise
        ] =
            K_TOD

        L_TOD_all[
            σ_η_internal_noise
        ] =
            L_TOD

        K_LCO_all[
            σ_η_internal_noise
        ] =
            K_LCO

        L_LCO_all[
            σ_η_internal_noise
        ] =
            L_LCO

        singular_values =
            observability.singular_values

        push!(
            df_results,
            (
                sigma_eta=
                    σ_η_internal_noise,
                completed=true,
                failure_message="",
                observability_rank=
                    observability.numerical_rank,
                observability_condition_number=
                    observability.condition_number,
                observability_nonzero_condition_number=
                    observability.nonzero_condition_number,
                observability_rank_tolerance=
                    observability.rank_tolerance,
                observability_singular_value_1=
                    singular_values[1],
                observability_singular_value_2=
                    singular_values[2],
                observability_singular_value_3=
                    singular_values[3],
                observability_singular_value_4=
                    singular_values[4],
                observability_singular_value_5=
                    singular_values[5],
                observability_singular_value_6=
                    singular_values[6],
                cost_TOD=cost_TOD,
                cost_LCO=cost_LCO,
                relative_LCO_improvement_over_TOD_percent=
                    relative_LCO_improvement,
                runtime_TOD_s=runtime_TOD,
                runtime_LCO_s=runtime_LCO,
                LCO_updates=LCO_updates,
                LCO_passes=LCO_passes,
                LCO_convergence_residual=
                    LCO_convergence_residual,
                LCO_converged=LCO_converged,
                LCO_maximum_positive_cost_increment=
                    maximum_positive_increment,
                LCO_monotonic=LCO_monotonic,
                LCO_selected_method=
                    result_LCO.best_method,
                LCO_selected_rho=
                    result_LCO.best_ρ,
                numerical_warnings=
                    numerical_warnings,
                all_costs_and_runtimes_finite=
                    all_finite
            )
        )

        sort!(
            df_results,
            :sigma_eta
        )

        save_progress()

        flags = String[]

        observability.numerical_rank <
        dimension_of_state &&
            push!(flags, "RANK-DEFICIENT O")

        !all_finite &&
            push!(flags, "NON-FINITE")

        !LCO_converged &&
            push!(flags, "FAILED CONVERGENCE")

        !LCO_monotonic &&
            push!(flags, "NON-MONOTONIC")

        occursin(
            "GD",
            result_LCO.best_method
        ) &&
            push!(flags, "LCO FALLBACK")

        flag_text =
            isempty(flags) ?
            "OK" :
            join(flags, "; ")

        println(
            @sprintf(
                "sigma_eta=%.2f | J_TOD=%.8g J_LCO=%.8g | improvement=%.5g%% | LCO time=%.3g s | passes=%d residual=%.3e | rho=%.3g | %s",
                σ_η_internal_noise,
                cost_TOD,
                cost_LCO,
                relative_LCO_improvement,
                runtime_LCO,
                LCO_passes,
                LCO_convergence_residual,
                result_LCO.best_ρ,
                flag_text
            )
        )

    catch error_object

        if error_object isa InterruptException
            save_progress()
            println(
                "\nInterrupted during sigma_eta = ",
                σ_η_internal_noise,
                ". Completed noise levels remain saved."
            )
            rethrow(error_object)
        end

        error_message =
            sprint(
                showerror,
                error_object
            )

        @warn(
            "sigma_eta = $σ_η_internal_noise failed: $error_message"
        )

        singular_values =
            observability.singular_values

        push!(
            df_results,
            (
                sigma_eta=
                    σ_η_internal_noise,
                completed=false,
                failure_message=
                    error_message,
                observability_rank=
                    observability.numerical_rank,
                observability_condition_number=
                    observability.condition_number,
                observability_nonzero_condition_number=
                    observability.nonzero_condition_number,
                observability_rank_tolerance=
                    observability.rank_tolerance,
                observability_singular_value_1=
                    singular_values[1],
                observability_singular_value_2=
                    singular_values[2],
                observability_singular_value_3=
                    singular_values[3],
                observability_singular_value_4=
                    singular_values[4],
                observability_singular_value_5=
                    singular_values[5],
                observability_singular_value_6=
                    singular_values[6],
                cost_TOD=NaN,
                cost_LCO=NaN,
                relative_LCO_improvement_over_TOD_percent=
                    NaN,
                runtime_TOD_s=NaN,
                runtime_LCO_s=NaN,
                LCO_updates=0,
                LCO_passes=0,
                LCO_convergence_residual=NaN,
                LCO_converged=false,
                LCO_maximum_positive_cost_increment=
                    NaN,
                LCO_monotonic=false,
                LCO_selected_method="FAILED",
                LCO_selected_rho=NaN,
                numerical_warnings=
                    error_message,
                all_costs_and_runtimes_finite=false
            )
        )

        sort!(
            df_results,
            :sigma_eta
        )

        save_progress()

    end

end


save_progress()
return JLD2.load(results_file)
end
