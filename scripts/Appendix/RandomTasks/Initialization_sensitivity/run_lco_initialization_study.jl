# Initialization study, final replication settings: task seed 4321, 500 complete passes.
using LinearAlgebra, Random, Statistics, Printf, DataFrames, JLD2

function generate_stable_matrix(
    dimension;
    ρ_max=0.9,
    rng=Random.default_rng()
)

    random_matrix =
        randn(
            rng,
            dimension,
            dimension
        )

    spectral_radius =
        maximum(
            abs,
            eigvals(random_matrix)
        )

    return random_matrix .* (
        ρ_max/
        spectral_radius
    )

end

function build_random_square_task(
    dimension;
    T=10,
    ρ_max=0.9,
    σ_ξ=0.2,
    σ_ω=0.2,
    σ_ρ=0.2,
    σ_ϵ=0.2,
    σ_η=0.5,
    rng=Random.default_rng()
)

    dimension_of_state = dimension
    dimension_of_control = dimension
    dimension_of_observation = dimension

    A =
        generate_stable_matrix(
            dimension;
            ρ_max=ρ_max,
            rng=rng
        )

    B =
        generate_stable_matrix(
            dimension;
            ρ_max=ρ_max,
            rng=rng
        )

    H =
        diagm(
            ones(dimension)
        )

    D_temporary =
        abs.(
            generate_stable_matrix(
                dimension;
                ρ_max=ρ_max,
                rng=rng
            )
        )

    D =
        σ_ρ .* D_temporary

    C =
        σ_ϵ .* B

    Ω_ω =
        Diagonal(
            fill(
                σ_ω^2,
                dimension_of_observation
            )
        )

    Ω_ξ =
        Diagonal(
            fill(
                σ_ξ^2,
                dimension_of_state
            )
        )

    Ω_η =
        Diagonal(
            fill(
                σ_η^2,
                dimension_of_state
            )
        )

    x0 =
        ones(
            dimension_of_state
        )

    z0 =
        copy(x0)

    x_1_mean =
        copy(x0)

    z_1_mean =
        copy(z0)

    Σ_1_x =
        zeros(
            dimension_of_state,
            dimension_of_state
        )

    Σ_1_z =
        zeros(
            dimension_of_state,
            dimension_of_state
        )

    Q_matrix =
        zeros(
            dimension_of_state,
            dimension_of_state,
            T
        )

    R_matrix =
        zeros(
            dimension_of_control,
            dimension_of_control,
            T
        )

    for t in 1:T
        Q_matrix[:,:,t] .=
            I(dimension_of_state)

        R_matrix[:,:,t] .=
            I(dimension_of_control)
    end

    R_matrix[:,:,T] .=
        zeros(
            dimension_of_control,
            dimension_of_control
        )

    Q_matrix[:,:,1] .=
        zeros(
            dimension_of_state,
            dimension_of_state
        )

    return (
        dimension_of_state,
        dimension_of_control,
        dimension_of_observation,
        T,
        A,
        B,
        H,
        C,
        D,
        Ω_ξ,
        Ω_ω,
        Ω_η,
        x0,
        z0,
        x_1_mean,
        z_1_mean,
        Σ_1_x,
        Σ_1_z,
        Q_matrix,
        R_matrix
    )

end

# GAIN AND COST HELPERS ------------------------------------------------------------------------------------------------

function pack_LK(
    L::Array{Float64,3},
    K::Array{Float64,3}
)
    return vcat(
        vec(L),
        vec(K)
    )
end

function unpack_LK(
    theta::Vector{Float64},
    dimension_of_control,
    dimension_of_state,
    dimension_of_observation,
    T
)

    number_of_L_values =
        dimension_of_control*
        dimension_of_state*
        (T-1)

    L =
        reshape(
            theta[
                1:number_of_L_values
            ],
            dimension_of_control,
            dimension_of_state,
            T-1
        )

    K =
        reshape(
            theta[
                number_of_L_values+1:end
            ],
            dimension_of_state,
            dimension_of_observation,
            T-1
        )

    return L, K

end

function analytical_cost(
    L,
    K,
    task
)

    return expected_cost_using_trace_formula_induction(
        L,
        K,
        task.A,
        task.B,
        task.H,
        task.C,
        task.D,
        task.T,
        task.dimension_of_state,
        task.x_1_mean,
        task.z_1_mean,
        task.Σ_1_x,
        task.Σ_1_z,
        task.Q_matrix,
        task.R_matrix,
        task.Ω_ξ,
        task.Ω_ω,
        task.Ω_η
    )

end

# COMPLETE LCO PASSES USING THE EXISTING ONE-PASS SOLVER -------------------------------------------------------------

function run_one_regularization(
    K_initial,
    L_initial,
    task;
    ρ,
    tolerance,
    maximum_iterations
)

    K_current =
        copy(K_initial)

    L_current =
        copy(L_initial)

    initial_cost =
        analytical_cost(
            L_current,
            K_current,
            task
        )

    cost_trajectory =
        Float64[
            initial_cost
        ]

    if !isfinite(initial_cost)
        return (
            K=K_current,
            L=L_current,
            cost_trajectory=
                cost_trajectory,
            iterations=0,
            converged=false,
            finite=false,
            reason=
                "non-finite initial analytical cost",
            ρ=ρ
        )
    end

    # The reported replication ran all 500 passes. Keep the stopping tolerance
    # as a diagnostic so the current core's best-checkpoint return cannot end
    # the requested fixed-length experiment at a temporary plateau.
    tolerance_reached = false
    for iteration in 1:maximum_iterations

        local K_updated
        local L_updated

        try

            K_updated,
            L_updated,
            _ =
                Optimal_Control_Estimation_with_Regularized_Inverse(
                    K_current,
                    L_current,
                    task.A,
                    task.B,
                    task.H,
                    task.C,
                    task.D,
                    task.T,
                    task.Ω_ξ,
                    task.Ω_ω,
                    task.Ω_η,
                    task.Q_matrix,
                    task.R_matrix,
                    task.Σ_1_x,
                    task.Σ_1_z,
                    task.x_1_mean,
                    task.z_1_mean;
                    N_iter_tot=1,
                    ρ=ρ,
                    iters_L=1,
                    iters_K=1
                )

        catch error_object

            return (
                K=K_current,
                L=L_current,
                cost_trajectory=
                    cost_trajectory,
                iterations=
                    iteration-1,
                converged=false,
                finite=true,
                reason=
                    "solver exception: "*
                    sprint(
                        showerror,
                        error_object
                    ),
                ρ=ρ
            )

        end

        K_current .=
            K_updated

        L_current .=
            L_updated

        local current_cost

        try

            current_cost =
                analytical_cost(
                    L_current,
                    K_current,
                    task
                )

        catch error_object

            return (
                K=K_current,
                L=L_current,
                cost_trajectory=
                    cost_trajectory,
                iterations=iteration,
                converged=false,
                finite=false,
                reason=
                    "analytical-cost exception at complete iteration $iteration: "*
                    sprint(
                        showerror,
                        error_object
                    ),
                ρ=ρ
            )

        end

        push!(
            cost_trajectory,
            current_cost
        )

        if !isfinite(current_cost)
            return (
                K=K_current,
                L=L_current,
                cost_trajectory=
                    cost_trajectory,
                iterations=iteration,
                converged=false,
                finite=false,
                reason=
                    "non-finite analytical cost at complete iteration $iteration",
                ρ=ρ
            )
        end

        relative_cost_change =
            abs(
                current_cost -
                cost_trajectory[end-1]
            )/
            max(
                1.0,
                abs(
                    cost_trajectory[end-1]
                )
            )

        if relative_cost_change <
           tolerance

            tolerance_reached = true

        end

    end

    return (
        K=K_current,
        L=L_current,
        cost_trajectory=
            cost_trajectory,
        iterations=
            maximum_iterations,
        converged=tolerance_reached,
        finite=true,
        reason=
            tolerance_reached ? "" : "relative cost-change tolerance not reached before maximum iterations",
        ρ=ρ
    )

end

function run_lco_with_repository_regularization(
    K_initial,
    L_initial,
    task;
    tolerance,
    maximum_iterations
)

    initial_cost =
        analytical_cost(
            L_initial,
            K_initial,
            task
        )

    regularization_groups = [
        [0.0],
        [1e-10, 1e-8, 1e-6]
    ]

    attempted_results =
        Any[]

    for regularization_group in regularization_groups

        group_results =
            Any[]

        for ρ in regularization_group

            result =
                run_one_regularization(
                    K_initial,
                    L_initial,
                    task;
                    ρ=ρ,
                    tolerance=tolerance,
                    maximum_iterations=
                        maximum_iterations
                )

            push!(
                attempted_results,
                result
            )

            push!(
                group_results,
                result
            )

        end

        valid_results =
            filter(
                result ->
                    result.finite &&
                    isfinite(
                        result.cost_trajectory[end]
                    ) &&
                    result.cost_trajectory[end] <=
                    initial_cost,
                group_results
            )

        if !isempty(valid_results)
            selected_index =
                argmin([
                    result.cost_trajectory[end]
                    for result in valid_results
                ])

            return valid_results[
                selected_index
            ]
        end

    end

    finite_results =
        filter(
            result ->
                result.finite &&
                isfinite(
                    result.cost_trajectory[end]
                ),
            attempted_results
        )

    if isempty(finite_results)
        last_result =
            attempted_results[end]

        failure_reasons =
            join(
                unique([
                    result.reason
                    for result in attempted_results
                    if !isempty(result.reason)
                ]),
                "; "
            )

        return merge(
            last_result,
            (
                converged=false,
                reason=
                    isempty(failure_reasons) ?
                    "all repository regularization values failed" :
                    failure_reasons
            )
        )
    end

    selected_index =
        argmin([
            result.cost_trajectory[end]
            for result in finite_results
        ])

    selected_result =
        finite_results[
            selected_index
        ]

    return merge(
        selected_result,
        (
            converged=false,
            reason=
                "all repository regularization values increased the initial cost"
        )
    )

end

function build_initializations(theta_TOD;number_directions=50,large_scale=10.0,direction_seed=4322)
    norm_theta_TOD = norm(theta_TOD)
initializations =
    NamedTuple[]

push!(
    initializations,
    (
        family="TOD",
        direction_index=0,
        sign=0,
        scale=1.0,
        theta=copy(theta_TOD)
    )
)

push!(
    initializations,
    (
        family="zero",
        direction_index=0,
        sign=0,
        scale=0.0,
        theta=zeros(length(theta_TOD))
    )
)

direction_rng =
    MersenneTwister(
        direction_seed
    )

directions =
    zeros(
        length(theta_TOD),
        number_directions
    )

for direction_index in 1:number_directions

    direction =
        randn(
            direction_rng,
            length(theta_TOD)
        )

    direction ./=
        norm(direction)

    directions[:,direction_index] .=
        direction

    for direction_sign in (-1, 1)

        push!(
            initializations,
            (
                family="random_unit_scale",
                direction_index=
                    direction_index,
                sign=direction_sign,
                scale=1.0,
                theta=
                    direction_sign*
                    norm_theta_TOD*
                    direction
            )
        )

        push!(
            initializations,
            (
                family="random_large_scale",
                direction_index=
                    direction_index,
                sign=direction_sign,
                scale=
                    large_scale,
                theta=
                    direction_sign*
                    large_scale*
                    norm_theta_TOD*
                    direction
            )
        )

    end

end


return (;initializations,directions)
end

function simulate_initialization_study(output_directory;validation_only=false)
    task_seed = 4321
    number_directions = validation_only ? 1 : 50
    large_scale = 10.0
    tolerance = 1e-10
    maximum_iterations = validation_only ? 3 : 500
dimension_of_state = 5
dimension_of_control = 5
dimension_of_observation = 5
T = 10

ρ_max = 0.9
σ_ξ = 0.2
σ_ω = 0.2
σ_ρ = 0.2
σ_ϵ = 0.2
σ_η = 0.5

N_iter_Todorov = validation_only ? 3 : 500

task_rng =
    MersenneTwister(
        task_seed
    )

(
    dimension_of_state,
    dimension_of_control,
    dimension_of_observation,
    T,
    A,
    B,
    H,
    C,
    D,
    Ω_ξ,
    Ω_ω,
    Ω_η,
    x0,
    z0,
    x_1_mean,
    z_1_mean,
    Σ_1_x,
    Σ_1_z,
    Q_matrix,
    R_matrix
) =
    build_random_square_task(
        dimension_of_state;
        T=T,
        ρ_max=ρ_max,
        σ_ξ=σ_ξ,
        σ_ω=σ_ω,
        σ_ρ=σ_ρ,
        σ_ϵ=σ_ϵ,
        σ_η=σ_η,
        rng=task_rng
    )

task = (
    dimension_of_state=
        dimension_of_state,
    dimension_of_control=
        dimension_of_control,
    dimension_of_observation=
        dimension_of_observation,
    T=T,
    A=A,
    B=B,
    H=H,
    C=C,
    D=D,
    Ω_ξ=Ω_ξ,
    Ω_ω=Ω_ω,
    Ω_η=Ω_η,
    x0=x0,
    z0=z0,
    x_1_mean=x_1_mean,
    z_1_mean=z_1_mean,
    Σ_1_x=Σ_1_x,
    Σ_1_z=Σ_1_z,
    Q_matrix=Q_matrix,
    R_matrix=R_matrix
)

# TOD ------------------------------------------------------------------------------------------------------------------

println(
    "\n================ LCO INITIALIZATION STUDY ================\n"
)

println(
    "Task seed: ",
    task_seed
)

println(
    "m = p = k = 5, T = 10"
)

println(
    "Directions: ",
    number_directions,
    " | large scale: ",
    large_scale
)

println(
    "Tolerance: ",
    tolerance,
    " | maximum complete iterations: ",
    maximum_iterations
)

local L_TOD
local K_TOD

TOD_runtime =
    @elapsed begin
        L_TOD,
        K_TOD =
            Todorov_optimization_multidimensional_case(
                A,
                B,
                H,
                C,
                D,
                T,
                dimension_of_state,
                dimension_of_control,
                dimension_of_observation,
                x0,
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
    analytical_cost(
        L_TOD,
        K_TOD,
        task
    )

theta_TOD =
    pack_LK(
        L_TOD,
        K_TOD
    )

norm_theta_TOD =
    norm(theta_TOD)

isfinite(cost_TOD) ||
    error(
        "TOD produced a non-finite analytical cost."
    )

norm_theta_TOD > 0 ||
    error(
        "TOD gain vector has zero norm."
    )

println(
    @sprintf(
        "TOD cost = %.10g | runtime = %.3g s | ||theta_TOD|| = %.6g",
        cost_TOD,
        TOD_runtime,
        norm_theta_TOD
    )
)

# INITIALIZATION LIST ---------------------------------------------------------------------------------------------------

(;initializations,directions) = build_initializations(theta_TOD;number_directions,large_scale,direction_seed=task_seed+1)

number_of_runs =
    length(initializations)

# RESULT CONTAINERS ----------------------------------------------------------------------------------------------------

df_results = DataFrame(
    run_id = Int[],
    task_seed = Int[],
    initialization_family = String[],
    direction_index = Int[],
    sign = Int[],
    scale = Float64[],
    status = String[],
    failure_reason = String[],
    tolerance_reached = Bool[],
    selected_rho = Float64[],
    iterations = Int[],
    theta0_norm_over_theta_TOD_norm = Float64[],
    initial_analytical_cost = Float64[],
    final_analytical_cost = Float64[],
    C_initial_over_C_TOD = Float64[],
    C_final_over_C_TOD = Float64[],
    improvement_vs_TOD_percent = Float64[],
    gap_vs_LCO_from_TOD_percent = Float64[],
    runtime_s = Float64[]
)

cost_trajectories =
    Vector{Vector{Float64}}(
        undef,
        number_of_runs
    )

final_L_all =
    fill(
        NaN,
        number_of_runs,
        size(L_TOD,1),
        size(L_TOD,2),
        size(L_TOD,3)
    )

final_K_all =
    fill(
        NaN,
        number_of_runs,
        size(K_TOD,1),
        size(K_TOD,2),
        size(K_TOD,3)
    )

cost_LCO_from_TOD =
    NaN

# RUN ALL INITIALIZATIONS -----------------------------------------------------------------------------------------------

for (
    run_id,
    initialization
) in enumerate(initializations)

    println(
        @sprintf(
            "\nRun %d/%d | family=%s direction=%d sign=%+d scale=%.4g",
            run_id,
            number_of_runs,
            initialization.family,
            initialization.direction_index,
            initialization.sign,
            initialization.scale
        )
    )

    L_initial,
    K_initial =
        unpack_LK(
            initialization.theta,
            dimension_of_control,
            dimension_of_state,
            dimension_of_observation,
            T
        )

    local result

    runtime =
        @elapsed begin
            try

                result =
                    run_lco_with_repository_regularization(
                        K_initial,
                        L_initial,
                        task;
                        tolerance=
                            tolerance,
                        maximum_iterations=
                            maximum_iterations
                    )

            catch error_object

                initial_cost_after_exception =
                    try
                        analytical_cost(
                            L_initial,
                            K_initial,
                            task
                        )
                    catch
                        NaN
                    end

                result = (
                    K=copy(K_initial),
                    L=copy(L_initial),
                    cost_trajectory=
                        Float64[
                            initial_cost_after_exception
                        ],
                    iterations=0,
                    converged=false,
                    finite=false,
                    reason=
                        "uncaught solver exception: "*
                        sprint(
                            showerror,
                            error_object
                        ),
                    ρ=NaN
                )

            end
        end

    initial_cost =
        result.cost_trajectory[1]

    final_cost =
        result.cost_trajectory[end]

    finite_decrease =
        result.finite &&
        isfinite(final_cost) &&
        isfinite(initial_cost) &&
        final_cost < initial_cost

    status =
        finite_decrease ?
        "finite_decrease" :
        "invalid"

    failure_reason =
        if finite_decrease
            ""
        elseif !isfinite(final_cost)
            isempty(result.reason) ?
            "non-finite final analytical cost" :
            result.reason
        elseif !isfinite(initial_cost)
            "non-finite initial analytical cost"
        else
            "final analytical cost did not decrease"
        end

    if run_id == 1

        cost_LCO_from_TOD =
            final_cost

        if !finite_decrease
            @warn(
                "The LCO-from-TOD reference did not produce a finite cost decrease: $failure_reason"
            )
        end

    end

    improvement_vs_TOD =
        isfinite(final_cost) ?
        100*(cost_TOD-final_cost)/
        cost_TOD :
        NaN

    gap_vs_LCO_from_TOD =
        isfinite(final_cost) &&
        isfinite(cost_LCO_from_TOD) ?
        100*(
            final_cost-
            cost_LCO_from_TOD
        )/
        cost_LCO_from_TOD :
        NaN

    cost_trajectories[run_id] =
        result.cost_trajectory

    final_L_all[
        run_id,
        Colon(),
        Colon(),
        Colon()
    ] .=
        result.L

    final_K_all[
        run_id,
        Colon(),
        Colon(),
        Colon()
    ] .=
        result.K

    push!(
        df_results,
        (
            run_id=run_id,
            task_seed=
                task_seed,
            initialization_family=
                initialization.family,
            direction_index=
                initialization.direction_index,
            sign=
                initialization.sign,
            scale=
                initialization.scale,
            status=status,
            failure_reason=
                failure_reason,
            tolerance_reached=
                result.converged,
            selected_rho=
                result.ρ,
            iterations=
                result.iterations,
            theta0_norm_over_theta_TOD_norm=
                norm(
                    initialization.theta
                )/
                norm_theta_TOD,
            initial_analytical_cost=
                initial_cost,
            final_analytical_cost=
                final_cost,
            C_initial_over_C_TOD=
                initial_cost/
                cost_TOD,
            C_final_over_C_TOD=
                final_cost/
                cost_TOD,
            improvement_vs_TOD_percent=
                improvement_vs_TOD,
            gap_vs_LCO_from_TOD_percent=
                gap_vs_LCO_from_TOD,
            runtime_s=runtime
        )
    )

    println(
        @sprintf(
            "status=%s | iterations=%d | initial=%.8g final=%.8g | improvement=%.5g%% | gap=%.5g%% | rho=%.3g | %.3g s",
            status,
            result.iterations,
            initial_cost,
            final_cost,
            improvement_vs_TOD,
            gap_vs_LCO_from_TOD,
            result.ρ,
            runtime
        )
    )

    if status == "invalid"
        println(
            "reason: ",
            failure_reason
        )
    end

end


jld2_file = joinpath(output_directory,"lco_initialization_study_results.jld2")
task_parameters = (
    dimension_of_state=
        dimension_of_state,
    dimension_of_control=
        dimension_of_control,
    dimension_of_observation=
        dimension_of_observation,
    T=T,
    ρ_max=ρ_max,
    σ_ξ=σ_ξ,
    σ_ω=σ_ω,
    σ_ρ=σ_ρ,
    σ_ϵ=σ_ϵ,
    σ_η=σ_η,
    N_iter_Todorov=
        N_iter_Todorov
)

study_metadata = (
    task_seed=
        task_seed,
    direction_seed=
        task_seed+1,
    number_directions=
        number_directions,
    both_direction_signs=true,
    small_scale=1.0,
    large_scale=
        large_scale,
    tolerance=
        tolerance,
    maximum_iterations=
        maximum_iterations,
    complete_iteration_costs=true,
    trajectory_includes_iteration_zero=true,
    cost_TOD=cost_TOD,
    cost_LCO_from_TOD=
        cost_LCO_from_TOD
)

JLD2.jldsave(
    jld2_file;
    generated_by_current_core=true,
    validation_only=validation_only,
    df_results=df_results,
    cost_trajectories=
        cost_trajectories,
    final_L_all=final_L_all,
    final_K_all=final_K_all,
    L_TOD=L_TOD,
    K_TOD=K_TOD,
    theta_TOD=theta_TOD,
    directions=directions,
    A=A,
    B=B,
    H=H,
    C=C,
    D=D,
    Ω_ξ=Ω_ξ,
    Ω_ω=Ω_ω,
    Ω_η=Ω_η,
    Q_matrix=Q_matrix,
    R_matrix=R_matrix,
    x_1_mean=x_1_mean,
    z_1_mean=z_1_mean,
    Σ_1_x=Σ_1_x,
    Σ_1_z=Σ_1_z,
    task_parameters=
        task_parameters,
    study_metadata=
        study_metadata
)


return JLD2.load(jld2_file)
end
