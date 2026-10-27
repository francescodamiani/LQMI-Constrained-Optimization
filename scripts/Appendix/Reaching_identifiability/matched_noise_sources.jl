# Migrated from run_elbow_reaching_matched_additive_noise_study.jl; numerical loops and statistics retained.
# Only --simulate loads this module. Full settings are the historical settings.
module MatchedNoiseSources
using LinearAlgebra, Statistics, Random, Printf, MultivariateStats, JLD2
using ..IdentifiabilitySimulation: Todorov_optimization_multidimensional_case,
    Optimal_Control_Estimation_FullProcedure, expected_cost_raw_moments_propagation,
    sample_identifiability, save_results

function terminal_elbow_SD(
    K_matrix, L_matrix, A, B, H, C, D, T,
    x_1_mean, z_1_mean, Σ_1_x, Σ_1_z, Ω_ξ, Ω_ω, Ω_η
)

    S_xx = Σ_1_x .+ x_1_mean*transpose(x_1_mean)
    S_zz = Σ_1_z .+ z_1_mean*transpose(z_1_mean)
    S_xz = x_1_mean*transpose(z_1_mean)
    S_zx = transpose(S_xz)

    x_mean = copy(x_1_mean)
    z_mean = copy(z_1_mean)

    for t in 2:T

        L_t = L_matrix[:,:,t-1]
        K_t = K_matrix[:,:,t-1]
        M_t = A .+ B*L_t .- K_t*H

        S_xx_new =
            A*S_xx*transpose(A) .+
            A*S_xz*transpose(B*L_t) .+
            B*L_t*S_zx*transpose(A) .+
            B*L_t*S_zz*transpose(B*L_t) .+
            C*L_t*S_zz*transpose(C*L_t) .+
            Ω_ξ

        S_zx_new =
            K_t*H*S_xx*transpose(A) .+
            K_t*H*S_xz*transpose(B*L_t) .+
            M_t*S_zx*transpose(A) .+
            M_t*S_zz*transpose(B*L_t)

        S_xz_new = transpose(S_zx_new)

        S_zz_new =
            K_t*H*S_xx*transpose(K_t*H) .+
            K_t*D*S_xx*transpose(K_t*D) .+
            K_t*H*S_xz*transpose(M_t) .+
            M_t*S_zx*transpose(K_t*H) .+
            M_t*S_zz*transpose(M_t) .+
            K_t*Ω_ω*transpose(K_t) .+
            Ω_η

        x_mean_new =
            A*x_mean .+
            B*L_t*z_mean

        z_mean_new =
            K_t*H*x_mean .+
            M_t*z_mean

        S_xx = S_xx_new
        S_zz = S_zz_new
        S_xz = S_xz_new
        S_zx = S_zx_new

        x_mean = x_mean_new
        z_mean = z_mean_new

    end

    terminal_elbow_variance =
        S_xx[2,2] - x_mean[2]^2

    if terminal_elbow_variance < -1e-10
        error("The analytical terminal elbow variance is negative.")
    end

    return rad2deg(sqrt(max(terminal_elbow_variance, 0.0)))

end

function match_noise_amplitude(noise_grid, endpoint_SD_grid, target_SD)

    exact_idx =
        argmin(abs.(endpoint_SD_grid .- target_SD))

    if isapprox(
        endpoint_SD_grid[exact_idx],
        target_SD;
        atol=1e-12,
        rtol=1e-10
    )
        return noise_grid[exact_idx]
    end

    crossing_indices = Int[]

    for idx in 1:length(noise_grid)-1
        difference_1 = endpoint_SD_grid[idx] - target_SD
        difference_2 = endpoint_SD_grid[idx+1] - target_SD

        if difference_1*difference_2 <= 0.0
            push!(crossing_indices, idx)
        end
    end

    if isempty(crossing_indices)
        closest_idx =
            argmin(abs.(endpoint_SD_grid .- target_SD))

        @warn(
            "The dense candidate grid does not bracket target terminal "*
            "SD = $target_SD degrees. Continuing with the closest "*
            "available condition: noise amplitude = "*
            "$(noise_grid[closest_idx]), analytical terminal SD = "*
            "$(endpoint_SD_grid[closest_idx]) degrees. This point is "*
            "not endpoint-variability matched."
        )

        return noise_grid[closest_idx]
    end

    crossing_scores = [
        abs(endpoint_SD_grid[idx] - target_SD) +
        abs(endpoint_SD_grid[idx+1] - target_SD)
        for idx in crossing_indices
    ]

    crossing_idx =
        crossing_indices[argmin(crossing_scores)]

    noise_1 = noise_grid[crossing_idx]
    noise_2 = noise_grid[crossing_idx+1]
    SD_1 = endpoint_SD_grid[crossing_idx]
    SD_2 = endpoint_SD_grid[crossing_idx+1]

    if isapprox(SD_1, SD_2; atol=1e-14, rtol=1e-12)
        return (noise_1 + noise_2)/2
    end

    interpolation_weight =
        (target_SD - SD_1)/(SD_2 - SD_1)

    return noise_1 +
           interpolation_weight*(noise_2 - noise_1)

end

function compute_MC_metrics(x_trajectories, time_vector)

    N_trials = size(x_trajectories,2)

    elbow_angles_deg =
        rad2deg.(x_trajectories[2,:,:])

    elbow_velocities_deg =
        rad2deg.(x_trajectories[5,:,:])

    peak_absolute_velocity = zeros(N_trials)
    peak_time = zeros(N_trials)

    for trial in 1:N_trials
        peak_absolute_velocity[trial],
        peak_idx =
            findmax(abs.(elbow_velocities_deg[trial,:]))

        peak_time[trial] =
            time_vector[peak_idx]
    end

    terminal_SD_deg =
        std(elbow_angles_deg[:,end])

    terminal_SD_percent =
        100*terminal_SD_deg/20.0

    mean_peak_velocity =
        mean(peak_absolute_velocity)

    sem_peak_velocity =
        std(peak_absolute_velocity)/sqrt(N_trials)

    mean_peak_time =
        mean(peak_time)

    sem_peak_time =
        std(peak_time)/sqrt(N_trials)

    mean_velocity_profile =
        vec(mean(elbow_velocities_deg, dims=1))

    angle_variance_profile =
        vec(var(elbow_angles_deg, dims=1))

    return terminal_SD_deg,
           terminal_SD_percent,
           mean_peak_velocity,
           sem_peak_velocity,
           mean_peak_time,
           sem_peak_time,
           mean_velocity_profile,
           angle_variance_profile

end

function compute_extended_metrics(
    x_trajectories, z_trajectories,
    K_matrix, L_matrix, time_indices
)

    n = size(x_trajectories,1)
    N_trials = size(x_trajectories,2)
    T = size(x_trajectories,3)
    m = size(L_matrix,1)

    u_trajectories =
        zeros(m, N_trials, T-1)

    for t in 1:T-1
        u_trajectories[:,:,t] .=
            L_matrix[:,:,t]*z_trajectories[:,:,t]
    end

    std_x = zeros(n)
    std_z = zeros(n)
    std_u = zeros(m)

    for component in 1:n
        std_x[component] =
            sqrt(mean([
                var(x_trajectories[component,:,t])
                for t in 1:T
            ]))

        std_z[component] =
            sqrt(mean([
                var(z_trajectories[component,:,t])
                for t in 1:T
            ]))
    end

    for component in 1:m
        std_u[component] =
            sqrt(mean([
                var(u_trajectories[component,:,t])
                for t in 1:T-1
            ]))
    end

    centered_z =
        zeros(n, N_trials, length(time_indices))

    for (time_idx, t) in enumerate(time_indices)
        centered_z[:,:,time_idx] .=
            z_trajectories[:,:,t] .-
            mean(z_trajectories[:,:,t], dims=2)
    end

    z_matrix =
        reshape(centered_z, n, :)

    pca_model = fit(PCA, z_matrix)
    eigenvalues = principalvars(pca_model)
    EV1 = eigenvalues[1]/sum(eigenvalues)

    PC1 = projection(pca_model)[:,1]
    PC1 ./= norm(PC1)

    PC1_readout_over_time =
        zeros(length(time_indices))

    for (time_idx, t) in enumerate(time_indices)
        PC1_readout_over_time[time_idx] =
            norm(L_matrix[:,:,t]*PC1)^2
    end

    PC1_readout_mean =
        mean(PC1_readout_over_time)

    PC1_readout_sem =
        std(PC1_readout_over_time)/
        sqrt(length(PC1_readout_over_time))

    K_Frobenius =
        mean([
            norm(K_matrix[:,:,t])
            for t in 2:size(K_matrix,3)-2
        ])

    L_Frobenius =
        mean([
            norm(L_matrix[:,:,t])
            for t in 2:size(L_matrix,3)-2
        ])

    velocity_norm_peaks = zeros(N_trials)

    for trial in 1:N_trials
        velocity_norm_over_time = [
            norm(x_trajectories[4:6,trial,t])
            for t in 1:T
        ]

        velocity_norm_peaks[trial] =
            maximum(velocity_norm_over_time)
    end

    mean_peak_velocity_norm =
        mean(velocity_norm_peaks)

    sem_peak_velocity_norm =
        std(velocity_norm_peaks)/sqrt(N_trials)

    elbow_endpoint_deg =
        rad2deg.(x_trajectories[2,:,T])

    endpoint_bias_deg =
        abs(mean(elbow_endpoint_deg))

    success_probability =
        mean(abs.(elbow_endpoint_deg) .< 2.0)

    terminal_absolute_velocity_deg_s =
        mean(abs.(rad2deg.(
            x_trajectories[5,:,T]
        )))

    control_effort =
        mean(vec(sum(
            abs2,
            u_trajectories,
            dims=(1,3)
        )))

    return (
        std_x=std_x,
        std_z=std_z,
        std_u=std_u,
        EV1=EV1,
        PC1=PC1,
        PC1_readout_mean=PC1_readout_mean,
        PC1_readout_sem=PC1_readout_sem,
        K_Frobenius=K_Frobenius,
        L_Frobenius=L_Frobenius,
        mean_peak_velocity_norm=mean_peak_velocity_norm,
        sem_peak_velocity_norm=sem_peak_velocity_norm,
        endpoint_bias_deg=endpoint_bias_deg,
        success_probability=success_probability,
        terminal_absolute_velocity_deg_s=
            terminal_absolute_velocity_deg_s,
        control_effort=control_effort
    )

end

function mean_log_pseudodeterminant(gain_matrix)

    log_pseudodeterminant_over_time =
        zeros(size(gain_matrix,3)-3)

    for (time_idx, t) in enumerate(
        2:size(gain_matrix,3)-2
    )
        singular_values =
            svd(gain_matrix[:,:,t]).S

        retained_singular_values =
            singular_values[
                singular_values .> 1e-6
            ]

        log_pseudodeterminant_over_time[time_idx] =
            sum(log, retained_singular_values; init=0.0)
    end

    return mean(log_pseudodeterminant_over_time)

end


function run_experiment(reference_results, output_dir; validation_only=false)
#REFERENCE REACHING SWEEP -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

println("\n------- ELBOW REACHING: ENDPOINT-VARIABILITY-MATCHED NOISE-SOURCE STUDY ----------------------\n")



Δt = reference_results["Δt"]
T_move = reference_results["T_move"]
T = reference_results["T"]
time_vector = reference_results["time_vector"]

dimension_of_state = reference_results["dimension_of_state"]
dimension_of_control = reference_results["dimension_of_control"]
dimension_of_observation = reference_results["dimension_of_observation"]

x_1_mean = reference_results["x_1_mean"]
z_1_mean = reference_results["z_1_mean"]
Σ_1_x = reference_results["Σ_1_x"]
Σ_1_z = reference_results["Σ_1_z"]

A = reference_results["A"]
B = reference_results["B"]
C = reference_results["C"]
H = reference_results["H"]
D = reference_results["D"]

Q_matrix = reference_results["Q_matrix"]
R_matrix = reference_results["R_matrix"]

Ω_ξ_baseline = reference_results["Ω_ξ"]
Ω_ω_baseline = reference_results["Ω_ω"]
Ω_η_zero = zeros(dimension_of_state, dimension_of_state)

σ_ξ_baseline = reference_results["σ_ξ"]
σ_ω_baseline =
    reference_results["σ_ω_add_sensory_noise"]

σ_ϵ_baseline =
    reference_results["σ_ϵ_mult_control_noise"]

σ_ρ_baseline =
    reference_results["σ_ρ_mult_sensory_noise"]

sigma_eta_list =
    reference_results["sigma_eta_list"]

K_LCO_internal =
    reference_results["K_LCO_all"]

L_LCO_internal =
    reference_results["L_LCO_all"]

x_store_LCO_internal =
    reference_results["x_store_LCO"]

z_store_LCO_internal =
    reference_results["z_store_LCO"]

trial_costs_LCO_internal =
    reference_results["trial_costs_LCO"]

analytical_cost_LCO_internal =
    reference_results["analytical_cost_LCO"]

N_reference = length(sigma_eta_list)
N_trials_MC = validation_only ? 64 : 2000
random_seed = reference_results["random_seed"]

baseline_noise_values = [
    σ_ξ_baseline,
    σ_ω_baseline,
    σ_ϵ_baseline,
    σ_ρ_baseline
]

if !all(isapprox.(baseline_noise_values, 0.01))
    error(
        "The saved reaching sweep does not use the required 0.01 external-noise baseline."
    )
end

if reference_results["N_trials_MC"] != N_trials_MC
    error("The saved internal-noise reference does not contain 2000 trials.")
end

#OPTIMIZATION PARAMETERS: IDENTICAL TO THE EXISTING REACHING SWEEP -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

#DENSE CANDIDATE GRIDS -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

maximum_candidate_noise = 0.10
maximum_candidate_multiplicative_noise = 0.30
N_candidate_grid = validation_only ? 3 : 31

sigma_xi_grid =
    sort(unique(vcat(
        collect(range(
            0.0,
            stop=maximum_candidate_noise,
            length=N_candidate_grid
        )),
        [σ_ξ_baseline]
    )))

#A tiny positive lower bound avoids a singular zero-sensory-noise condition.
sigma_omega_grid =
    sort(unique(vcat(
        collect(range(
            1e-4,
            stop=maximum_candidate_noise,
            length=N_candidate_grid
        )),
        [σ_ω_baseline]
    )))

sigma_epsilon_grid =
    sort(unique(vcat(
        collect(range(
            0.0,
            stop=maximum_candidate_multiplicative_noise,
            length=N_candidate_grid
        )),
        [σ_ϵ_baseline]
    )))

sigma_rho_grid =
    sort(unique(vcat(
        collect(range(
            0.0,
            stop=maximum_candidate_multiplicative_noise,
            length=N_candidate_grid
        )),
        [σ_ρ_baseline]
    )))

source_names = [
    "internal_noise",
    "additive_process_noise",
    "additive_sensory_noise",
    "multiplicative_control_noise",
    "multiplicative_sensory_noise"
]

#SEPARATE OUTPUT FOLDER -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

grid_results_file = joinpath(output_dir,"candidate_noise_grids.jld2")
matched_results_file = joinpath(output_dir,"matched_noise_sources.jld2")

#ANALYTICAL REFERENCE TERMINAL VARIABILITY -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

reference_terminal_SD_analytical_deg =
    zeros(N_reference)

for reference_idx in 1:N_reference

    Ω_η =
        Diagonal(fill(
            sigma_eta_list[reference_idx]^2,
            dimension_of_state
        ))

    reference_terminal_SD_analytical_deg[reference_idx] =
        terminal_elbow_SD(
            K_LCO_internal[reference_idx],
            L_LCO_internal[reference_idx],
            A, B, H, C, D, T,
            x_1_mean, z_1_mean,
            Σ_1_x, Σ_1_z,
            Ω_ξ_baseline,
            Ω_ω_baseline,
            Ω_η
        )

end

#DENSE GRID: REOPTIMIZE LCO FOR EACH EXTERNAL-NOISE VALUE -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

N_xi_grid = length(sigma_xi_grid)
N_omega_grid = length(sigma_omega_grid)
N_epsilon_grid = length(sigma_epsilon_grid)
N_rho_grid = length(sigma_rho_grid)

K_LCO_xi_grid =
    zeros(
        dimension_of_state,
        dimension_of_observation,
        T-1,
        N_xi_grid
    )

L_LCO_xi_grid =
    zeros(
        dimension_of_control,
        dimension_of_state,
        T-1,
        N_xi_grid
    )

K_LCO_omega_grid =
    zeros(
        dimension_of_state,
        dimension_of_observation,
        T-1,
        N_omega_grid
    )

L_LCO_omega_grid =
    zeros(
        dimension_of_control,
        dimension_of_state,
        T-1,
        N_omega_grid
    )

K_LCO_epsilon_grid =
    zeros(
        dimension_of_state,
        dimension_of_observation,
        T-1,
        N_epsilon_grid
    )

L_LCO_epsilon_grid =
    zeros(
        dimension_of_control,
        dimension_of_state,
        T-1,
        N_epsilon_grid
    )

K_LCO_rho_grid =
    zeros(
        dimension_of_state,
        dimension_of_observation,
        T-1,
        N_rho_grid
    )

L_LCO_rho_grid =
    zeros(
        dimension_of_control,
        dimension_of_state,
        T-1,
        N_rho_grid
    )

terminal_SD_xi_grid = zeros(N_xi_grid)
terminal_SD_omega_grid = zeros(N_omega_grid)
terminal_SD_epsilon_grid = zeros(N_epsilon_grid)
terminal_SD_rho_grid = zeros(N_rho_grid)

analytical_cost_xi_grid = zeros(N_xi_grid)
analytical_cost_omega_grid = zeros(N_omega_grid)
analytical_cost_epsilon_grid = zeros(N_epsilon_grid)
analytical_cost_rho_grid = zeros(N_rho_grid)

completed_xi_grid = falses(N_xi_grid)
completed_omega_grid = falses(N_omega_grid)
completed_epsilon_grid = falses(N_epsilon_grid)
completed_rho_grid = falses(N_rho_grid)

if isfile(grid_results_file)

    println("Loading resumable dense-grid results from:")
    println(grid_results_file)

    saved_grid = JLD2.load(grid_results_file)

    if saved_grid["sigma_xi_grid"] != sigma_xi_grid ||
       saved_grid["sigma_omega_grid"] != sigma_omega_grid ||
       saved_grid["sigma_epsilon_grid"] != sigma_epsilon_grid ||
       saved_grid["sigma_rho_grid"] != sigma_rho_grid
        error(
            "The saved dense grids do not match the present grids. "*
            "The existing results were not overwritten."
        )
    end

    K_LCO_xi_grid .= saved_grid["K_LCO_xi_grid"]
    L_LCO_xi_grid .= saved_grid["L_LCO_xi_grid"]
    K_LCO_omega_grid .= saved_grid["K_LCO_omega_grid"]
    L_LCO_omega_grid .= saved_grid["L_LCO_omega_grid"]
    K_LCO_epsilon_grid .= saved_grid["K_LCO_epsilon_grid"]
    L_LCO_epsilon_grid .= saved_grid["L_LCO_epsilon_grid"]
    K_LCO_rho_grid .= saved_grid["K_LCO_rho_grid"]
    L_LCO_rho_grid .= saved_grid["L_LCO_rho_grid"]

    terminal_SD_xi_grid .= saved_grid["terminal_SD_xi_grid"]
    terminal_SD_omega_grid .= saved_grid["terminal_SD_omega_grid"]
    terminal_SD_epsilon_grid .=
        saved_grid["terminal_SD_epsilon_grid"]
    terminal_SD_rho_grid .=
        saved_grid["terminal_SD_rho_grid"]

    analytical_cost_xi_grid .=
        saved_grid["analytical_cost_xi_grid"]

    analytical_cost_omega_grid .=
        saved_grid["analytical_cost_omega_grid"]
    analytical_cost_epsilon_grid .=
        saved_grid["analytical_cost_epsilon_grid"]
    analytical_cost_rho_grid .=
        saved_grid["analytical_cost_rho_grid"]

    completed_xi_grid .= saved_grid["completed_xi_grid"]
    completed_omega_grid .= saved_grid["completed_omega_grid"]
    completed_epsilon_grid .= saved_grid["completed_epsilon_grid"]
    completed_rho_grid .= saved_grid["completed_rho_grid"]

end

for source_idx in 1:4

    if source_idx == 1
        noise_grid = sigma_xi_grid
        completed_grid = completed_xi_grid
        source_label = "sigma_xi"
    elseif source_idx == 2
        noise_grid = sigma_omega_grid
        completed_grid = completed_omega_grid
        source_label = "sigma_omega"
    elseif source_idx == 3
        noise_grid = sigma_epsilon_grid
        completed_grid = completed_epsilon_grid
        source_label = "sigma_epsilon"
    else
        noise_grid = sigma_rho_grid
        completed_grid = completed_rho_grid
        source_label = "sigma_rho"
    end

    for grid_idx in eachindex(noise_grid)

        if completed_grid[grid_idx]
            continue
        end

        candidate_noise = noise_grid[grid_idx]

        C_condition = C
        D_condition = D
        Ω_ξ_condition = Ω_ξ_baseline
        Ω_ω_condition = Ω_ω_baseline

        if source_idx == 1
            Ω_ξ_condition =
                Diagonal(fill(
                    candidate_noise^2,
                    dimension_of_state
                ))
        elseif source_idx == 2
            Ω_ω_condition =
                Diagonal(fill(
                    candidate_noise^2,
                    dimension_of_observation
                ))
        elseif source_idx == 3
            C_condition =
                candidate_noise .* B
        else
            D_condition =
                candidate_noise .* H
        end

        println(
            "\nDense grid ",
            source_label,
            " ",
            grid_idx, "/", length(noise_grid),
            ": ", candidate_noise
        )

        #The baseline grid point is exactly the saved sigma_eta=0 condition.
        if candidate_noise == 0.01
            baseline_idx =
                argmin(abs.(sigma_eta_list))

            K_LCO = K_LCO_internal[baseline_idx]
            L_LCO = L_LCO_internal[baseline_idx]
        else
            L_TOD, K_TOD =
                Todorov_optimization_multidimensional_case(
                    A, B, H,
                    C_condition, D_condition, T,
                    dimension_of_state,
                    dimension_of_control,
                    dimension_of_observation,
                    x_1_mean,
                    Σ_1_x, Σ_1_z,
                    Q_matrix, R_matrix,
                    Ω_ξ_condition,
                    Ω_ω_condition,
                    Ω_η_zero,
                    N_iter_Todorov
                )

            K_LCO, L_LCO, _, _, _, _, _ =
                Optimal_Control_Estimation_FullProcedure(
                    copy(K_TOD), copy(L_TOD),
                    A, B, H,
                    C_condition, D_condition, T,
                    Ω_ξ_condition,
                    Ω_ω_condition,
                    Ω_η_zero,
                    Q_matrix, R_matrix,
                    Σ_1_x, Σ_1_z,
                    x_1_mean, z_1_mean;
                    ρ_list = ρ_list,
                    ρ_list_large = ρ_list_large,
                    N_iter_tot = N_iter_tot,
                    iters_K = iters_K,
                    iters_L = iters_L,
                    N_iter_GD_refine = 0,
                    N_iter_GD_from_IC = 0,
                    gd_method = gd_method,
                    η_GD = η_GD,
                    maxiter_GD = maxiter_GD,
                    tol_GD = tol_GD,
                    β1_Adam = β1_Adam,
                    β2_Adam = β2_Adam,
                    ϵ_Adam = ϵ_Adam,
                    use_GD_refinement = false,
                    use_GD_from_IC = false,
                    verbose = false
                )
        end

        if !all(isfinite, K_LCO) ||
           !all(isfinite, L_LCO)
            error(
                "Non-finite LCO gains at $source_label = $candidate_noise"
            )
        end

        terminal_SD =
            terminal_elbow_SD(
                K_LCO, L_LCO,
                A, B, H,
                C_condition, D_condition, T,
                x_1_mean, z_1_mean,
                Σ_1_x, Σ_1_z,
                Ω_ξ_condition,
                Ω_ω_condition,
                Ω_η_zero
            )

        analytical_cost =
            expected_cost_raw_moments_propagation(
                T,
                dimension_of_state,
                dimension_of_control,
                dimension_of_observation,
                K_LCO, L_LCO,
                A, B, H,
                C_condition, D_condition,
                x_1_mean, z_1_mean,
                Σ_1_x, Σ_1_z,
                Q_matrix, R_matrix,
                Ω_ξ_condition,
                Ω_ω_condition,
                Ω_η_zero
            )

        if source_idx == 1
            K_LCO_xi_grid[:,:,:,grid_idx] .= K_LCO
            L_LCO_xi_grid[:,:,:,grid_idx] .= L_LCO
            terminal_SD_xi_grid[grid_idx] = terminal_SD
            analytical_cost_xi_grid[grid_idx] = analytical_cost
            completed_xi_grid[grid_idx] = true
        elseif source_idx == 2
            K_LCO_omega_grid[:,:,:,grid_idx] .= K_LCO
            L_LCO_omega_grid[:,:,:,grid_idx] .= L_LCO
            terminal_SD_omega_grid[grid_idx] = terminal_SD
            analytical_cost_omega_grid[grid_idx] = analytical_cost
            completed_omega_grid[grid_idx] = true
        elseif source_idx == 3
            K_LCO_epsilon_grid[:,:,:,grid_idx] .= K_LCO
            L_LCO_epsilon_grid[:,:,:,grid_idx] .= L_LCO
            terminal_SD_epsilon_grid[grid_idx] = terminal_SD
            analytical_cost_epsilon_grid[grid_idx] = analytical_cost
            completed_epsilon_grid[grid_idx] = true
        else
            K_LCO_rho_grid[:,:,:,grid_idx] .= K_LCO
            L_LCO_rho_grid[:,:,:,grid_idx] .= L_LCO
            terminal_SD_rho_grid[grid_idx] = terminal_SD
            analytical_cost_rho_grid[grid_idx] = analytical_cost
            completed_rho_grid[grid_idx] = true
        end

        save_results(validation_only,
            grid_results_file,
            "sigma_xi_grid", sigma_xi_grid,
            "sigma_omega_grid", sigma_omega_grid,
            "sigma_epsilon_grid", sigma_epsilon_grid,
            "sigma_rho_grid", sigma_rho_grid,
            "K_LCO_xi_grid", K_LCO_xi_grid,
            "L_LCO_xi_grid", L_LCO_xi_grid,
            "K_LCO_omega_grid", K_LCO_omega_grid,
            "L_LCO_omega_grid", L_LCO_omega_grid,
            "K_LCO_epsilon_grid", K_LCO_epsilon_grid,
            "L_LCO_epsilon_grid", L_LCO_epsilon_grid,
            "K_LCO_rho_grid", K_LCO_rho_grid,
            "L_LCO_rho_grid", L_LCO_rho_grid,
            "terminal_SD_xi_grid", terminal_SD_xi_grid,
            "terminal_SD_omega_grid", terminal_SD_omega_grid,
            "terminal_SD_epsilon_grid", terminal_SD_epsilon_grid,
            "terminal_SD_rho_grid", terminal_SD_rho_grid,
            "analytical_cost_xi_grid", analytical_cost_xi_grid,
            "analytical_cost_omega_grid",
                analytical_cost_omega_grid,
            "analytical_cost_epsilon_grid",
                analytical_cost_epsilon_grid,
            "analytical_cost_rho_grid",
                analytical_cost_rho_grid,
            "completed_xi_grid", completed_xi_grid,
            "completed_omega_grid", completed_omega_grid,
            "completed_epsilon_grid", completed_epsilon_grid,
            "completed_rho_grid", completed_rho_grid
        )

    end
end

#INTERPOLATE THE MATCHED EXTERNAL-NOISE VALUES -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

matched_sigma_xi = zeros(N_reference)
matched_sigma_omega = zeros(N_reference)
matched_sigma_epsilon = zeros(N_reference)
matched_sigma_rho = zeros(N_reference)

for reference_idx in 1:N_reference
    target_SD =
        reference_terminal_SD_analytical_deg[reference_idx]

    matched_sigma_xi[reference_idx] =
        match_noise_amplitude(
            sigma_xi_grid,
            terminal_SD_xi_grid,
            target_SD
        )

    matched_sigma_omega[reference_idx] =
        match_noise_amplitude(
            sigma_omega_grid,
            terminal_SD_omega_grid,
            target_SD
        )

    matched_sigma_epsilon[reference_idx] =
        match_noise_amplitude(
            sigma_epsilon_grid,
            terminal_SD_epsilon_grid,
            target_SD
        )

    matched_sigma_rho[reference_idx] =
        match_noise_amplitude(
            sigma_rho_grid,
            terminal_SD_rho_grid,
            target_SD
        )
end

#MONTE CARLO METRICS FOR INTERNAL AND MATCHED EXTERNAL NOISE -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

N_sources = length(source_names)

matched_noise_values =
    zeros(N_reference, N_sources)

matched_noise_values[:,1] .= sigma_eta_list
matched_noise_values[:,2] .= matched_sigma_xi
matched_noise_values[:,3] .= matched_sigma_omega
matched_noise_values[:,4] .= matched_sigma_epsilon
matched_noise_values[:,5] .= matched_sigma_rho

analytical_terminal_SD_deg =
    zeros(N_reference, N_sources)

analytical_terminal_SD_deg[:,1] .=
    reference_terminal_SD_analytical_deg

MC_terminal_SD_deg =
    zeros(N_reference, N_sources)

MC_terminal_SD_percent =
    zeros(N_reference, N_sources)

mean_peak_velocity_deg_s =
    zeros(N_reference, N_sources)

sem_peak_velocity_deg_s =
    zeros(N_reference, N_sources)

mean_peak_time_s =
    zeros(N_reference, N_sources)

sem_peak_time_s =
    zeros(N_reference, N_sources)

mean_elbow_velocity_profile_deg_s =
    zeros(T, N_reference, N_sources)

elbow_angle_variance_profile_deg2 =
    zeros(T, N_reference, N_sources)

analytical_cost =
    zeros(N_reference, N_sources)

mean_MC_cost =
    zeros(N_reference, N_sources)

sem_MC_cost =
    zeros(N_reference, N_sources)

K_Frobenius =
    zeros(N_reference, N_sources)

L_Frobenius =
    zeros(N_reference, N_sources)

EV1_z =
    zeros(N_reference, N_sources)

PC1_z =
    zeros(
        dimension_of_state,
        N_reference,
        N_sources
    )

PC1_readout_mean =
    zeros(N_reference, N_sources)

PC1_readout_sem =
    zeros(N_reference, N_sources)

mean_peak_joint_velocity_norm_rad_s =
    zeros(N_reference, N_sources)

sem_peak_joint_velocity_norm_rad_s =
    zeros(N_reference, N_sources)

std_x =
    zeros(
        dimension_of_state,
        N_reference,
        N_sources
    )

std_z = similar(std_x)

std_u =
    zeros(
        dimension_of_control,
        N_reference,
        N_sources
    )

endpoint_bias_deg =
    zeros(N_reference, N_sources)

success_probability =
    zeros(N_reference, N_sources)

terminal_absolute_elbow_velocity_deg_s =
    zeros(N_reference, N_sources)

control_effort =
    zeros(N_reference, N_sources)

K_LCO_matched =
    zeros(
        dimension_of_state,
        dimension_of_observation,
        T-1,
        N_reference,
        N_sources-1
    )

L_LCO_matched =
    zeros(
        dimension_of_control,
        dimension_of_state,
        T-1,
        N_reference,
        N_sources-1
    )

completed_MC = falses(N_reference, N_sources)

pca_time_indices = 2:T-2

#The internal-noise reference already contains the requested 2000 trials.
for reference_idx in 1:N_reference

    MC_terminal_SD_deg[reference_idx,1],
    MC_terminal_SD_percent[reference_idx,1],
    mean_peak_velocity_deg_s[reference_idx,1],
    sem_peak_velocity_deg_s[reference_idx,1],
    mean_peak_time_s[reference_idx,1],
    sem_peak_time_s[reference_idx,1],
    mean_elbow_velocity_profile_deg_s[:,reference_idx,1],
    elbow_angle_variance_profile_deg2[:,reference_idx,1] =
        compute_MC_metrics(
            x_store_LCO_internal[reference_idx],
            time_vector
        )

    extended_metrics =
        compute_extended_metrics(
            x_store_LCO_internal[reference_idx],
            z_store_LCO_internal[reference_idx],
            K_LCO_internal[reference_idx],
            L_LCO_internal[reference_idx],
            pca_time_indices
        )

    analytical_cost[reference_idx,1] =
        analytical_cost_LCO_internal[reference_idx]

    mean_MC_cost[reference_idx,1] =
        mean(trial_costs_LCO_internal[reference_idx])

    sem_MC_cost[reference_idx,1] =
        std(trial_costs_LCO_internal[reference_idx])/
        sqrt(N_trials_MC)

    K_Frobenius[reference_idx,1] =
        extended_metrics.K_Frobenius

    L_Frobenius[reference_idx,1] =
        extended_metrics.L_Frobenius

    EV1_z[reference_idx,1] =
        extended_metrics.EV1

    PC1_z[:,reference_idx,1] .=
        extended_metrics.PC1

    PC1_readout_mean[reference_idx,1] =
        extended_metrics.PC1_readout_mean

    PC1_readout_sem[reference_idx,1] =
        extended_metrics.PC1_readout_sem

    mean_peak_joint_velocity_norm_rad_s[
        reference_idx,1
    ] =
        extended_metrics.mean_peak_velocity_norm

    sem_peak_joint_velocity_norm_rad_s[
        reference_idx,1
    ] =
        extended_metrics.sem_peak_velocity_norm

    std_x[:,reference_idx,1] .=
        extended_metrics.std_x

    std_z[:,reference_idx,1] .=
        extended_metrics.std_z

    std_u[:,reference_idx,1] .=
        extended_metrics.std_u

    endpoint_bias_deg[reference_idx,1] =
        extended_metrics.endpoint_bias_deg

    success_probability[reference_idx,1] =
        extended_metrics.success_probability

    terminal_absolute_elbow_velocity_deg_s[
        reference_idx,1
    ] =
        extended_metrics.terminal_absolute_velocity_deg_s

    control_effort[reference_idx,1] =
        extended_metrics.control_effort

    completed_MC[reference_idx,1] = true

end

if isfile(matched_results_file)

    println("Loading resumable matched Monte Carlo results from:")
    println(matched_results_file)

    saved_matched = JLD2.load(matched_results_file)

    if saved_matched["sigma_eta_list"] != sigma_eta_list ||
       !all(isapprox.(
           saved_matched["matched_noise_values"],
           matched_noise_values;
           atol=1e-12,
           rtol=1e-10
       ))
        error(
            "The saved matched noise values are inconsistent. "*
            "The existing results were not overwritten."
        )
    end

    analytical_terminal_SD_deg .=
        saved_matched["analytical_terminal_SD_deg"]

    MC_terminal_SD_deg .=
        saved_matched["MC_terminal_SD_deg"]

    MC_terminal_SD_percent .=
        saved_matched["MC_terminal_SD_percent"]

    mean_peak_velocity_deg_s .=
        saved_matched["mean_peak_velocity_deg_s"]

    sem_peak_velocity_deg_s .=
        saved_matched["sem_peak_velocity_deg_s"]

    mean_peak_time_s .=
        saved_matched["mean_peak_time_s"]

    sem_peak_time_s .=
        saved_matched["sem_peak_time_s"]

    mean_elbow_velocity_profile_deg_s .=
        saved_matched["mean_elbow_velocity_profile_deg_s"]

    elbow_angle_variance_profile_deg2 .=
        saved_matched["elbow_angle_variance_profile_deg2"]

    analytical_cost .=
        saved_matched["analytical_cost"]

    mean_MC_cost .=
        saved_matched["mean_MC_cost"]

    sem_MC_cost .=
        saved_matched["sem_MC_cost"]

    K_Frobenius .=
        saved_matched["K_Frobenius"]

    L_Frobenius .=
        saved_matched["L_Frobenius"]

    EV1_z .=
        saved_matched["EV1_z"]

    PC1_z .=
        saved_matched["PC1_z"]

    PC1_readout_mean .=
        saved_matched["PC1_readout_mean"]

    PC1_readout_sem .=
        saved_matched["PC1_readout_sem"]

    mean_peak_joint_velocity_norm_rad_s .=
        saved_matched[
            "mean_peak_joint_velocity_norm_rad_s"
        ]

    sem_peak_joint_velocity_norm_rad_s .=
        saved_matched[
            "sem_peak_joint_velocity_norm_rad_s"
        ]

    std_x .=
        saved_matched["std_x"]

    std_z .=
        saved_matched["std_z"]

    std_u .=
        saved_matched["std_u"]

    endpoint_bias_deg .=
        saved_matched["endpoint_bias_deg"]

    success_probability .=
        saved_matched["success_probability"]

    terminal_absolute_elbow_velocity_deg_s .=
        saved_matched[
            "terminal_absolute_elbow_velocity_deg_s"
        ]

    control_effort .=
        saved_matched["control_effort"]

    K_LCO_matched .= saved_matched["K_LCO_matched"]
    L_LCO_matched .= saved_matched["L_LCO_matched"]
    completed_MC .= saved_matched["completed_MC"]

end

for source_idx in 2:N_sources
    for reference_idx in 1:N_reference

        if completed_MC[reference_idx,source_idx]
            continue
        end

        matched_noise =
            matched_noise_values[reference_idx,source_idx]

        C_condition = C
        D_condition = D
        Ω_ξ_condition = Ω_ξ_baseline
        Ω_ω_condition = Ω_ω_baseline

        if source_idx == 2
            Ω_ξ_condition =
                Diagonal(fill(
                    matched_noise^2,
                    dimension_of_state
                ))

            source_label = "sigma_xi"
        elseif source_idx == 3
            Ω_ω_condition =
                Diagonal(fill(
                    matched_noise^2,
                    dimension_of_observation
                ))

            source_label = "sigma_omega"
        elseif source_idx == 4
            C_condition =
                matched_noise .* B

            source_label = "sigma_epsilon"
        else
            D_condition =
                matched_noise .* H

            source_label = "sigma_rho"
        end

        println(
            "\nMatched simulation ",
            source_label,
            ", reference ",
            reference_idx, "/", N_reference,
            ": ", matched_noise
        )

        if isapprox(
            matched_noise,
            0.01;
            atol=1e-12,
            rtol=1e-10
        )
            baseline_idx =
                argmin(abs.(sigma_eta_list))

            K_LCO = K_LCO_internal[baseline_idx]
            L_LCO = L_LCO_internal[baseline_idx]
        else
            L_TOD, K_TOD =
                Todorov_optimization_multidimensional_case(
                    A, B, H,
                    C_condition, D_condition, T,
                    dimension_of_state,
                    dimension_of_control,
                    dimension_of_observation,
                    x_1_mean,
                    Σ_1_x, Σ_1_z,
                    Q_matrix, R_matrix,
                    Ω_ξ_condition,
                    Ω_ω_condition,
                    Ω_η_zero,
                    N_iter_Todorov
                )

            K_LCO, L_LCO, _, _, _, _, _ =
                Optimal_Control_Estimation_FullProcedure(
                    copy(K_TOD), copy(L_TOD),
                    A, B, H,
                    C_condition, D_condition, T,
                    Ω_ξ_condition,
                    Ω_ω_condition,
                    Ω_η_zero,
                    Q_matrix, R_matrix,
                    Σ_1_x, Σ_1_z,
                    x_1_mean, z_1_mean;
                    ρ_list = ρ_list,
                    ρ_list_large = ρ_list_large,
                    N_iter_tot = N_iter_tot,
                    iters_K = iters_K,
                    iters_L = iters_L,
                    N_iter_GD_refine = 0,
                    N_iter_GD_from_IC = 0,
                    gd_method = gd_method,
                    η_GD = η_GD,
                    maxiter_GD = maxiter_GD,
                    tol_GD = tol_GD,
                    β1_Adam = β1_Adam,
                    β2_Adam = β2_Adam,
                    ϵ_Adam = ϵ_Adam,
                    use_GD_refinement = false,
                    use_GD_from_IC = false,
                    verbose = false
                )
        end

        K_LCO_matched[
            :,:,:,reference_idx,source_idx-1
        ] .= K_LCO

        L_LCO_matched[
            :,:,:,reference_idx,source_idx-1
        ] .= L_LCO

        analytical_terminal_SD_deg[
            reference_idx,source_idx
        ] =
            terminal_elbow_SD(
                K_LCO, L_LCO,
                A, B, H,
                C_condition, D_condition, T,
                x_1_mean, z_1_mean,
                Σ_1_x, Σ_1_z,
                Ω_ξ_condition,
                Ω_ω_condition,
                Ω_η_zero
            )

        target_terminal_SD =
            reference_terminal_SD_analytical_deg[reference_idx]

        relative_matching_error_percent =
            100*abs(
                analytical_terminal_SD_deg[
                    reference_idx,source_idx
                ] -
                target_terminal_SD
            )/target_terminal_SD

        println(
            "Analytical terminal-SD matching error: ",
            relative_matching_error_percent,
            "%"
        )

        if relative_matching_error_percent > 1.0
            @warn(
                "The exact reoptimized matched condition differs from "*
                "the interpolated target by more than 1%. Check the "*
                "saved numerical endpoint-SD values."
            )
        end

        M_LCO =
            zeros(
                dimension_of_state,
                dimension_of_state,
                T-1
            )

        for t in 1:T-1
            M_LCO[:,:,t] .=
                A .+
                B*L_LCO[:,:,t] .-
                K_LCO[:,:,t]*H
        end

        x_LCO, z_LCO, cost_over_time_LCO,
        mean_MC_cost_condition,
        sem_MC_cost_condition =
            sample_identifiability(
                random_seed,
                N_trials_MC,
                dimension_of_state,
                dimension_of_state,
                dimension_of_observation,
                T,
                x_1_mean, z_1_mean,
                Σ_1_x, Σ_1_z,
                A, B, C_condition, H, D_condition,
                K_LCO, M_LCO, L_LCO,
                Ω_ξ_condition,
                Ω_ω_condition,
                Ω_η_zero,
                Q_matrix, R_matrix
            )

        analytical_cost[reference_idx,source_idx] =
            expected_cost_raw_moments_propagation(
                T,
                dimension_of_state,
                dimension_of_control,
                dimension_of_observation,
                K_LCO, L_LCO,
                A, B, H,
                C_condition, D_condition,
                x_1_mean, z_1_mean,
                Σ_1_x, Σ_1_z,
                Q_matrix, R_matrix,
                Ω_ξ_condition,
                Ω_ω_condition,
                Ω_η_zero
            )

        mean_MC_cost[reference_idx,source_idx] =
            mean_MC_cost_condition

        sem_MC_cost[reference_idx,source_idx] =
            sem_MC_cost_condition

        MC_terminal_SD_deg[reference_idx,source_idx],
        MC_terminal_SD_percent[reference_idx,source_idx],
        mean_peak_velocity_deg_s[reference_idx,source_idx],
        sem_peak_velocity_deg_s[reference_idx,source_idx],
        mean_peak_time_s[reference_idx,source_idx],
        sem_peak_time_s[reference_idx,source_idx],
        mean_elbow_velocity_profile_deg_s[
            :,reference_idx,source_idx
        ],
        elbow_angle_variance_profile_deg2[
            :,reference_idx,source_idx
        ] =
            compute_MC_metrics(
                x_LCO,
                time_vector
            )

        extended_metrics =
            compute_extended_metrics(
                x_LCO, z_LCO,
                K_LCO, L_LCO,
                pca_time_indices
            )

        K_Frobenius[reference_idx,source_idx] =
            extended_metrics.K_Frobenius

        L_Frobenius[reference_idx,source_idx] =
            extended_metrics.L_Frobenius

        EV1_z[reference_idx,source_idx] =
            extended_metrics.EV1

        PC1_z[:,reference_idx,source_idx] .=
            extended_metrics.PC1

        PC1_readout_mean[reference_idx,source_idx] =
            extended_metrics.PC1_readout_mean

        PC1_readout_sem[reference_idx,source_idx] =
            extended_metrics.PC1_readout_sem

        mean_peak_joint_velocity_norm_rad_s[
            reference_idx,source_idx
        ] =
            extended_metrics.mean_peak_velocity_norm

        sem_peak_joint_velocity_norm_rad_s[
            reference_idx,source_idx
        ] =
            extended_metrics.sem_peak_velocity_norm

        std_x[:,reference_idx,source_idx] .=
            extended_metrics.std_x

        std_z[:,reference_idx,source_idx] .=
            extended_metrics.std_z

        std_u[:,reference_idx,source_idx] .=
            extended_metrics.std_u

        endpoint_bias_deg[reference_idx,source_idx] =
            extended_metrics.endpoint_bias_deg

        success_probability[reference_idx,source_idx] =
            extended_metrics.success_probability

        terminal_absolute_elbow_velocity_deg_s[
            reference_idx,source_idx
        ] =
            extended_metrics.terminal_absolute_velocity_deg_s

        control_effort[reference_idx,source_idx] =
            extended_metrics.control_effort

        completed_MC[reference_idx,source_idx] = true

        save_results(validation_only,
            matched_results_file,
            "Δt", Δt,
            "T_move", T_move,
            "T", T,
            "time_vector", time_vector,
            "N_trials_MC", N_trials_MC,
            "random_seed", random_seed,
            "sigma_eta_list", sigma_eta_list,
            "source_names", source_names,
            "matched_noise_values", matched_noise_values,
            "reference_terminal_SD_analytical_deg",
                reference_terminal_SD_analytical_deg,
            "analytical_terminal_SD_deg",
                analytical_terminal_SD_deg,
            "MC_terminal_SD_deg", MC_terminal_SD_deg,
            "MC_terminal_SD_percent", MC_terminal_SD_percent,
            "mean_peak_velocity_deg_s",
                mean_peak_velocity_deg_s,
            "sem_peak_velocity_deg_s",
                sem_peak_velocity_deg_s,
            "mean_peak_time_s", mean_peak_time_s,
            "sem_peak_time_s", sem_peak_time_s,
            "mean_elbow_velocity_profile_deg_s",
                mean_elbow_velocity_profile_deg_s,
            "elbow_angle_variance_profile_deg2",
                elbow_angle_variance_profile_deg2,
            "analytical_cost", analytical_cost,
            "mean_MC_cost", mean_MC_cost,
            "sem_MC_cost", sem_MC_cost,
            "K_Frobenius", K_Frobenius,
            "L_Frobenius", L_Frobenius,
            "EV1_z", EV1_z,
            "PC1_z", PC1_z,
            "PC1_readout_mean", PC1_readout_mean,
            "PC1_readout_sem", PC1_readout_sem,
            "mean_peak_joint_velocity_norm_rad_s",
                mean_peak_joint_velocity_norm_rad_s,
            "sem_peak_joint_velocity_norm_rad_s",
                sem_peak_joint_velocity_norm_rad_s,
            "std_x", std_x,
            "std_z", std_z,
            "std_u", std_u,
            "endpoint_bias_deg", endpoint_bias_deg,
            "success_probability", success_probability,
            "terminal_absolute_elbow_velocity_deg_s",
                terminal_absolute_elbow_velocity_deg_s,
            "control_effort", control_effort,
            "K_LCO_matched", K_LCO_matched,
            "L_LCO_matched", L_LCO_matched,
            "completed_MC", completed_MC,
            "sigma_xi_grid", sigma_xi_grid,
            "sigma_omega_grid", sigma_omega_grid,
            "terminal_SD_xi_grid", terminal_SD_xi_grid,
            "terminal_SD_omega_grid", terminal_SD_omega_grid
        )

    end
end


return JLD2.load(matched_results_file)
end
end
