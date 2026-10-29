# Sensory-model mismatch sweep: optimize under assumed noise and evaluate under true noise.
# Loaded only for explicit simulation reruns.
module Overestimation
using LinearAlgebra, Statistics, Random, Printf, MultivariateStats, JLD2
using ..IdentifiabilitySimulation: Todorov_optimization_multidimensional_case,
    Optimal_Control_Estimation_FullProcedure, expected_cost_raw_moments_propagation,
    sample_identifiability, save_results, mean_log_pdet_over_time

function terminal_elbow_SD(
    K_matrix, L_matrix, A, B, H, C, D, T,
    x_1_mean, z_1_mean, Σ_1_x, Σ_1_z,
    Ω_ξ, Ω_ω, Ω_η
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

    terminal_variance =
        S_xx[2,2] - x_mean[2]^2

    if terminal_variance < -1e-10
        error("The analytical terminal elbow variance is negative.")
    end

    return rad2deg(sqrt(max(terminal_variance, 0.0)))

end

function compute_MC_metrics(x_trajectories, time_vector)

    N_trials = size(x_trajectories,2)

    elbow_angle_deg =
        rad2deg.(x_trajectories[2,:,:])

    elbow_velocity_deg_s =
        rad2deg.(x_trajectories[5,:,:])

    peak_velocity = zeros(N_trials)
    peak_time = zeros(N_trials)

    for trial in 1:N_trials
        peak_velocity[trial],
        peak_idx =
            findmax(abs.(elbow_velocity_deg_s[trial,:]))

        peak_time[trial] =
            time_vector[peak_idx]
    end

    terminal_error_deg =
        elbow_angle_deg[:,end]

    terminal_SD_deg =
        std(terminal_error_deg)

    return (
        terminal_SD_deg=terminal_SD_deg,
        terminal_SD_percent=100*terminal_SD_deg/20.0,
        mean_peak_velocity_deg_s=mean(peak_velocity),
        sem_peak_velocity_deg_s=
            std(peak_velocity)/sqrt(N_trials),
        mean_peak_time_s=mean(peak_time),
        sem_peak_time_s=
            std(peak_time)/sqrt(N_trials),
        signed_terminal_mean_error_deg=
            mean(terminal_error_deg),
        absolute_terminal_mean_error_deg=
            abs(mean(terminal_error_deg)),
        mean_absolute_terminal_error_deg=
            mean(abs.(terminal_error_deg)),
        sem_terminal_mean_error_deg=
            std(terminal_error_deg)/sqrt(N_trials)
    )

end

mean_gain_Frobenius(gain_matrix) =
    mean([
        norm(gain_matrix[:,:,t])
        for t in 2:size(gain_matrix,3)-2
    ])


function run_experiment(reference_results, output_dir; validation_only=false)
#LOAD THE EXISTING BEHAVIORALLY CALIBRATED REACHING TASK -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

println("\n------- ELBOW REACHING: WRONG SENSORY-MODEL SWEEP ----------------------\n")



Δt = reference_results["Δt"]
T_move = reference_results["T_move"]
T = reference_results["T"]
time_vector = reference_results["time_vector"]

n = reference_results["dimension_of_state"]
m = reference_results["dimension_of_control"]
p = reference_results["dimension_of_observation"]

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

Ω_ξ = reference_results["Ω_ξ"]

sigma_omega_true =
    reference_results["σ_ω_add_sensory_noise"]

Ω_ω_true =
    Diagonal(fill(sigma_omega_true^2, p))

Ω_η_zero = zeros(n,n)

N_trials_MC = validation_only ? 64 : 2000
random_seed = reference_results["random_seed"]

if reference_results["N_trials_MC"] != N_trials_MC
    error("The saved internal-noise sweep must contain 2000 MC trials.")
end

#INTERNAL-NOISE REFERENCE: REUSE SAVED GAINS AND TRAJECTORIES -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

sigma_eta_list =
    reference_results["sigma_eta_list"]

K_LCO_internal =
    reference_results["K_LCO_all"]

L_LCO_internal =
    reference_results["L_LCO_all"]

x_store_LCO_internal =
    reference_results["x_store_LCO"]

analytical_cost_internal =
    reference_results["analytical_cost_LCO"]

mean_MC_cost_internal =
    reference_results["mean_MC_cost_LCO"]

sem_MC_cost_internal =
    reference_results["sem_MC_cost_LCO"]

N_internal = length(sigma_eta_list)

terminal_SD_analytical_internal_deg =
    zeros(N_internal)

terminal_SD_analytical_internal_percent =
    zeros(N_internal)

terminal_SD_MC_internal_deg =
    zeros(N_internal)

terminal_SD_MC_internal_percent =
    zeros(N_internal)

mean_peak_velocity_internal_deg_s =
    zeros(N_internal)

sem_peak_velocity_internal_deg_s =
    zeros(N_internal)

mean_peak_time_internal_s =
    zeros(N_internal)

sem_peak_time_internal_s =
    zeros(N_internal)

signed_terminal_mean_error_internal_deg =
    zeros(N_internal)

absolute_terminal_mean_error_internal_deg =
    zeros(N_internal)

mean_absolute_terminal_error_internal_deg =
    zeros(N_internal)

sem_terminal_mean_error_internal_deg =
    zeros(N_internal)

K_Frobenius_internal =
    zeros(N_internal)

L_Frobenius_internal =
    zeros(N_internal)

K_log_pdet_internal = zeros(N_internal)
L_log_pdet_internal = zeros(N_internal)

for noise_idx in 1:N_internal

    Ω_η =
        Diagonal(fill(
            sigma_eta_list[noise_idx]^2,
            n
        ))

    terminal_SD_analytical_internal_deg[noise_idx] =
        terminal_elbow_SD(
            K_LCO_internal[noise_idx],
            L_LCO_internal[noise_idx],
            A, B, H, C, D, T,
            x_1_mean, z_1_mean,
            Σ_1_x, Σ_1_z,
            Ω_ξ, Ω_ω_true, Ω_η
        )

    terminal_SD_analytical_internal_percent[noise_idx] =
        100*
        terminal_SD_analytical_internal_deg[noise_idx]/
        20.0

    internal_MC_metrics =
        compute_MC_metrics(
            x_store_LCO_internal[noise_idx],
            time_vector
        )

    terminal_SD_MC_internal_deg[noise_idx] =
        internal_MC_metrics.terminal_SD_deg

    terminal_SD_MC_internal_percent[noise_idx] =
        internal_MC_metrics.terminal_SD_percent

    mean_peak_velocity_internal_deg_s[noise_idx] =
        internal_MC_metrics.mean_peak_velocity_deg_s

    sem_peak_velocity_internal_deg_s[noise_idx] =
        internal_MC_metrics.sem_peak_velocity_deg_s

    mean_peak_time_internal_s[noise_idx] =
        internal_MC_metrics.mean_peak_time_s

    sem_peak_time_internal_s[noise_idx] =
        internal_MC_metrics.sem_peak_time_s

    signed_terminal_mean_error_internal_deg[noise_idx] =
        internal_MC_metrics.signed_terminal_mean_error_deg

    absolute_terminal_mean_error_internal_deg[noise_idx] =
        internal_MC_metrics.absolute_terminal_mean_error_deg

    mean_absolute_terminal_error_internal_deg[noise_idx] =
        internal_MC_metrics.mean_absolute_terminal_error_deg

    sem_terminal_mean_error_internal_deg[noise_idx] =
        internal_MC_metrics.sem_terminal_mean_error_deg

    K_Frobenius_internal[noise_idx] =
        mean_gain_Frobenius(
            K_LCO_internal[noise_idx]
        )

    L_Frobenius_internal[noise_idx] =
        mean_gain_Frobenius(
            L_LCO_internal[noise_idx]
        )

    K_log_pdet_internal[noise_idx] = mean_log_pdet_over_time(K_LCO_internal[noise_idx])
    L_log_pdet_internal[noise_idx] = mean_log_pdet_over_time(L_LCO_internal[noise_idx])

end

#WRONG ASSUMED SENSORY-NOISE GRID -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

N_assumed_noise = validation_only ? 3 : 21

sigma_omega_hat_list = sigma_omega_true .* 10.0 .^
    collect(range(-1.0,stop=1.0,length=N_assumed_noise))

baseline_idx =
    argmin(abs.(
        sigma_omega_hat_list .-
        sigma_omega_true
    ))

sigma_omega_hat_list[baseline_idx] =
    sigma_omega_true

if !isapprox(
    sigma_omega_hat_list[baseline_idx],
    sigma_omega_true;
    atol=0.0,
    rtol=0.0
)
    error("The assumed sensory-noise grid must contain the true amplitude exactly.")
end

#OPTIMIZATION PARAMETERS -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

#OUTPUT AND STORAGE -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

results_file = joinpath(output_dir,"sensory_overestimation.jld2")

K_LCO_wrong =
    zeros(n, p, T-1, N_assumed_noise)

L_LCO_wrong =
    zeros(m, n, T-1, N_assumed_noise)

assumed_model_cost =
    zeros(N_assumed_noise)

true_analytical_cost =
    zeros(N_assumed_noise)

mean_true_MC_cost =
    zeros(N_assumed_noise)

sem_true_MC_cost =
    zeros(N_assumed_noise)

terminal_SD_analytical_true_deg =
    zeros(N_assumed_noise)

terminal_SD_analytical_true_percent =
    zeros(N_assumed_noise)

terminal_SD_MC_true_deg =
    zeros(N_assumed_noise)

terminal_SD_MC_true_percent =
    zeros(N_assumed_noise)

mean_peak_velocity_wrong_deg_s =
    zeros(N_assumed_noise)

sem_peak_velocity_wrong_deg_s =
    zeros(N_assumed_noise)

mean_peak_time_wrong_s =
    zeros(N_assumed_noise)

sem_peak_time_wrong_s =
    zeros(N_assumed_noise)

signed_terminal_mean_error_wrong_deg =
    zeros(N_assumed_noise)

absolute_terminal_mean_error_wrong_deg =
    zeros(N_assumed_noise)

mean_absolute_terminal_error_wrong_deg =
    zeros(N_assumed_noise)

sem_terminal_mean_error_wrong_deg =
    zeros(N_assumed_noise)

K_Frobenius_wrong =
    zeros(N_assumed_noise)

L_Frobenius_wrong =
    zeros(N_assumed_noise)

K_log_pdet_wrong = zeros(N_assumed_noise)
L_log_pdet_wrong = zeros(N_assumed_noise)

completed_wrong_conditions =
    falses(N_assumed_noise)

if isfile(results_file)

    saved_results = JLD2.load(results_file)

    if get(saved_results,"model",nothing) != reference_results["model"] ||
       saved_results["sigma_omega_hat_list"] != sigma_omega_hat_list ||
       saved_results["N_trials_MC"] != N_trials_MC ||
       saved_results["sigma_eta_list"] != sigma_eta_list
        error(
            "The saved wrong-model grid differs from the requested grid."
        )
    end

    K_LCO_wrong .= saved_results["K_LCO_wrong"]
    L_LCO_wrong .= saved_results["L_LCO_wrong"]
    assumed_model_cost .= saved_results["assumed_model_cost"]
    true_analytical_cost .= saved_results["true_analytical_cost"]
    mean_true_MC_cost .= saved_results["mean_true_MC_cost"]
    sem_true_MC_cost .= saved_results["sem_true_MC_cost"]

    terminal_SD_analytical_true_deg .=
        saved_results["terminal_SD_analytical_true_deg"]

    terminal_SD_analytical_true_percent .=
        saved_results["terminal_SD_analytical_true_percent"]

    terminal_SD_MC_true_deg .=
        saved_results["terminal_SD_MC_true_deg"]

    terminal_SD_MC_true_percent .=
        saved_results["terminal_SD_MC_true_percent"]

    mean_peak_velocity_wrong_deg_s .=
        saved_results["mean_peak_velocity_wrong_deg_s"]

    sem_peak_velocity_wrong_deg_s .=
        saved_results["sem_peak_velocity_wrong_deg_s"]

    mean_peak_time_wrong_s .=
        saved_results["mean_peak_time_wrong_s"]

    sem_peak_time_wrong_s .=
        saved_results["sem_peak_time_wrong_s"]

    signed_terminal_mean_error_wrong_deg .=
        saved_results["signed_terminal_mean_error_wrong_deg"]

    absolute_terminal_mean_error_wrong_deg .=
        saved_results["absolute_terminal_mean_error_wrong_deg"]

    mean_absolute_terminal_error_wrong_deg .=
        saved_results["mean_absolute_terminal_error_wrong_deg"]

    sem_terminal_mean_error_wrong_deg .=
        saved_results["sem_terminal_mean_error_wrong_deg"]

    K_Frobenius_wrong .=
        saved_results["K_Frobenius_wrong"]

    L_Frobenius_wrong .=
        saved_results["L_Frobenius_wrong"]

    K_log_pdet_wrong .= saved_results["K_log_pdet_wrong"]
    L_log_pdet_wrong .= saved_results["L_log_pdet_wrong"]

    completed_wrong_conditions .=
        saved_results["completed_wrong_conditions"]

    println(
        "Resuming from ",
        count(completed_wrong_conditions),
        "/",
        N_assumed_noise,
        " completed conditions."
    )

end

function save_current_results()

    save_results(validation_only,
        results_file,
        "Δt", Δt,
        "T_move", T_move,
        "T", T,
        "time_vector", time_vector,
        "N_trials_MC", N_trials_MC,
        "random_seed", random_seed,
        "model", reference_results["model"],
        "reference_dataset", "data/published/Main/Figure2/",
        "pseudodeterminant_tolerance", 1e-6,
        "sigma_omega_true", sigma_omega_true,
        "sigma_omega_hat_list", sigma_omega_hat_list,
        "sigma_eta_list", sigma_eta_list,
        "K_LCO_wrong", K_LCO_wrong,
        "L_LCO_wrong", L_LCO_wrong,
        "assumed_model_cost", assumed_model_cost,
        "true_analytical_cost", true_analytical_cost,
        "mean_true_MC_cost", mean_true_MC_cost,
        "sem_true_MC_cost", sem_true_MC_cost,
        "terminal_SD_analytical_true_deg",
            terminal_SD_analytical_true_deg,
        "terminal_SD_analytical_true_percent",
            terminal_SD_analytical_true_percent,
        "terminal_SD_MC_true_deg",
            terminal_SD_MC_true_deg,
        "terminal_SD_MC_true_percent",
            terminal_SD_MC_true_percent,
        "mean_peak_velocity_wrong_deg_s",
            mean_peak_velocity_wrong_deg_s,
        "sem_peak_velocity_wrong_deg_s",
            sem_peak_velocity_wrong_deg_s,
        "mean_peak_time_wrong_s",
            mean_peak_time_wrong_s,
        "sem_peak_time_wrong_s",
            sem_peak_time_wrong_s,
        "signed_terminal_mean_error_wrong_deg",
            signed_terminal_mean_error_wrong_deg,
        "absolute_terminal_mean_error_wrong_deg",
            absolute_terminal_mean_error_wrong_deg,
        "mean_absolute_terminal_error_wrong_deg",
            mean_absolute_terminal_error_wrong_deg,
        "sem_terminal_mean_error_wrong_deg",
            sem_terminal_mean_error_wrong_deg,
        "K_Frobenius_wrong", K_Frobenius_wrong,
        "L_Frobenius_wrong", L_Frobenius_wrong,
        "K_log_pdet_wrong", K_log_pdet_wrong,
        "L_log_pdet_wrong", L_log_pdet_wrong,
        "completed_wrong_conditions",
            completed_wrong_conditions,
        "terminal_SD_analytical_internal_deg",
            terminal_SD_analytical_internal_deg,
        "terminal_SD_analytical_internal_percent",
            terminal_SD_analytical_internal_percent,
        "terminal_SD_MC_internal_deg",
            terminal_SD_MC_internal_deg,
        "terminal_SD_MC_internal_percent",
            terminal_SD_MC_internal_percent,
        "analytical_cost_internal",
            analytical_cost_internal,
        "mean_MC_cost_internal",
            mean_MC_cost_internal,
        "sem_MC_cost_internal",
            sem_MC_cost_internal,
        "mean_peak_velocity_internal_deg_s",
            mean_peak_velocity_internal_deg_s,
        "sem_peak_velocity_internal_deg_s",
            sem_peak_velocity_internal_deg_s,
        "mean_peak_time_internal_s",
            mean_peak_time_internal_s,
        "sem_peak_time_internal_s",
            sem_peak_time_internal_s,
        "signed_terminal_mean_error_internal_deg",
            signed_terminal_mean_error_internal_deg,
        "absolute_terminal_mean_error_internal_deg",
            absolute_terminal_mean_error_internal_deg,
        "mean_absolute_terminal_error_internal_deg",
            mean_absolute_terminal_error_internal_deg,
        "sem_terminal_mean_error_internal_deg",
            sem_terminal_mean_error_internal_deg,
        "K_Frobenius_internal",
            K_Frobenius_internal,
        "L_Frobenius_internal",
            L_Frobenius_internal,
        "K_log_pdet_internal", K_log_pdet_internal,
        "L_log_pdet_internal", L_log_pdet_internal
    )

end

#WRONG-MODEL OPTIMIZATION AND TRUE-SYSTEM EVALUATION -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

baseline_internal_idx =
    argmin(abs.(sigma_eta_list))

for noise_idx in 1:N_assumed_noise

    if completed_wrong_conditions[noise_idx]
        println(
            "Already completed sigma_omega_hat = ",
            sigma_omega_hat_list[noise_idx]
        )
        continue
    end

    sigma_omega_hat =
        sigma_omega_hat_list[noise_idx]

    println(
        "\nAssumed sensory noise ",
        noise_idx,
        "/",
        N_assumed_noise,
        ": sigma_omega_hat = ",
        sigma_omega_hat
    )

    Ω_ω_assumed =
        Diagonal(fill(sigma_omega_hat^2, p))

    if noise_idx == baseline_idx

        K_LCO =
            copy(K_LCO_internal[baseline_internal_idx])

        L_LCO =
            copy(L_LCO_internal[baseline_internal_idx])

        x_true =
            x_store_LCO_internal[baseline_internal_idx]

        mean_true_MC_cost[noise_idx] =
            mean_MC_cost_internal[baseline_internal_idx]

        sem_true_MC_cost[noise_idx] =
            sem_MC_cost_internal[baseline_internal_idx]

        println(
            "Reused the identical saved sigma_eta = 0 baseline solution."
        )

    else

        L_TOD, K_TOD =
            Todorov_optimization_multidimensional_case(
                A, B, H, C, D, T,
                n, m, p,
                x_1_mean, Σ_1_x, Σ_1_z,
                Q_matrix, R_matrix,
                Ω_ξ, Ω_ω_assumed, Ω_η_zero,
                N_iter_Todorov
            )

        K_LCO, L_LCO, _, _, _, _, _ =
            Optimal_Control_Estimation_FullProcedure(
                copy(K_TOD), copy(L_TOD),
                A, B, H, C, D, T,
                Ω_ξ, Ω_ω_assumed, Ω_η_zero,
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

        M_LCO =
            zeros(n,n,T-1)

        for t in 1:T-1
            M_LCO[:,:,t] .=
                A .+
                B*L_LCO[:,:,t] .-
                K_LCO[:,:,t]*H
        end

        x_true, _, _, mean_true_MC_cost[noise_idx],
            sem_true_MC_cost[noise_idx] =
                sample_identifiability(
                    random_seed,
                    N_trials_MC,
                    n, n, p, T,
                    x_1_mean, z_1_mean,
                    Σ_1_x, Σ_1_z,
                    A, B, C, H, D,
                    K_LCO, M_LCO, L_LCO,
                    Ω_ξ, Ω_ω_true, Ω_η_zero,
                    Q_matrix, R_matrix
                )

    end

    if !all(isfinite, K_LCO) ||
       !all(isfinite, L_LCO)
        error(
            "Non-finite gains at sigma_omega_hat = $sigma_omega_hat"
        )
    end

    K_LCO_wrong[:,:,:,noise_idx] .= K_LCO
    L_LCO_wrong[:,:,:,noise_idx] .= L_LCO

    assumed_model_cost[noise_idx] =
        expected_cost_raw_moments_propagation(
            T, n, m, p,
            K_LCO, L_LCO,
            A, B, H, C, D,
            x_1_mean, z_1_mean,
            Σ_1_x, Σ_1_z,
            Q_matrix, R_matrix,
            Ω_ξ, Ω_ω_assumed, Ω_η_zero
        )

    true_analytical_cost[noise_idx] =
        expected_cost_raw_moments_propagation(
            T, n, m, p,
            K_LCO, L_LCO,
            A, B, H, C, D,
            x_1_mean, z_1_mean,
            Σ_1_x, Σ_1_z,
            Q_matrix, R_matrix,
            Ω_ξ, Ω_ω_true, Ω_η_zero
        )

    terminal_SD_analytical_true_deg[noise_idx] =
        terminal_elbow_SD(
            K_LCO, L_LCO,
            A, B, H, C, D, T,
            x_1_mean, z_1_mean,
            Σ_1_x, Σ_1_z,
            Ω_ξ, Ω_ω_true, Ω_η_zero
        )

    terminal_SD_analytical_true_percent[noise_idx] =
        100*
        terminal_SD_analytical_true_deg[noise_idx]/
        20.0

    wrong_MC_metrics =
        compute_MC_metrics(x_true, time_vector)

    terminal_SD_MC_true_deg[noise_idx] =
        wrong_MC_metrics.terminal_SD_deg

    terminal_SD_MC_true_percent[noise_idx] =
        wrong_MC_metrics.terminal_SD_percent

    mean_peak_velocity_wrong_deg_s[noise_idx] =
        wrong_MC_metrics.mean_peak_velocity_deg_s

    sem_peak_velocity_wrong_deg_s[noise_idx] =
        wrong_MC_metrics.sem_peak_velocity_deg_s

    mean_peak_time_wrong_s[noise_idx] =
        wrong_MC_metrics.mean_peak_time_s

    sem_peak_time_wrong_s[noise_idx] =
        wrong_MC_metrics.sem_peak_time_s

    signed_terminal_mean_error_wrong_deg[noise_idx] =
        wrong_MC_metrics.signed_terminal_mean_error_deg

    absolute_terminal_mean_error_wrong_deg[noise_idx] =
        wrong_MC_metrics.absolute_terminal_mean_error_deg

    mean_absolute_terminal_error_wrong_deg[noise_idx] =
        wrong_MC_metrics.mean_absolute_terminal_error_deg

    sem_terminal_mean_error_wrong_deg[noise_idx] =
        wrong_MC_metrics.sem_terminal_mean_error_deg

    K_Frobenius_wrong[noise_idx] =
        mean_gain_Frobenius(K_LCO)

    L_Frobenius_wrong[noise_idx] =
        mean_gain_Frobenius(L_LCO)

    K_log_pdet_wrong[noise_idx] = mean_log_pdet_over_time(K_LCO)
    L_log_pdet_wrong[noise_idx] = mean_log_pdet_over_time(L_LCO)

    completed_wrong_conditions[noise_idx] = true

    save_current_results()

    println(
        "True analytical endpoint SD (% of 20 deg): ",
        terminal_SD_analytical_true_percent[noise_idx]
    )

    println(
        "Assumed-model cost: ",
        assumed_model_cost[noise_idx],
        " | true-system cost: ",
        true_analytical_cost[noise_idx]
    )

end

if !all(completed_wrong_conditions)
    error("At least one wrong-model condition is incomplete.")
end


return JLD2.load(results_file)
end
end
