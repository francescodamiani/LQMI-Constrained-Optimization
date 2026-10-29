# Sensory-underestimation and internal-noise raw sweeps.
# Loaded only for explicit simulation reruns.
module Underestimation
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

        x_mean_new = A*x_mean .+ B*L_t*z_mean
        z_mean_new = K_t*H*x_mean .+ M_t*z_mean

        S_xx = S_xx_new
        S_zz = S_zz_new
        S_xz = S_xz_new
        S_zx = S_zx_new

        x_mean = x_mean_new
        z_mean = z_mean_new

    end

    terminal_variance = S_xx[2,2] - x_mean[2]^2

    if terminal_variance < -1e-10
        error("The analytical terminal elbow variance is negative.")
    end

    return rad2deg(sqrt(max(terminal_variance, 0.0)))

end

function compute_MC_metrics(x_trajectories, time_vector)

    N_trials = size(x_trajectories,2)

    elbow_angle_deg = rad2deg.(x_trajectories[2,:,:])
    elbow_velocity_deg_s = rad2deg.(x_trajectories[5,:,:])

    peak_velocity = zeros(N_trials)
    peak_time = zeros(N_trials)

    for trial in 1:N_trials
        peak_velocity[trial], peak_idx =
            findmax(abs.(elbow_velocity_deg_s[trial,:]))

        peak_time[trial] = time_vector[peak_idx]
    end

    terminal_error_deg = elbow_angle_deg[:,end]

    return (
        terminal_SD_deg=std(terminal_error_deg),
        terminal_SD_percent=100*std(terminal_error_deg)/20.0,
        mean_peak_velocity_deg_s=mean(peak_velocity),
        sem_peak_velocity_deg_s=
            std(peak_velocity)/sqrt(N_trials),
        mean_peak_time_s=mean(peak_time),
        sem_peak_time_s=std(peak_time)/sqrt(N_trials),
        signed_terminal_mean_error_deg=mean(terminal_error_deg),
        absolute_terminal_mean_error_deg=abs(mean(terminal_error_deg)),
        terminal_RMSE_deg=sqrt(mean(terminal_error_deg.^2)),
        mean_elbow_angle_deg=
            vec(mean(elbow_angle_deg, dims=1)),
        mean_elbow_velocity_deg_s=
            vec(mean(elbow_velocity_deg_s, dims=1))
    )

end

mean_gain_Frobenius(gain_matrix) =
    mean([
        norm(gain_matrix[:,:,t])
        for t in 2:size(gain_matrix,3)-2
    ])



function run_experiment(reference_results, output_dir; validation_only=false)
#LOAD THE EXISTING REACHING TASK -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

println("\n------- ELBOW REACHING: SENSORY-NOISE UNDERESTIMATION TEST ----------------------\n")



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

sigma_omega_true = reference_results["σ_ω_add_sensory_noise"]
Ω_ω_true = Diagonal(fill(sigma_omega_true^2, p))

# Preserve the ten-point relative range: true amplitude down to 0.32 times true.
sigma_omega_hat_list = sigma_omega_true .*
    collect(range(1.0, stop=0.32, length=validation_only ? 3 : 10))
sigma_eta_list = reference_results["sigma_eta_list"]

N_wrong = length(sigma_omega_hat_list)
N_internal = length(sigma_eta_list)

baseline_wrong_idx =
    findfirst(==(sigma_omega_true), sigma_omega_hat_list)

baseline_internal_idx =
    findfirst(==(0.0), sigma_eta_list)

#OPTIMIZATION AND MONTE CARLO SETTINGS -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

#Use 1000 Monte Carlo trials, matching the published dataset.
N_trials_MC = validation_only ? 64 : 1000
random_seed = reference_results["random_seed"]
pseudodeterminant_tolerance = 1e-6

function optimize_LCO_solution(Ω_ω_assumed, Ω_η_assumed)

    L_TOD, K_TOD =
        Todorov_optimization_multidimensional_case(
            A, B, H, C, D, T,
            n, m, p,
            x_1_mean, Σ_1_x, Σ_1_z,
            Q_matrix, R_matrix,
            Ω_ξ, Ω_ω_assumed, Ω_η_assumed,
            N_iter_Todorov
        )

    K_LCO, L_LCO, _, _, _, _, _ =
        Optimal_Control_Estimation_FullProcedure(
            copy(K_TOD), copy(L_TOD),
            A, B, H, C, D, T,
            Ω_ξ, Ω_ω_assumed, Ω_η_assumed,
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

    if !all(isfinite, K_LCO) || !all(isfinite, L_LCO)
        error("LCO returned non-finite gains.")
    end

    return K_LCO, L_LCO

end


function simulate_true_system(K_matrix, L_matrix, Ω_η_true)

    M_matrix = zeros(n,n,T-1)

    for t in 1:T-1
        M_matrix[:,:,t] .=
            A .+ B*L_matrix[:,:,t] .- K_matrix[:,:,t]*H
    end

    x_trajectories, _, _, mean_MC_cost, sem_MC_cost =
        sample_identifiability(
            random_seed, N_trials_MC,
            n, n, p, T,
            x_1_mean, z_1_mean,
            Σ_1_x, Σ_1_z,
            A, B, C, H, D,
            K_matrix, M_matrix, L_matrix,
            Ω_ξ, Ω_ω_true, Ω_η_true,
            Q_matrix, R_matrix
        )

    return x_trajectories, mean_MC_cost, sem_MC_cost

end

#OUTPUT ARRAYS -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

results_file = joinpath(output_dir,"sensory_underestimation.jld2")

K_internal = zeros(n,p,T-1,N_internal)
L_internal = zeros(m,n,T-1,N_internal)
K_wrong = zeros(n,p,T-1,N_wrong)
L_wrong = zeros(m,n,T-1,N_wrong)

internal_cost = fill(NaN,N_internal)
internal_MC_cost = fill(NaN,N_internal)
internal_MC_cost_SEM = fill(NaN,N_internal)
internal_SD_analytical_deg = fill(NaN,N_internal)
internal_SD_MC_deg = fill(NaN,N_internal)
internal_peak_velocity = fill(NaN,N_internal)
internal_peak_velocity_SEM = fill(NaN,N_internal)
internal_peak_time = fill(NaN,N_internal)
internal_peak_time_SEM = fill(NaN,N_internal)
internal_bias = fill(NaN,N_internal)
internal_absolute_bias = fill(NaN,N_internal)
internal_RMSE = fill(NaN,N_internal)
internal_K_norm = fill(NaN,N_internal)
internal_L_norm = fill(NaN,N_internal)
internal_K_log_pdet = fill(NaN,N_internal)
internal_L_log_pdet = fill(NaN,N_internal)
internal_mean_angle = fill(NaN,N_internal,T)
internal_mean_velocity = fill(NaN,N_internal,T)
completed_internal = falses(N_internal)

wrong_assumed_cost = fill(NaN,N_wrong)
wrong_true_cost = fill(NaN,N_wrong)
wrong_MC_cost = fill(NaN,N_wrong)
wrong_MC_cost_SEM = fill(NaN,N_wrong)
wrong_SD_analytical_deg = fill(NaN,N_wrong)
wrong_SD_MC_deg = fill(NaN,N_wrong)
wrong_peak_velocity = fill(NaN,N_wrong)
wrong_peak_velocity_SEM = fill(NaN,N_wrong)
wrong_peak_time = fill(NaN,N_wrong)
wrong_peak_time_SEM = fill(NaN,N_wrong)
wrong_bias = fill(NaN,N_wrong)
wrong_absolute_bias = fill(NaN,N_wrong)
wrong_RMSE = fill(NaN,N_wrong)
wrong_K_norm = fill(NaN,N_wrong)
wrong_L_norm = fill(NaN,N_wrong)
wrong_K_log_pdet = fill(NaN,N_wrong)
wrong_L_log_pdet = fill(NaN,N_wrong)
wrong_mean_angle = fill(NaN,N_wrong,T)
wrong_mean_velocity = fill(NaN,N_wrong,T)
completed_wrong = falses(N_wrong)

if isfile(results_file)

    saved_results = JLD2.load(results_file)

    if get(saved_results,"model",nothing) != reference_results["model"] ||
       saved_results["sigma_omega_hat_list"] != sigma_omega_hat_list ||
       saved_results["sigma_eta_list"] != sigma_eta_list ||
       saved_results["N_trials_MC"] != N_trials_MC
        error("The saved configuration differs from the requested dense-sweep configuration.")
    end

    for (variable_name, target) in [
        "K_internal"=>K_internal, "L_internal"=>L_internal, "K_wrong"=>K_wrong, "L_wrong"=>L_wrong,
        "internal_cost"=>internal_cost, "internal_MC_cost"=>internal_MC_cost, "internal_MC_cost_SEM"=>internal_MC_cost_SEM,
        "internal_SD_analytical_deg"=>internal_SD_analytical_deg, "internal_SD_MC_deg"=>internal_SD_MC_deg,
        "internal_peak_velocity"=>internal_peak_velocity, "internal_peak_velocity_SEM"=>internal_peak_velocity_SEM,
        "internal_peak_time"=>internal_peak_time, "internal_peak_time_SEM"=>internal_peak_time_SEM,
        "internal_bias"=>internal_bias, "internal_absolute_bias"=>internal_absolute_bias, "internal_RMSE"=>internal_RMSE,
        "internal_K_norm"=>internal_K_norm, "internal_L_norm"=>internal_L_norm,
        "internal_mean_angle"=>internal_mean_angle, "internal_mean_velocity"=>internal_mean_velocity,
        "completed_internal"=>completed_internal, "wrong_assumed_cost"=>wrong_assumed_cost, "wrong_true_cost"=>wrong_true_cost,
        "wrong_MC_cost"=>wrong_MC_cost, "wrong_MC_cost_SEM"=>wrong_MC_cost_SEM,
        "wrong_SD_analytical_deg"=>wrong_SD_analytical_deg, "wrong_SD_MC_deg"=>wrong_SD_MC_deg,
        "wrong_peak_velocity"=>wrong_peak_velocity, "wrong_peak_velocity_SEM"=>wrong_peak_velocity_SEM,
        "wrong_peak_time"=>wrong_peak_time, "wrong_peak_time_SEM"=>wrong_peak_time_SEM,
        "wrong_bias"=>wrong_bias, "wrong_absolute_bias"=>wrong_absolute_bias, "wrong_RMSE"=>wrong_RMSE,
        "wrong_K_norm"=>wrong_K_norm, "wrong_L_norm"=>wrong_L_norm,
        "wrong_mean_angle"=>wrong_mean_angle, "wrong_mean_velocity"=>wrong_mean_velocity, "completed_wrong"=>completed_wrong
    ]
        target .= saved_results[variable_name]
    end

    saved_pseudodeterminant_tolerance = NaN

    if haskey(saved_results,"pseudodeterminant_tolerance")
        saved_pseudodeterminant_tolerance =
            saved_results["pseudodeterminant_tolerance"]
    end

    if haskey(saved_results,"internal_K_log_pdet") &&
       isapprox(
           saved_pseudodeterminant_tolerance,
           pseudodeterminant_tolerance;
           atol=0.0,
           rtol=0.0
       )
        internal_K_log_pdet .= saved_results["internal_K_log_pdet"]
        internal_L_log_pdet .= saved_results["internal_L_log_pdet"]
        wrong_K_log_pdet .= saved_results["wrong_K_log_pdet"]
        wrong_L_log_pdet .= saved_results["wrong_L_log_pdet"]
    else
        for noise_idx in findall(completed_internal)
            internal_K_log_pdet[noise_idx] =
                mean_log_pdet_over_time(
                    K_internal[:,:,:,noise_idx]
                )

            internal_L_log_pdet[noise_idx] =
                mean_log_pdet_over_time(
                    L_internal[:,:,:,noise_idx]
                )
        end

        for noise_idx in findall(completed_wrong)
            wrong_K_log_pdet[noise_idx] =
                mean_log_pdet_over_time(
                    K_wrong[:,:,:,noise_idx]
                )

            wrong_L_log_pdet[noise_idx] =
                mean_log_pdet_over_time(
                    L_wrong[:,:,:,noise_idx]
                )
        end
    end

    println(
        "Resuming: internal ", count(completed_internal), "/", N_internal,
        ", wrong model ", count(completed_wrong), "/", N_wrong
    )

end


function save_current_results()

    save_results(validation_only,
        results_file,
        "Delta_t", Δt,
        "T_move", T_move,
        "T", T,
        "time_vector", time_vector,
        "N_trials_MC", N_trials_MC,
        "random_seed", random_seed,
        "model", reference_results["model"],
        "reference_dataset", "data/published/Main/Figure2/",
        "pseudodeterminant_tolerance", pseudodeterminant_tolerance,
        "sigma_omega_true", sigma_omega_true,
        "sigma_omega_hat_list", sigma_omega_hat_list,
        "sigma_eta_list", sigma_eta_list,
        "K_internal", K_internal,
        "L_internal", L_internal,
        "K_wrong", K_wrong,
        "L_wrong", L_wrong,
        "internal_cost", internal_cost,
        "internal_MC_cost", internal_MC_cost,
        "internal_MC_cost_SEM", internal_MC_cost_SEM,
        "internal_SD_analytical_deg", internal_SD_analytical_deg,
        "internal_SD_MC_deg", internal_SD_MC_deg,
        "internal_peak_velocity", internal_peak_velocity,
        "internal_peak_velocity_SEM", internal_peak_velocity_SEM,
        "internal_peak_time", internal_peak_time,
        "internal_peak_time_SEM", internal_peak_time_SEM,
        "internal_bias", internal_bias,
        "internal_absolute_bias", internal_absolute_bias,
        "internal_RMSE", internal_RMSE,
        "internal_K_norm", internal_K_norm,
        "internal_L_norm", internal_L_norm,
        "internal_K_log_pdet", internal_K_log_pdet,
        "internal_L_log_pdet", internal_L_log_pdet,
        "internal_mean_angle", internal_mean_angle,
        "internal_mean_velocity", internal_mean_velocity,
        "completed_internal", completed_internal,
        "wrong_assumed_cost", wrong_assumed_cost,
        "wrong_true_cost", wrong_true_cost,
        "wrong_MC_cost", wrong_MC_cost,
        "wrong_MC_cost_SEM", wrong_MC_cost_SEM,
        "wrong_SD_analytical_deg", wrong_SD_analytical_deg,
        "wrong_SD_MC_deg", wrong_SD_MC_deg,
        "wrong_peak_velocity", wrong_peak_velocity,
        "wrong_peak_velocity_SEM", wrong_peak_velocity_SEM,
        "wrong_peak_time", wrong_peak_time,
        "wrong_peak_time_SEM", wrong_peak_time_SEM,
        "wrong_bias", wrong_bias,
        "wrong_absolute_bias", wrong_absolute_bias,
        "wrong_RMSE", wrong_RMSE,
        "wrong_K_norm", wrong_K_norm,
        "wrong_L_norm", wrong_L_norm,
        "wrong_K_log_pdet", wrong_K_log_pdet,
        "wrong_L_log_pdet", wrong_L_log_pdet,
        "wrong_mean_angle", wrong_mean_angle,
        "wrong_mean_velocity", wrong_mean_velocity,
        "completed_wrong", completed_wrong
    )

end

#CORRECT SENSORY MODEL WITH INTERNAL NOISE -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

println("\n--- Correct sensory model: internal-noise sweep ---")

for noise_idx in 1:N_internal

    if completed_internal[noise_idx]
        println("Already completed sigma_eta = ", sigma_eta_list[noise_idx])
        continue
    end

    sigma_eta = sigma_eta_list[noise_idx]
    Ω_η = Diagonal(fill(sigma_eta^2,n))

    println(
        "\nInternal noise ", noise_idx, "/", N_internal,
        ": sigma_eta = ", sigma_eta
    )

    # Reuse Figure 2 policies and the first N_trials_MC saved trials.
    # This preserves this analysis's 1000-trial convention without reoptimization.
    K_LCO = reference_results["K_LCO_all"][noise_idx]
    L_LCO = reference_results["L_LCO_all"][noise_idx]
    x_trajectories = view(reference_results["x_store_LCO"][noise_idx],:,1:N_trials_MC,:)
    reference_costs = reference_results["trial_costs_LCO"][noise_idx][1:N_trials_MC]
    internal_MC_cost[noise_idx] = mean(reference_costs)
    internal_MC_cost_SEM[noise_idx] = std(reference_costs)/sqrt(N_trials_MC)

    metrics = compute_MC_metrics(x_trajectories,time_vector)

    K_internal[:,:,:,noise_idx] .= K_LCO
    L_internal[:,:,:,noise_idx] .= L_LCO

    internal_cost[noise_idx] =
        expected_cost_raw_moments_propagation(
            T, n, m, p,
            K_LCO, L_LCO,
            A, B, H, C, D,
            x_1_mean, z_1_mean,
            Σ_1_x, Σ_1_z,
            Q_matrix, R_matrix,
            Ω_ξ, Ω_ω_true, Ω_η
        )

    internal_SD_analytical_deg[noise_idx] =
        terminal_elbow_SD(
            K_LCO, L_LCO,
            A, B, H, C, D, T,
            x_1_mean, z_1_mean,
            Σ_1_x, Σ_1_z,
            Ω_ξ, Ω_ω_true, Ω_η
        )

    internal_SD_MC_deg[noise_idx] = metrics.terminal_SD_deg
    internal_peak_velocity[noise_idx] = metrics.mean_peak_velocity_deg_s
    internal_peak_velocity_SEM[noise_idx] = metrics.sem_peak_velocity_deg_s
    internal_peak_time[noise_idx] = metrics.mean_peak_time_s
    internal_peak_time_SEM[noise_idx] = metrics.sem_peak_time_s
    internal_bias[noise_idx] = metrics.signed_terminal_mean_error_deg
    internal_absolute_bias[noise_idx] = metrics.absolute_terminal_mean_error_deg
    internal_RMSE[noise_idx] = metrics.terminal_RMSE_deg
    internal_K_norm[noise_idx] = mean_gain_Frobenius(K_LCO)
    internal_L_norm[noise_idx] = mean_gain_Frobenius(L_LCO)
    internal_K_log_pdet[noise_idx] =
        mean_log_pdet_over_time(
            K_LCO
        )

    internal_L_log_pdet[noise_idx] =
        mean_log_pdet_over_time(
            L_LCO
        )
    internal_mean_angle[noise_idx,:] .= metrics.mean_elbow_angle_deg
    internal_mean_velocity[noise_idx,:] .= metrics.mean_elbow_velocity_deg_s

    completed_internal[noise_idx] = true
    save_current_results()

    @printf(
        "sigma_eta=%.4f | endpoint SD=%.5f deg | peak=%.5f deg/s | cost=%.7g\n",
        sigma_eta,
        internal_SD_analytical_deg[noise_idx],
        internal_peak_velocity[noise_idx],
        internal_cost[noise_idx]
    )

end


#WRONG ASSUMED SENSORY MODEL, EVALUATED IN THE TRUE SYSTEM -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

println("\n--- Sensory-noise underestimation: true sigma_omega remains ",sigma_omega_true," ---")

for noise_idx in 1:N_wrong

    if completed_wrong[noise_idx]
        println(
            "Already completed sigma_omega_assumed = ",
            sigma_omega_hat_list[noise_idx]
        )
        continue
    end

    sigma_omega_hat = sigma_omega_hat_list[noise_idx]
    Ω_ω_assumed = Diagonal(fill(sigma_omega_hat^2,p))
    Ω_η_zero = zeros(n,n)

    println(
        "\nAssumed sensory noise ", noise_idx, "/", N_wrong,
        ": sigma_omega_assumed = ", sigma_omega_hat,
        " | sigma_omega_true = ", sigma_omega_true
    )

    if noise_idx == baseline_wrong_idx
        K_LCO = copy(K_internal[:,:,:,baseline_internal_idx])
        L_LCO = copy(L_internal[:,:,:,baseline_internal_idx])
    else
        K_LCO, L_LCO =
            optimize_LCO_solution(Ω_ω_assumed, Ω_η_zero)
    end

    if noise_idx == baseline_wrong_idx
        x_trajectories = view(reference_results["x_store_LCO"][baseline_internal_idx],:,1:N_trials_MC,:)
        wrong_MC_cost[noise_idx] = internal_MC_cost[baseline_internal_idx]
        wrong_MC_cost_SEM[noise_idx] = internal_MC_cost_SEM[baseline_internal_idx]
    else
        x_trajectories,
        wrong_MC_cost[noise_idx],
        wrong_MC_cost_SEM[noise_idx] =
            simulate_true_system(K_LCO, L_LCO, Ω_η_zero)
    end

    metrics = compute_MC_metrics(x_trajectories,time_vector)

    K_wrong[:,:,:,noise_idx] .= K_LCO
    L_wrong[:,:,:,noise_idx] .= L_LCO

    wrong_assumed_cost[noise_idx] =
        expected_cost_raw_moments_propagation(
            T, n, m, p,
            K_LCO, L_LCO,
            A, B, H, C, D,
            x_1_mean, z_1_mean,
            Σ_1_x, Σ_1_z,
            Q_matrix, R_matrix,
            Ω_ξ, Ω_ω_assumed, Ω_η_zero
        )

    wrong_true_cost[noise_idx] =
        expected_cost_raw_moments_propagation(
            T, n, m, p,
            K_LCO, L_LCO,
            A, B, H, C, D,
            x_1_mean, z_1_mean,
            Σ_1_x, Σ_1_z,
            Q_matrix, R_matrix,
            Ω_ξ, Ω_ω_true, Ω_η_zero
        )

    wrong_SD_analytical_deg[noise_idx] =
        terminal_elbow_SD(
            K_LCO, L_LCO,
            A, B, H, C, D, T,
            x_1_mean, z_1_mean,
            Σ_1_x, Σ_1_z,
            Ω_ξ, Ω_ω_true, Ω_η_zero
        )

    wrong_SD_MC_deg[noise_idx] = metrics.terminal_SD_deg
    wrong_peak_velocity[noise_idx] = metrics.mean_peak_velocity_deg_s
    wrong_peak_velocity_SEM[noise_idx] = metrics.sem_peak_velocity_deg_s
    wrong_peak_time[noise_idx] = metrics.mean_peak_time_s
    wrong_peak_time_SEM[noise_idx] = metrics.sem_peak_time_s
    wrong_bias[noise_idx] = metrics.signed_terminal_mean_error_deg
    wrong_absolute_bias[noise_idx] = metrics.absolute_terminal_mean_error_deg
    wrong_RMSE[noise_idx] = metrics.terminal_RMSE_deg
    wrong_K_norm[noise_idx] = mean_gain_Frobenius(K_LCO)
    wrong_L_norm[noise_idx] = mean_gain_Frobenius(L_LCO)
    wrong_K_log_pdet[noise_idx] =
        mean_log_pdet_over_time(
            K_LCO
        )

    wrong_L_log_pdet[noise_idx] =
        mean_log_pdet_over_time(
            L_LCO
        )
    wrong_mean_angle[noise_idx,:] .= metrics.mean_elbow_angle_deg
    wrong_mean_velocity[noise_idx,:] .= metrics.mean_elbow_velocity_deg_s

    completed_wrong[noise_idx] = true
    save_current_results()

    @printf(
        "assumed=%.4f | true endpoint SD=%.5f deg | peak=%.5f deg/s | assumed cost=%.7g | true cost=%.7g\n",
        sigma_omega_hat,
        wrong_SD_analytical_deg[noise_idx],
        wrong_peak_velocity[noise_idx],
        wrong_assumed_cost[noise_idx],
        wrong_true_cost[noise_idx]
    )

end


save_current_results()

return JLD2.load(results_file)
end
end
