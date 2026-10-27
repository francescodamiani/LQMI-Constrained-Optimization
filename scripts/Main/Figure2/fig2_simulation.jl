# Optional rerun of the final calibrated sweep; no plotting in this helper.
using LinearAlgebra, Statistics, Random

function simulate_fig2(data_dir; problem_data=reaching_problem(),
    sigma_eta_list=round.(range(0.0,stop=0.04,length=20),digits=3),
    N_trials_MC=2000, random_seed=1, N_iter_Todorov=500, N_iter_tot=500)

p = problem_data
N_noise = length(sigma_eta_list)
results = Dict{String,Any}(string(key)=>value for (key,value) in pairs(p))
results["time_vector"] = collect(0:p.Δt:p.T_move)
results["sigma_eta_list"] = collect(sigma_eta_list)
results["N_trials_MC"] = N_trials_MC
results["random_seed"] = random_seed
results["N_iter_Todorov"] = N_iter_Todorov
results["N_iter_tot_LCO"] = N_iter_tot
results["pseudodeterminant_tolerance"] = 1e-6
# Freeze the recovered final amplitude rather than refitting to new random draws.
results["calibrated_sigma_xi"] = p.σ_ξ
results["process_noise_origin"] = "fixed published calibration"
results["validation_only"] = !(p.T == 41 && p.T_move == 0.40 && p.Δt == 0.01 &&
    p.σ_ξ == 0.01 && N_trials_MC == 2000 && random_seed == 1 &&
    N_iter_Todorov == 500 && N_iter_tot == 500 &&
    sigma_eta_list == round.(range(0.0,stop=0.04,length=20),digits=3))
results["solver_settings"] = (;ρ_list=[0.0],ρ_list_large=[1e-10,1e-8,1e-6],
    iters_K=1,iters_L=1,N_iter_GD_refine=0,N_iter_GD_from_IC=0,
    gd_method=:adam,η_GD=1e-3,maxiter_GD=50,tol_GD=1e-6,
    β1_Adam=0.9,β2_Adam=0.999,ϵ_Adam=1e-8,
    use_GD_refinement=false,use_GD_from_IC=false)

for key in vcat(FIG2_TRAJECTORY_FIELDS,["K_TOD_all","L_TOD_all","K_LCO_all","L_LCO_all"])
    results[key] = Array{Float64,3}[]
end
for method in ["TOD","LCO"]
    for name in ["trial_costs","cost_history"]
        results["$(name)_$(method)"] = Vector{Float64}[]
    end
    for name in ["analytical_cost","mean_MC_cost","sem_MC_cost","mean_peak_elbow_velocity"]
        key = name == "mean_peak_elbow_velocity" ? "$(name)_$(method)_deg_s" : "$(name)_$(method)"
        results[key] = zeros(N_noise)
    end
    results["$(method)_converged"] = falses(N_noise)
end
results["cost_history_iterations_TOD"] = Vector{Int}[]
results["initial_cost"] = zeros(N_noise)
results["method_LCO"] = String[]
results["ρ_LCO"] = zeros(N_noise)
results["terminal_elbow_SD_pct_LCO"] = zeros(N_noise)
results["log_pseudodeterminant_K_LCO"] = zeros(N_noise)
results["log_pseudodeterminant_L_LCO"] = zeros(N_noise)

for (idx,σ_η_internal_noise) in enumerate(sigma_eta_list)
    println("Internal noise ",idx,"/",N_noise,": sigma_eta = ",σ_η_internal_noise)
    K_zero = zeros(p.dimension_of_state,p.dimension_of_observation,p.T-1)
    L_zero = zeros(p.dimension_of_control,p.dimension_of_state,p.T-1)
    results["initial_cost"][idx] = reaching_cost(p,L_zero,K_zero,σ_η_internal_noise)
    solution = solve_reaching(p,σ_η_internal_noise;N_iter_Todorov,N_iter_tot)
    for method in ["TOD","LCO"]
        L = getproperty(solution,Symbol("L_$(method)"))
        K = getproperty(solution,Symbol("K_$(method)"))
        cost = reaching_cost(p,L,K,σ_η_internal_noise)
        push!(results["K_$(method)_all"],K)
        push!(results["L_$(method)_all"],L)
        results["analytical_cost_$(method)"][idx] = cost
        # Reset the same explicit RNG for each method/level, as in the source.
        mc = simulate_reaching(p,L,K,σ_η_internal_noise,N_trials_MC;random_seed)
        # Preserve the source data schema; the core's native axes are dimension,time,trial.
        x = permutedims(mc.x_traces,(1,3,2))
        z = permutedims(mc.z_traces,(1,3,2))
        u = permutedims(mc.u_traces,(1,3,2))
        push!(results["x_store_$(method)"],x)
        push!(results["z_store_$(method)"],z)
        push!(results["u_store_$(method)"],u)
        push!(results["trial_costs_$(method)"],mc.cost_per_trial)
        results["mean_MC_cost_$(method)"][idx] = mc.empirical_cost
        results["sem_MC_cost_$(method)"][idx] = mc.empirical_cost_SE
        results["mean_peak_elbow_velocity_$(method)_deg_s"][idx] = compute_peak_elbow_velocity(x).mean_peak
        results["$(method)_converged"][idx] = isfinite(cost) &&
            (method == "TOD" || cost <= solution.cost_init_LCO)
    end
    push!(results["cost_history_TOD"],[results["initial_cost"][idx],results["analytical_cost_TOD"][idx]])
    push!(results["cost_history_iterations_TOD"],[0,N_iter_Todorov])
    push!(results["cost_history_LCO"],solution.cost_history_LCO)
    push!(results["method_LCO"],string(solution.method_LCO))
    results["ρ_LCO"][idx] = solution.ρ_LCO
    results["terminal_elbow_SD_pct_LCO"][idx] = 100*std(results["x_store_LCO"][idx][2,:,end])/(pi/9)
    results["log_pseudodeterminant_K_LCO"][idx] = mean_log_pdet_over_time(solution.K_LCO)
    results["log_pseudodeterminant_L_LCO"][idx] = mean_log_pdet_over_time(solution.L_LCO)
end
results["relative_cost_difference_percent"] = 100 .* (results["analytical_cost_TOD"] .-
    results["analytical_cost_LCO"]) ./ results["analytical_cost_TOD"]
save_fig2(data_dir,results)
return results
end
