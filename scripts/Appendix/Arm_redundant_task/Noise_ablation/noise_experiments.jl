# Numerical reruns only; the entry point loads this file exclusively with --simulate.
using LinearAlgebra, Random, Statistics, JLD2
include(joinpath(@__DIR__,"..","arm_task.jl"))
include(joinpath(@__DIR__,"component_variances.jl"))

const ORIGINAL_NOISE_CONDITIONS = ["full_noise","no_sensory_noise",
    "no_motor_multiplicative_noise","no_process_noise","no_internal_noise"]
const HIGH_NOISE_ABLATIONS = ["no_process_noise","no_additive_sensory_noise",
    "no_motor_multiplicative_noise","no_sensory_multiplicative_noise","no_internal_noise"]

function noise_problem(;T=300,s=1.0)
    p = arm_problem(;T)
    # The process amplitude remains 0.1. Scale sensory/motor amplitudes,
    # squaring only additive sensory amplitude for its covariance.
    σ_ω_add_sensory_noise = 0.1s
    σ_ρ_mult_sensory_noise = 0.1s
    σ_ϵ_mult_control_noise = 0.1s
    C = σ_ϵ_mult_control_noise .* p.B
    D = σ_ρ_mult_sensory_noise .* p.H
    Ω_ω = Diagonal(fill(σ_ω_add_sensory_noise^2,p.dimension_of_observation))
    return merge(p,(;σ_ω_add_sensory_noise,σ_ρ_mult_sensory_noise,σ_ϵ_mult_control_noise,C,D,Ω_ω))
end

function ablation_noise_parameters(p,σ_η,condition)
    C,D,Ω_ξ,Ω_ω = p.C,p.D,p.Ω_ξ,p.Ω_ω
    Ω_η = Diagonal(fill(σ_η^2,p.dimension_of_state))
    if condition == "no_sensory_noise"
        D = zeros(size(p.D))
        Ω_ω = zeros(size(p.Ω_ω))
    elseif condition == "no_additive_sensory_noise"
        Ω_ω = zeros(size(p.Ω_ω))
    elseif condition == "no_sensory_multiplicative_noise"
        D = zeros(size(p.D))
    elseif condition == "no_motor_multiplicative_noise"
        C = zeros(size(p.C))
    elseif condition == "no_process_noise"
        Ω_ξ = zeros(size(p.Ω_ξ))
    elseif condition == "no_internal_noise"
        Ω_η = zeros(p.dimension_of_state,p.dimension_of_state)
    elseif condition != "full_noise"
        error("Unknown noise condition: $condition")
    end
    return (;C,D,Ω_ξ,Ω_ω,Ω_η)
end

function fixed_policy_component_variances(p,K_matrix,L_matrix,σ_η,condition;temporal_window=1:p.T-2)
    noise = ablation_noise_parameters(p,σ_η,condition)
    # K/L and the observation matrix H remain exactly the full-model policy.
    variance_x,variance_z = compute_analytical_component_variances(
        K_matrix,L_matrix,p.A,p.B,p.H,noise.C,noise.D,p.T,temporal_window,
        p.x_1_mean,p.z_1_mean,p.Σ_1_x,p.Σ_1_z,noise.Ω_ξ,noise.Ω_ω,noise.Ω_η)
    min(minimum(variance_x),minimum(variance_z)) >= -1e-10 || error("Negative propagated variance")
    return max.(variance_x,0.0),max.(variance_z,0.0)
end

function monte_carlo_component_metrics(p,K_matrix,L_matrix,σ_η,N_trials_MC;
    random_seed=1234,temporal_window=1:p.T-2,condition="full_noise")
    noise = ablation_noise_parameters(p,σ_η,condition)
    mc = monte_carlo_simulation(L_matrix,K_matrix,p.A,p.B,p.H,noise.C,noise.D,p.T,
        p.dimension_of_state,p.dimension_of_control,p.dimension_of_observation,
        p.x_1_mean,p.z_1_mean,p.Σ_1_x,p.Σ_1_z,p.Q_matrix,p.R_matrix,
        Matrix(noise.Ω_ξ),Matrix(noise.Ω_ω),Matrix(noise.Ω_η),N_trials_MC;
        rng=MersenneTwister(random_seed))
    # Core axes: dimension,time,trial. Statistics use sample variance
    # across trials at each time, then average those variances over 1:T-2.
    variance_x_MC = component_sample_variances(mc.x_traces,temporal_window)
    variance_z_MC = component_sample_variances(mc.z_traces,temporal_window)
    trial_peak_velocity = [maximum(norm(mc.x_traces[4:6,t,trial]) for t in 1:p.T) for trial in 1:N_trials_MC]
    return (;variance_x_MC,variance_z_MC,mean_MC_cost=mc.empirical_cost,
        sem_MC_cost=mc.empirical_cost_SE,mean_peak_velocity=mean(trial_peak_velocity),
        sem_peak_velocity=std(trial_peak_velocity)/sqrt(N_trials_MC))
end

function full_model_policies(p,sigma_eta_list;N_iter_Todorov=500,N_iter_tot=500)
    N_noise = length(sigma_eta_list)
    data = Dict{String,Any}("sigma_eta_list"=>sigma_eta_list,"completed_solutions"=>falses(N_noise),
        "K_TOD_all"=>zeros(6,3,p.T-1,N_noise),"L_TOD_all"=>zeros(9,6,p.T-1,N_noise),
        "K_LCO_all"=>zeros(6,3,p.T-1,N_noise),"L_LCO_all"=>zeros(9,6,p.T-1,N_noise),
        "T"=>p.T,"model"=>p,"N_iter_Todorov"=>N_iter_Todorov,"N_iter_LCO"=>N_iter_tot,
        "rho_list"=>[0.0],"rho_list_large"=>[1e-10,1e-8,1e-6],"solver_diagnostics"=>NamedTuple[])
    for (idx,σ_η) in enumerate(sigma_eta_list)
        # Each noise level starts from its own independently recomputed TOD.
        solution = solve_arm(p,σ_η;N_iter_Todorov,N_iter_tot,verbose=false)
        for (method,L,K) in [("TOD",solution.L_TOD,solution.K_TOD),
            ("LCO",solution.L_Full_Solver,solution.K_Full_Solver)]
            data["K_$(method)_all"][:,:,:,idx] .= K
            data["L_$(method)_all"][:,:,:,idx] .= L
        end
        push!(data["solver_diagnostics"],(;sigma_eta=σ_η,solution.best_method,solution.best_ρ,
            solution.cost_init,solution.cost_Full_Solver,solution.best_cost_list,solution.t_TOD,solution.t_LCO))
        data["completed_solutions"][idx] = true
    end
    return data
end

function simulate_behavioral_variability(p,policies;N_trials_MC=5000)
    sigma_eta_list = policies["sigma_eta_list"]
    data = Dict{String,Any}("sigma_eta_list"=>sigma_eta_list,
        "temporal_window"=>collect(1:p.T-2),"N_trials_MC"=>N_trials_MC,"random_seed"=>1234)
    for method in ["TOD","LCO"]
        for state in ["x","z"],mode in ["analytical","MC"]
            data["variance_$(state)_$(mode)_$(method)"] = zeros(6,length(sigma_eta_list))
        end
        data["mean_MC_cost_$(method)"] = zeros(length(sigma_eta_list))
        data["sem_MC_cost_$(method)"] = zeros(length(sigma_eta_list))
        for (idx,σ_η) in enumerate(sigma_eta_list)
            K_matrix,L_matrix = policies["K_$(method)_all"][:,:,:,idx],policies["L_$(method)_all"][:,:,:,idx]
            variance_x,variance_z = fixed_policy_component_variances(p,K_matrix,L_matrix,σ_η,"full_noise")
            mc = monte_carlo_component_metrics(p,K_matrix,L_matrix,σ_η,N_trials_MC)
            data["variance_x_analytical_$(method)"][:,idx] .= variance_x
            data["variance_z_analytical_$(method)"][:,idx] .= variance_z
            data["variance_x_MC_$(method)"][:,idx] .= mc.variance_x_MC
            data["variance_z_MC_$(method)"][:,idx] .= mc.variance_z_MC
            data["mean_MC_cost_$(method)"][idx] = mc.mean_MC_cost
            data["sem_MC_cost_$(method)"][idx] = mc.sem_MC_cost
        end
        baseline = findfirst(iszero,sigma_eta_list)
        baseline === nothing && error("Behavioral sweep requires zero internal noise")
        for state in ["x","z"],mode in ["analytical","MC"]
            variance = data["variance_$(state)_$(mode)_$(method)"]
            # The archive stores variance increases; SD increases
            # for the printed tables are derived after square roots.
            data["relative_$(state)_$(mode)_$(method)"] = 100 .* (variance ./ variance[:,baseline] .- 1)
        end
    end
    return data
end

function simulate_original_ablation(p,policies)
    sigma_eta_list = policies["sigma_eta_list"]
    conditions = ORIGINAL_NOISE_CONDITIONS
    variance_x = zeros(6,length(sigma_eta_list),2,length(conditions))
    variance_z = similar(variance_x)
    for (idx,σ_η) in enumerate(sigma_eta_list),(method_idx,method) in enumerate(["TOD","LCO"])
        K_matrix,L_matrix = policies["K_$(method)_all"][:,:,:,idx],policies["L_$(method)_all"][:,:,:,idx]
        for (condition_idx,condition) in enumerate(conditions)
            vx,vz = fixed_policy_component_variances(p,K_matrix,L_matrix,σ_η,condition)
            variance_x[:,idx,method_idx,condition_idx] .= vx
            variance_z[:,idx,method_idx,condition_idx] .= vz
        end
    end
    data = Dict{String,Any}("sigma_eta_list"=>sigma_eta_list,"temporal_window"=>collect(1:p.T-2),
        "method_names"=>["TOD","LCO"],"condition_names"=>copy(conditions),
        "ablation_names"=>conditions[2:end],"component_names"=>copy(ARM_COMPONENT_NAMES),
        "variance_x"=>variance_x,"variance_z"=>variance_z)
    for (state,variance) in [("x",variance_x),("z",variance_z)]
        full = variance[:,:,:,1:1]
        all(>(0),full) || error("Nonpositive full-model variance")
        # Keep the source's raw VARIANCE reductions alongside the derived SD table.
        data["relative_reduction_$(state)"] = 100 .* (full .- variance[:,:,:,2:end]) ./ full
    end
    return data
end

# These gain metrics average over interior time slices.
function mean_log_pseudodeterminant(gain_matrix;tolerance=1e-6)
    return mean(sum(log,filter(>(tolerance),svdvals(gain_matrix[:,:,t]));init=0.0)
        for t in 2:size(gain_matrix,3)-1)
end
mean_frobenius_norm(gain_matrix) = mean(norm(gain_matrix[:,:,t]) for t in 2:size(gain_matrix,3)-1)

function simulate_high_noise(p,sigma_eta_list;N_trials_MC=2000,N_iter_Todorov=500,N_iter_tot=500)
    policies = full_model_policies(p,sigma_eta_list;N_iter_Todorov,N_iter_tot)
    N_noise = length(sigma_eta_list)
    data = Dict{String,Any}("external_noise_multipliers"=>[10.0],"sigma_eta_list"=>sigma_eta_list,
        "method_names"=>["TOD","LCO"],"ablation_names"=>copy(HIGH_NOISE_ABLATIONS),
        "N_trials_MC"=>N_trials_MC,"random_seed"=>1234,"temporal_window"=>collect(1:p.T-2),
        "pseudodeterminant_tolerance"=>1e-6,"K_all"=>zeros(6,3,p.T-1,2,1,N_noise),
        "L_all"=>zeros(9,6,p.T-1,2,1,N_noise),"T"=>p.T,"Delta_t"=>p.Δt,
        "N_iter_Todorov"=>N_iter_Todorov,"N_iter_LCO"=>N_iter_tot,
        "sigma_xi"=>0.1,"sigma_omega_baseline"=>0.1,"sigma_rho_baseline"=>0.1,"sigma_epsilon_baseline"=>0.1,
        "rho_list"=>[0.0],"rho_list_large"=>[1e-10,1e-8,1e-6],"model"=>p)
    for key in ["analytical_cost","mean_MC_cost","sem_MC_cost","gain_K_log_pdet","gain_L_log_pdet",
        "gain_K_frobenius","gain_L_frobenius","mean_peak_velocity","sem_peak_velocity","runtime"]
        data[key] = zeros(2,1,N_noise)
    end
    for key in ["variance_x_MC","variance_z_MC","variance_x_analytical","variance_z_analytical"]
        data[key] = zeros(6,2,1,N_noise)
    end
    for key in ["ablation_variance_x","ablation_reduction_x"]
        data[key] = zeros(6,5,2,1,N_noise)
    end
    for key in ["LCO_initial_cost","LCO_selected_rho","LCO_final_residual","LCO_maximum_positive_increment"]
        data[key] = zeros(1,N_noise)
    end
    data["LCO_number_updates"],data["LCO_number_passes"] = zeros(Int,1,N_noise),zeros(Int,1,N_noise)
    data["LCO_solver_method"] = fill("",1,N_noise)
    data["LCO_cost_history"] = zeros(2N_iter_tot+1,1,N_noise)
    data["optimization_status"],data["analysis_status"] = fill("OK",2,1,N_noise),fill("OK",2,1,N_noise)
    data["failure_reason"] = fill("",2,1,N_noise)
    data["completed_conditions"] = falses(1,N_noise)
    for (idx,σ_η) in enumerate(sigma_eta_list)
        diagnostics = policies["solver_diagnostics"][idx]
        history = diagnostics.best_cost_list
        data["LCO_initial_cost"][1,idx] = diagnostics.cost_init
        data["LCO_selected_rho"][1,idx] = diagnostics.best_ρ
        data["LCO_solver_method"][1,idx] = diagnostics.best_method
        data["LCO_number_updates"][1,idx] = length(history)-1
        data["LCO_number_passes"][1,idx] = (length(history)-1)÷2
        data["LCO_cost_history"][:,1,idx] .= history
        data["LCO_final_residual"][1,idx] = abs(history[end]-history[end-2])/max(1.0,abs(history[end-2]))
        data["LCO_maximum_positive_increment"][1,idx] = max(maximum(diff(history)),0.0)
        data["runtime"][:,1,idx] .= [diagnostics.t_TOD,diagnostics.t_LCO]
        for (method_idx,method) in enumerate(["TOD","LCO"])
            K_matrix,L_matrix = policies["K_$(method)_all"][:,:,:,idx],policies["L_$(method)_all"][:,:,:,idx]
            data["K_all"][:,:,:,method_idx,1,idx] .= K_matrix
            data["L_all"][:,:,:,method_idx,1,idx] .= L_matrix
            data["analytical_cost"][method_idx,1,idx] = arm_analytical_cost(p,L_matrix,K_matrix,σ_η)
            data["gain_K_log_pdet"][method_idx,1,idx] = mean_log_pseudodeterminant(K_matrix)
            data["gain_L_log_pdet"][method_idx,1,idx] = mean_log_pseudodeterminant(L_matrix)
            data["gain_K_frobenius"][method_idx,1,idx] = mean_frobenius_norm(K_matrix)
            data["gain_L_frobenius"][method_idx,1,idx] = mean_frobenius_norm(L_matrix)
            mc = monte_carlo_component_metrics(p,K_matrix,L_matrix,σ_η,N_trials_MC)
            for key in ["mean_MC_cost","sem_MC_cost","mean_peak_velocity","sem_peak_velocity"]
                data[key][method_idx,1,idx] = getproperty(mc,Symbol(key))
            end
            data["variance_x_MC"][:,method_idx,1,idx] .= mc.variance_x_MC
            data["variance_z_MC"][:,method_idx,1,idx] .= mc.variance_z_MC
            full_x,full_z = fixed_policy_component_variances(p,K_matrix,L_matrix,σ_η,"full_noise")
            data["variance_x_analytical"][:,method_idx,1,idx] .= full_x
            data["variance_z_analytical"][:,method_idx,1,idx] .= full_z
            for (a,condition) in enumerate(HIGH_NOISE_ABLATIONS)
                ablated_x,_ = fixed_policy_component_variances(p,K_matrix,L_matrix,σ_η,condition)
                data["ablation_variance_x"][:,a,method_idx,1,idx] .= ablated_x
                data["ablation_reduction_x"][:,a,method_idx,1,idx] .= 100 .* (full_x .- ablated_x) ./ full_x
            end
        end
        data["completed_conditions"][1,idx] = true
    end
    costs = data["analytical_cost"]
    data["relative_LCO_improvement"] = 100 .* (costs[1,:,:] .- costs[2,:,:]) ./ costs[1,:,:]
    return data
end

function save_noise_result(output_dir,filename,data;validation_only=false)
    data["generated_by_current_core"] = true
    data["validation_only"] = validation_only
    JLD2.jldopen(joinpath(output_dir,filename),"w") do file
        for (key,value) in data
            file[key] = value
        end
    end
    return data
end

function simulate_noise_experiments(output_dir;experiments=["behavioral","original_ablation","high_noise"],
    validation_only=false,reoptimize=false)
    p = noise_problem(T=validation_only ? 12 : 300)
    N_iter_Todorov,N_iter_tot = validation_only ? (3,3) : (500,500)
    results = Dict{String,Any}()
    if any(in(["behavioral","original_ablation"]),experiments)
        sigma_eta_list = validation_only ? [0.0,0.1] : round.(range(0.0,stop=0.4,length=20),digits=2)
        if validation_only || reoptimize
            policies = full_model_policies(p,sigma_eta_list;N_iter_Todorov,N_iter_tot)
        else
            published = joinpath(PROJECT_ROOT,"data","published","Appendix","Arm_redundant_task","Noise_ablation")
            policies = load_benchmark(published,NOISE_DATA_FILES["policies"])
            policies["sigma_eta_list"] == sigma_eta_list || error("Saved policy grid differs")
            all(policies["completed_solutions"]) || error("Saved policies are incomplete")
        end
        policies["policy_provenance"] = validation_only || reoptimize ? "recomputed full model" : "archived full model"
        save_noise_result(output_dir,NOISE_DATA_FILES["policies"],policies;validation_only)
        if "behavioral" in experiments
            data = simulate_behavioral_variability(p,policies;N_trials_MC=validation_only ? 64 : 5000)
            results["behavioral"] = save_noise_result(output_dir,NOISE_DATA_FILES["behavioral"],data;validation_only)
        end
        if "original_ablation" in experiments
            data = simulate_original_ablation(p,policies)
            results["original_ablation"] = save_noise_result(output_dir,NOISE_DATA_FILES["original_ablation"],data;validation_only)
        end
    end
    if "high_noise" in experiments
        p_high = noise_problem(T=p.T,s=10.0)
        data = simulate_high_noise(p_high,validation_only ? [0.0,0.1] : [0.0,0.1,0.2,0.3,0.4];
            N_trials_MC=validation_only ? 64 : 2000,N_iter_Todorov,N_iter_tot)
        results["high_noise"] = save_noise_result(output_dir,NOISE_DATA_FILES["high_noise"],data;validation_only)
    end
    return results
end
