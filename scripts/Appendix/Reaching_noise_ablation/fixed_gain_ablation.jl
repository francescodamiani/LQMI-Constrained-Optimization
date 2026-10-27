# Migrated from analyze_elbow_reaching_calibrated_fixed_gain_noise_ablation.jl.
# Reuse the actual public calibrated model, original policies and existing moments.
using LinearAlgebra, Statistics, SHA
include(joinpath(@__DIR__,"..","..","Main","Figure2","reaching_task.jl"))
include(joinpath(@__DIR__,"..","Arm_redundant_task","Noise_ablation","component_variances.jl"))

function calibrated_reaching_reference()
    metadata=joinpath(PROJECT_ROOT,"data","published","Main","Figure2","sweep.jld2")
    isfile(metadata) || error("Missing published Figure 2 metadata/gains: $metadata")
    saved_results=JLD2.load(metadata) # No trajectory archive is loaded or copied.
    p=reaching_problem()
    # Refuse a scientific change rather than silently substituting other settings.
    for key in [:Δt,:T_move,:T,:dimension_of_state,:dimension_of_control,:dimension_of_observation,
        :x_1_mean,:z_1_mean,:Σ_1_x,:Σ_1_z,:A,:B,:H,:C,:D,:Ω_ξ,:Ω_ω,:Q_matrix,:R_matrix,
        :σ_ξ,:σ_ω_add_sensory_noise,:σ_ϵ_mult_control_noise,:σ_ρ_mult_sensory_noise]
        isequal(saved_results[string(key)],getproperty(p,key)) || error("Figure 2 model mismatch: $key")
    end
    return p,saved_results,metadata
end

function reaching_ablation_noise(p,σ_η,condition)
    condition in REACHING_ABLATION_CONDITIONS || error("Unknown reaching ablation: $condition")
    C_full,D_full,Ω_ξ_full,Ω_ω_full=p.C,p.D,p.Ω_ξ,p.Ω_ω
    Ω_η_full=Diagonal(fill(σ_η^2,p.dimension_of_state))
    C_test=condition=="no_motor_multiplicative_noise" ? zeros(size(C_full)) : C_full
    D_test=condition in ["no_sensory_multiplicative_noise","no_all_sensory_noise"] ? zeros(size(D_full)) : D_full
    Ω_ξ_test=condition=="no_process_noise" ? zeros(size(Ω_ξ_full)) : Ω_ξ_full
    Ω_ω_test=condition in ["no_additive_sensory_noise","no_all_sensory_noise"] ? zeros(size(Ω_ω_full)) : Ω_ω_full
    Ω_η_test=condition=="no_internal_noise" ? zeros(size(Ω_η_full)) : Ω_η_full
    return (;C=C_test,D=D_test,Ω_ξ=Ω_ξ_test,Ω_ω=Ω_ω_test,Ω_η=Ω_η_test)
end

function physical_variances_over_time(p,K_matrix,L_matrix,noise)
    # The existing exact analysis helper returns temporal means of variances.
    # A singleton window returns the original instantaneous variance. Reusing
    # it avoids another implementation of moment propagation or a core change.
    variance_x_over_time=zeros(p.dimension_of_state,p.T)
    for t in 1:p.T
        variance_x,_=compute_analytical_component_variances(K_matrix,L_matrix,p.A,p.B,p.H,
            noise.C,noise.D,p.T,t:t,p.x_1_mean,p.z_1_mean,p.Σ_1_x,p.Σ_1_z,
            noise.Ω_ξ,noise.Ω_ω,noise.Ω_η)
        variance_x_over_time[:,t].=variance_x
    end
    minimum(variance_x_over_time)>=-1e-10 || error("A propagated physical-state variance is negative beyond numerical tolerance")
    return max.(variance_x_over_time,0.0)
end

function recompute_reaching_ablation(output_dir)
    p,saved_results,metadata=calibrated_reaching_reference()
    sigma_eta_list=saved_results["sigma_eta_list"]
    target_sigma_eta=[0.0,0.01,0.02,0.03,0.04]
    selected_indices=[argmin(abs.(sigma_eta_list.-value)) for value in target_sigma_eta]
    sigma_eta_list[selected_indices]==REACHING_ABLATION_LEVELS || error("Selected Figure 2 levels differ from the historical analysis")
    ablation_names=REACHING_ABLATION_CONDITIONS
    results=DataFrame(sigma_eta=Float64[],condition=String[],terminal_elbow_SD_pct=Float64[],
        mean_angle_SD=Float64[],mean_velocity_SD=Float64[],terminal_elbow_SD_reduction_pct=Float64[],
        mean_angle_SD_reduction_pct=Float64[],mean_velocity_SD_reduction_pct=Float64[])
    variance_x_all=zeros(6,p.T,length(ablation_names),length(selected_indices))
    movement_amplitude=pi/9

    for (level,idx) in enumerate(selected_indices)
        sigma_eta=sigma_eta_list[idx]
        K_LCO=saved_results["K_LCO_all"][idx]
        L_LCO=saved_results["L_LCO_all"][idx]
        # Reuse Figure 2's validated core cost call to verify the loaded policy
        # against its original full-model objective before any ablation.
        isapprox(reaching_cost(p,L_LCO,K_LCO,sigma_eta),saved_results["analytical_cost_LCO"][idx];rtol=1e-10) ||
            error("Archived Figure 2 policy cost mismatch at sigma_eta=$sigma_eta")
        condition_metrics=zeros(length(ablation_names),3)
        for (condition_idx,condition) in enumerate(ablation_names)
            noise=reaching_ablation_noise(p,sigma_eta,condition)
            variance_x_over_time=physical_variances_over_time(p,K_LCO,L_LCO,noise)
            variance_x_all[:,:,condition_idx,level].=variance_x_over_time
            terminal_elbow_SD_pct=100*sqrt(variance_x_over_time[2,end])/movement_amplitude
            angle_SD=sqrt.(variance_x_over_time[1:3,2:p.T])
            velocity_SD=sqrt.(variance_x_over_time[4:6,2:p.T])
            condition_metrics[condition_idx,1]=terminal_elbow_SD_pct
            condition_metrics[condition_idx,2]=mean(angle_SD)
            condition_metrics[condition_idx,3]=mean(velocity_SD)
        end
        full_metrics=condition_metrics[1,:]
        for condition_idx in eachindex(ablation_names)
            reductions=100 .* (full_metrics.-condition_metrics[condition_idx,:])./full_metrics
            push!(results,(sigma_eta,ablation_names[condition_idx],condition_metrics[condition_idx,:]...,reductions...))
        end
    end
    data=Dict("df_results"=>results,"variance_x_over_time"=>variance_x_all,
        "sigma_eta_list"=>sigma_eta_list[selected_indices],"target_sigma_eta"=>target_sigma_eta,
        "selected_indices"=>selected_indices,"ablation_names"=>ablation_names,
        "temporal_window"=>collect(2:p.T),"movement_amplitude_rad"=>movement_amplitude,
        "variance_tolerance"=>1e-10,"generated_by_current_core"=>true,
        "fixed_policy"=>true,"optimization_rerun"=>false,"monte_carlo_rerun"=>false,
        "figure2_metadata_sha256"=>open(io->bytes2hex(sha256(io)),metadata))
    JLD2.jldopen(joinpath(output_dir,REACHING_ABLATION_RUN),"w";compress=true) do file
        for (key,value) in data;file[key]=value;end
    end
    return data
end
