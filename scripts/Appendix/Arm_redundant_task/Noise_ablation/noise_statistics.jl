# Stored-data statistics only. Variances are averaged before taking square roots.
using LinearAlgebra, Statistics, DataFrames, JLD2

const NOISE_DATA_FILES = Dict(
    "behavioral"=>"behavioral_variability.jld2",
    "original_ablation"=>"original_fixed_policy_ablation.jld2",
    "policies"=>"original_full_model_policies.jld2",
    "high_noise"=>"external_noise_robustness.jld2")
const ARM_COMPONENT_NAMES = ["theta_1","theta_2","theta_3","omega_1","omega_2","omega_3"]
const NOISE_CHANNEL_LABELS = Dict("no_sensory_noise"=>"Sensory (both)",
    "no_process_noise"=>"Process","no_additive_sensory_noise"=>"Sensory add.",
    "no_motor_multiplicative_noise"=>"Motor mult.",
    "no_sensory_multiplicative_noise"=>"Sensory mult.","no_internal_noise"=>"Internal")

function component_sample_variances(traces,temporal_window)
    size(traces,3) > 1 || error("Sample variances require at least two trials")
    return vec(mean(var(traces[:,temporal_window,:];dims=3,corrected=true);dims=2))
end

function load_noise_experiments(data_dir,experiments)
    return Dict(name=>load_benchmark(data_dir,NOISE_DATA_FILES[name]) for name in experiments)
end

function component_SD_reduction(full_variance,ablated_variance)
    all(>(0),full_variance) || error("Full component variances must be positive")
    minimum(ablated_variance) >= -1e-10 || error("Negative ablated variance")
    full_SD = sqrt.(full_variance)
    return 100 .* (full_SD .- sqrt.(max.(ablated_variance,0.0))) ./ full_SD
end

function group_SD_reduction(full_variance,ablated_variance,components)
    # A ratio of mean component SDs, not a mean of component percentages or
    # sqrt of a mean across coordinates. Angle/velocity units stay separate.
    full_SD = mean(sqrt.(full_variance[components]))
    full_SD > 0 || error("Full mean component SD must be positive")
    return 100*(full_SD-mean(sqrt.(max.(ablated_variance[components],0.0))))/full_SD
end

function behavioral_tables(data;method="LCO")
    sigma_eta_list = data["sigma_eta_list"]
    baseline = findfirst(iszero,sigma_eta_list)
    baseline === nothing && error("Behavioral change requires zero internal noise")
    SD = sqrt.(data["variance_x_MC_$(method)"])
    rows = [merge((;sigma_eta=sigma_eta_list[i]),
        NamedTuple{Tuple(Symbol.(ARM_COMPONENT_NAMES))}(Tuple(SD[:,i]))) for i in eachindex(sigma_eta_list)]
    endpoints = DataFrame(component=ARM_COMPONENT_NAMES,
        baseline_SD=SD[:,baseline],final_SD=SD[:,end],
        change_pct=100 .* (SD[:,end] ./ SD[:,baseline] .- 1),
        maximum_change_pct=[maximum(100 .* (SD[j,:] ./ SD[j,baseline] .- 1)) for j in 1:6],
        sigma_at_maximum=[sigma_eta_list[argmax(SD[j,:])] for j in 1:6])
    return (;sweep=DataFrame(rows),endpoints)
end

function original_ablation_tables(data;method="LCO")
    method_idx = findfirst(==(method),data["method_names"])
    reductions = zeros(6,length(data["sigma_eta_list"]),length(data["ablation_names"]))
    groups = NamedTuple[]
    for (idx,sigma_eta) in enumerate(data["sigma_eta_list"])
        full = data["variance_x"][:,idx,method_idx,1]
        for a in eachindex(data["ablation_names"])
            reductions[:,idx,a] .= component_SD_reduction(full,data["variance_x"][:,idx,method_idx,a+1])
        end
        for (group,components) in [("angles",1:3),("velocities",4:6)]
            values = [group_SD_reduction(full,data["variance_x"][:,idx,method_idx,a+1],components)
                for a in eachindex(data["ablation_names"])]
            push!(groups,merge((;sigma_eta,group,full_mean_SD=mean(sqrt.(full[components]))),
                NamedTuple{Tuple(Symbol.(data["ablation_names"]))}(Tuple(values))))
        end
    end
    ranges = [(;component=ARM_COMPONENT_NAMES[j],ablation=name,
        minimum_reduction_pct=minimum(reductions[j,:,a]),maximum_reduction_pct=maximum(reductions[j,:,a]))
        for (a,name) in enumerate(data["ablation_names"]) for j in 1:6]
    return (;groups=DataFrame(groups),ranges=DataFrame(ranges),component_reductions=reductions)
end

function high_noise_tables(data;method="LCO")
    s_idx = findfirst(==(10.0),data["external_noise_multipliers"])
    s_idx === nothing && error("The s=10 noise regime is missing")
    method_idx = findfirst(==(method),data["method_names"])
    costs = DataFrame(sigma_eta=data["sigma_eta_list"],
        TOD_cost=data["analytical_cost"][1,s_idx,:],LCO_cost=data["analytical_cost"][2,s_idx,:])
    costs.improvement_pct = 100 .* (costs.TOD_cost .- costs.LCO_cost) ./ costs.TOD_cost
    rows = NamedTuple[]
    for (idx,sigma_eta) in enumerate(data["sigma_eta_list"])
        full = data["variance_x_analytical"][:,method_idx,s_idx,idx]
        for (group,components) in [("angles",1:3),("velocities",4:6)]
            values = [group_SD_reduction(full,data["ablation_variance_x"][:,a,method_idx,s_idx,idx],components)
                for a in eachindex(data["ablation_names"])]
            push!(rows,merge((;sigma_eta,group,full_mean_SD=mean(sqrt.(full[components]))),
                NamedTuple{Tuple(Symbol.(data["ablation_names"]))}(Tuple(values))))
        end
    end
    return (;costs,groups=DataFrame(rows))
end

function print_noise_experiments(results;io=stdout)
    if haskey(results,"behavioral")
        data = results["behavioral"]
        tables = behavioral_tables(data)
        println(io,"1. Behavioral variability — ",benchmark_run_label(data),"; LCO, ",data["N_trials_MC"]," trials/condition")
        println(io,"SD = sqrt(mean_t sample_var_trials(x_i,t)); window ",first(data["temporal_window"]),":",last(data["temporal_window"]),".")
        print_benchmark_table(io,tables.sweep,vcat(["sigma_eta"=>:sigma_eta],[name=>Symbol(name) for name in ARM_COMPONENT_NAMES]))
        println(io,"\nChanges relative to zero internal noise (final grid value and maximum over the sweep):")
        print_benchmark_table(io,tables.endpoints,["Component"=>:component,"SD at zero"=>:baseline_SD,
            "SD at final noise"=>:final_SD,"Final change (%)"=>:change_pct,
            "Max change (%)"=>:maximum_change_pct,"Noise at max"=>:sigma_at_maximum])
    end
    if haskey(results,"original_ablation")
        data = results["original_ablation"]
        tables = original_ablation_tables(data)
        println(io,"\n2. Baseline-noise fixed-policy ablation — ",benchmark_run_label(data),"; LCO")
        println(io,"Sensory removal zeros BOTH additive sensory covariance and multiplicative sensory D. No reoptimization.")
        println(io,"Component SD-reduction ranges across the complete internal-noise sweep:")
        ranges = copy(tables.ranges)
        ranges.ablation = [NOISE_CHANNEL_LABELS[name] for name in ranges.ablation]
        print_benchmark_table(io,ranges,["Component"=>:component,"Ablation"=>:ablation,
            "Minimum (%)"=>:minimum_reduction_pct,"Maximum (%)"=>:maximum_reduction_pct])
        println(io,"\nReductions in mean component SD, angles and velocities separately (not mean component percentages):")
        print_benchmark_table(io,tables.groups,vcat(["sigma_eta"=>:sigma_eta,"Group"=>:group,"Full mean SD"=>:full_mean_SD],
            [NOISE_CHANNEL_LABELS[name]*" (%)"=>Symbol(name) for name in data["ablation_names"]]))
    end
    if haskey(results,"high_noise")
        data = results["high_noise"]
        tables = high_noise_tables(data)
        println(io,"\n3. High sensory/motor regime — ",benchmark_run_label(data),"; s=10, process amplitude stays 0.1")
        println(io,"TOD recomputed and LCO initialized from TOD at each full-model internal-noise level:")
        # Show the archived cost comparison with full Float64 round-trip precision.
        costs = copy(tables.costs)
        for column in [:TOD_cost,:LCO_cost,:improvement_pct]
            costs[!,column] = string.(costs[!,column])
        end
        print_benchmark_table(io,costs,["sigma_eta"=>:sigma_eta,"TOD cost"=>:TOD_cost,
            "LCO cost"=>:LCO_cost,"Improvement (%)"=>:improvement_pct])
        println(io,"\nHigh-noise fixed-policy ablations — LCO mean component SD reductions:")
        print_benchmark_table(io,tables.groups,vcat(["sigma_eta"=>:sigma_eta,"Group"=>:group,"Full mean SD"=>:full_mean_SD],
            [NOISE_CHANNEL_LABELS[name]*" (%)"=>Symbol(name) for name in data["ablation_names"]]))
        println(io,"Additive and multiplicative sensory ablations are SEPARATE; their percentages are not an additive decomposition.")
    end
    println(io,"\nAngles are in rad, velocities in rad/s. SD reductions are fixed-closed-loop effects, not additive variance fractions.")
end
