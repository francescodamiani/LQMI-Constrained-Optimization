# Reaching-identifiability matching statistics and numerical tables; no plots.
function interpolate_internal_result(
    target_SD, internal_SD, internal_values
)

    tolerance = 1e-10

    if target_SD < minimum(internal_SD)-tolerance ||
       target_SD > maximum(internal_SD)+tolerance
        return NaN, false
    end

    exact_idx =
        argmin(abs.(internal_SD .- target_SD))

    if isapprox(
        internal_SD[exact_idx],
        target_SD;
        atol=1e-12,
        rtol=1e-10
    )
        return internal_values[exact_idx], true
    end

    for idx in 1:length(internal_SD)-1

        SD_1 = internal_SD[idx]
        SD_2 = internal_SD[idx+1]

        if (SD_1-target_SD)*(SD_2-target_SD) <= 0.0

            if isapprox(
                SD_1, SD_2;
                atol=1e-14,
                rtol=1e-12
            )
                return (
                    internal_values[idx] +
                    internal_values[idx+1]
                )/2, true
            end

            weight =
                (target_SD-SD_1)/(SD_2-SD_1)

            interpolated_value =
                internal_values[idx] +
                weight*(
                    internal_values[idx+1] -
                    internal_values[idx]
                )

            return interpolated_value, true

        end
    end

    return NaN, false

end

mean_gain_norm(gain) = mean(norm(gain[:,:,t]) for t in 2:size(gain,3)-2)
function mean_gain_logvolume(gain; tolerance=1e-4)
    return mean(sum(log,filter(>(tolerance),svdvals(gain[:,:,t]));init=0.0)
        for t in 2:size(gain,3)-2)
end

function overestimation_table(d)
    SD=d["terminal_SD_analytical_internal_percent"]
    any(diff(SD).< -1e-10) && error("Internal endpoint variability is not monotonic")
    base=argmin(abs.(d["sigma_omega_hat_list"].-d["sigma_omega_true"]))
    return DataFrame([begin
        target=d["terminal_SD_analytical_true_percent"][i]
        eta,matched=interpolate_internal_result(target,SD,d["sigma_eta_list"])
        internal_cost,_=interpolate_internal_result(target,SD,d["analytical_cost_internal"])
        internal_peak,_=interpolate_internal_result(target,SD,d["mean_peak_velocity_internal_deg_s"])
        (;assumed=d["sigma_omega_hat_list"][i],true_noise=d["sigma_omega_true"],
            cost=d["true_analytical_cost"][i],assumed_cost=d["assumed_model_cost"][i],
            change=100*(d["true_analytical_cost"][i]/d["true_analytical_cost"][base]-1),
            SD_percent=target,MC_SD=d["terminal_SD_MC_true_deg"][i],
            peak=d["mean_peak_velocity_wrong_deg_s"][i],SEM=d["sem_peak_velocity_wrong_deg_s"][i],
            peak_time=d["mean_peak_time_wrong_s"][i],K=d["K_Frobenius_wrong"][i],L=d["L_Frobenius_wrong"][i],
            logK=d["K_log_pdet_wrong"][i],logL=d["L_log_pdet_wrong"][i],
            matched,eta,internal_cost,internal_peak)
    end for i in eachindex(d["sigma_omega_hat_list"])])
end

function underestimation_table(d,branch)
    wrong=branch=="wrong"
    noise=wrong ? d["sigma_omega_hat_list"] : d["sigma_eta_list"]
    cost=d[wrong ? "wrong_true_cost" : "internal_cost"]
    return DataFrame([(;noise=noise[i],cost=cost[i],change=100*(cost[i]/cost[1]-1),
        SD_percent=5*d[branch*"_SD_analytical_deg"][i],
        MC_SD=d[branch*"_SD_MC_deg"][i],peak=d[branch*"_peak_velocity"][i],
        SEM=d[branch*"_peak_velocity_SEM"][i],peak_time=d[branch*"_peak_time"][i],
        K=d[branch*"_K_norm"][i],L=d[branch*"_L_norm"][i],
        logK=d[branch*"_K_log_pdet"][i],logL=d[branch*"_L_log_pdet"][i])
        for i in eachindex(noise)])
end

function noise_source_table(d,reference=nothing)
    return DataFrame([(;source=d["source_names"][j],eta=d["sigma_eta_list"][i],
        amplitude=d["matched_noise_values"][i,j],target=d["reference_terminal_SD_analytical_deg"][i],
        SD=d["analytical_terminal_SD_deg"][i,j],
        match_error=100*(d["analytical_terminal_SD_deg"][i,j]/d["reference_terminal_SD_analytical_deg"][i]-1),
        cost=d["analytical_cost"][i,j],peak=d["mean_peak_velocity_deg_s"][i,j],
        SEM=d["sem_peak_velocity_deg_s"][i,j],K=d["K_Frobenius"][i,j],L=d["L_Frobenius"][i,j],
        logK=j==1 ? (isnothing(reference) ? NaN : mean_gain_logvolume(reference["K_LCO_all"][i];tolerance=1e-6)) :
            mean_gain_logvolume(d["K_LCO_matched"][:,:,:,i,j-1];tolerance=1e-6),
        logL=j==1 ? (isnothing(reference) ? NaN : mean_gain_logvolume(reference["L_LCO_all"][i];tolerance=1e-6)) :
            mean_gain_logvolume(d["L_LCO_matched"][:,:,:,i,j-1];tolerance=1e-6))
        for j in eachindex(d["source_names"]) for i in eachindex(d["sigma_eta_list"])])
end

function print_identifiability(data; io=stdout, full=false)
    println(io,"Internal-noise comparisons within the tested reaching model and parameter ranges.")
    println(io,"These results do not establish general identifiability from behavior alone.")
    if haskey(data,"overestimation")
        d=data["overestimation"];table=overestimation_table(d)
        println(io,"\nA: sensory-model mismatch — ",benchmark_run_label(d),"; true sigma_omega=",d["sigma_omega_true"],"; eta=0; N_MC=",d["N_trials_MC"])
        full || (table=table[table.assumed .>= d["sigma_omega_true"],:])
        print_benchmark_table(io,table,["Assumed omega"=>:assumed,"True cost"=>:cost,
            "Cost change (%)"=>:change,"Endpoint SD (%20deg)"=>:SD_percent,
            "Peak (deg/s)"=>:peak,"Peak SEM"=>:SEM,"Mean ‖K‖F"=>:K,"Mean ‖L‖F"=>:L,
            "Mean log pdet(K)"=>:logK,"Mean log pdet(L)"=>:logL,
            "Overlap"=>:matched,"Matched eta"=>:eta,"Matched internal peak"=>:internal_peak])
        println(io,"Gains are frozen for true-system evaluation. Matching interpolates archived internal statistics only within endpoint-SD overlap; no extrapolation.")
    end
    if haskey(data,"noise_sources")
        d=data["noise_sources"];table=noise_source_table(d,get(data,"reference",nothing))
        println(io,"\nB: separately reoptimized noise sources — ",benchmark_run_label(d),"; N_MC=",d["N_trials_MC"])
        if full
            print_benchmark_table(io,table,["Source"=>:source,"Reference eta"=>:eta,"Amplitude"=>:amplitude,
                "Target SD (deg)"=>:target,"Achieved SD (deg)"=>:SD,"SD error (%)"=>:match_error,
                "Cost"=>:cost,"Peak (deg/s)"=>:peak,"Peak SEM"=>:SEM,"Mean ‖K‖F"=>:K,"Mean ‖L‖F"=>:L,
                "Mean log pdet(K)"=>:logK,"Mean log pdet(L)"=>:logL])
        else
            summary=DataFrame([begin
                g=table[table.source .== source,:]
                (;source,first=g.peak[1],maximum=maximum(g.peak),last=g.peak[end],
                    min_SD=minimum(5*g.SD),max_SD=maximum(5*g.SD),
                    max_match_error=maximum(abs.(g.match_error)),unmatched=count(>(1.0),abs.(g.match_error)),
                    rises=count(>(0),diff(g.peak)),falls=count(<(0),diff(g.peak)))
            end for source in d["source_names"]])
            print_benchmark_table(io,summary,["Source"=>:source,"First peak"=>:first,"Max peak"=>:maximum,
                "Last peak"=>:last,"Min SD (%20deg)"=>:min_SD,"Max SD (%20deg)"=>:max_SD,
                "Max abs SD error (%)"=>:max_match_error,"SD error >1%"=>:unmatched,"Rises"=>:rises,"Falls"=>:falls])
        end
        println(io,"Candidate-grid crossings interpolate noise amplitude; gains are then reoptimized at that amplitude. SD errors report residual matching error after reoptimization.")
        if any(abs.(table.match_error).>1.0)
            println(io,"LIMITATION: some candidate grids do not bracket the targets. Closest-grid fallback conditions with >1% SD error are not endpoint matched.")
        end
        println(io,"Rises/falls count sampled adjacent changes; small Monte Carlo reversals do not establish a systematic non-monotonic effect.")
    end
    if haskey(data,"underestimation")
        d=data["underestimation"]
        for branch in ["wrong","internal"]
            println(io,"\nC: ",branch=="wrong" ? "sensory underestimation (frozen gains, eta=0)" : "correct sensory model with internal noise",
                " — ",benchmark_run_label(d),"; true omega=",d["sigma_omega_true"],"; N_MC=",d["N_trials_MC"])
            print_benchmark_table(io,underestimation_table(d,branch),[
                (branch=="wrong" ? "Assumed omega" : "eta")=>:noise,"True cost"=>:cost,
                "Cost change (%)"=>:change,"Endpoint SD (%20deg)"=>:SD_percent,
                "Peak (deg/s)"=>:peak,"Peak SEM"=>:SEM,"Mean log pdet(K)"=>:logK,"Mean log pdet(L)"=>:logL])
        end
        println(io,"C uses raw sweeps, without endpoint matching. Gain volumes average log singular-value products over t=2:",d["T"]-3,", threshold ",d["pseudodeterminant_tolerance"],"; changes are not entrywise.")
    end
    println(io,"\nCosts and trial statistics come from the selected archives. Tables show 12 significant digits; JLD2 retains full precision. --full includes A's underestimation branch and every B condition.")
end

