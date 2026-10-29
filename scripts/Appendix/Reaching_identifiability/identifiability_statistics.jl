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

function print_identifiability(data; io=stdout, full=false)
    println(io,"Internal-noise comparisons within the tested reaching model and parameter ranges.")
    println(io,"These results do not establish general identifiability from behavior alone.")
    if haskey(data,"overestimation")
        d=data["overestimation"];table=overestimation_table(d)
        println(io,"\nSensory-model mismatch — ",benchmark_run_label(d),"; true sigma_omega=",d["sigma_omega_true"],"; eta=0; N_MC=",d["N_trials_MC"])
        full || (table=table[table.assumed .>= d["sigma_omega_true"],:])
        print_benchmark_table(io,table,["Assumed omega"=>:assumed,"True cost"=>:cost,
            "Cost change (%)"=>:change,"Endpoint SD (%20deg)"=>:SD_percent,
            "Peak (deg/s)"=>:peak,"Peak SEM"=>:SEM,"Mean ‖K‖F"=>:K,"Mean ‖L‖F"=>:L,
            "Mean log pdet(K)"=>:logK,"Mean log pdet(L)"=>:logL,
            "Overlap"=>:matched,"Matched eta"=>:eta,"Matched internal peak"=>:internal_peak])
        println(io,"Gains are frozen for true-system evaluation. Matching interpolates saved internal statistics only within endpoint-SD overlap; no extrapolation.")
    end
    if haskey(data,"underestimation")
        d=data["underestimation"]
        for branch in ["wrong","internal"]
            println(io,"\n",branch=="wrong" ? "sensory underestimation (frozen gains, eta=0)" : "correct sensory model with internal noise",
                " — ",benchmark_run_label(d),"; true omega=",d["sigma_omega_true"],"; N_MC=",d["N_trials_MC"])
            print_benchmark_table(io,underestimation_table(d,branch),[
                (branch=="wrong" ? "Assumed omega" : "eta")=>:noise,"True cost"=>:cost,
                "Cost change (%)"=>:change,"Endpoint SD (%20deg)"=>:SD_percent,
                "Peak (deg/s)"=>:peak,"Peak SEM"=>:SEM,"Mean log pdet(K)"=>:logK,"Mean log pdet(L)"=>:logL])
        end
        println(io,"Underestimation uses raw sweeps, without endpoint matching. Gain volumes average log singular-value products over t=2:",d["T"]-3,", threshold ",d["pseudodeterminant_tolerance"],"; changes are not entrywise.")
    end
    println(io,"\nCosts and trial statistics come from the saved results. Tables show 12 significant digits; JLD2 retains full precision. --full includes the lower half of the sensory-model mismatch grid.")
end

