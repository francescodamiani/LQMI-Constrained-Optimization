# Source definitions use archived dimension,trial,time axes.
using LinearAlgebra, Statistics

log_pseudodet(A; tol=1e-6) = sum(log, filter(>(tol), svdvals(A)); init=0.0)

function mean_log_pdet_over_time(tensor)
    T_gain = size(tensor,3)
    return mean(log_pseudodet.(eachslice(tensor[:,:,2:T_gain-2], dims=3)))
end

function compute_peak_elbow_velocity(x)
    # Mean of per-trial peaks, not the peak of the mean trajectory or joint norm.
    peaks = rad2deg.(vec(maximum(abs.(x[5,:,:]), dims=2)))
    return (;peaks, mean_peak=mean(peaks), sem_peak=std(peaks)/sqrt(length(peaks)))
end

function compute_control_variability(z, L; time_indices=1:size(z,3)-2)
    N_trials = size(z,2)
    projected_variance = zeros(length(time_indices))
    internal_variance = zeros(length(time_indices))
    for (idx,t) in enumerate(time_indices)
        centered_z = z[:,:,t] .- mean(z[:,:,t], dims=2)
        covariance = centered_z * transpose(centered_z) / (N_trials-1)
        covariance = 0.5 .* (covariance .+ transpose(covariance))
        projected_variance[idx] = tr(L[:,:,t] * covariance * transpose(L[:,:,t]))
        internal_variance[idx] = tr(covariance)
    end
    V_u = sum(projected_variance)
    # Internal-state variance transmission per unit internal-state variance.
    variance_transmission_ratio = V_u / sum(internal_variance)
    return (;V_u, projected_variance, internal_variance, variance_transmission_ratio)
end

function normalize_control_variability(V_u, sigma_eta_list)
    length(V_u) == length(sigma_eta_list) || error("Control variance and noise grid have inconsistent lengths")
    baseline = findfirst(iszero, sigma_eta_list)
    baseline === nothing && error("Control-variance normalization requires a zero-internal-noise condition")
    isfinite(V_u[baseline]) && V_u[baseline] > 0 ||
        error("Zero-internal-noise control variance must be finite and positive")
    return V_u ./ V_u[baseline]
end

function fig2_statistics(results; representative_sigma_eta=0.019)
    sigma_eta_list = results["sigma_eta_list"]
    idx = findfirst(==(representative_sigma_eta), sigma_eta_list)
    idx === nothing && error("Representative internal noise $representative_sigma_eta is absent from the saved grid")
    statistics = Dict{String,Any}("representative_sigma_eta"=>representative_sigma_eta)
    for method in ["TOD","LCO"]
        X = results["x_store_$(method)"]
        statistics["endpoint_SD_pct_$(method)"] = [100*std(x[2,:,end])/(pi/9) for x in X]
        peaks = compute_peak_elbow_velocity.(X)
        statistics["mean_peak_$(method)"] = getproperty.(peaks,:mean_peak)
        statistics["sem_peak_$(method)"] = getproperty.(peaks,:sem_peak)
        statistics["gain_K_$(method)"] = mean_log_pdet_over_time.(results["K_$(method)_all"])
        statistics["gain_L_$(method)"] = mean_log_pdet_over_time.(results["L_$(method)_all"])
        transmission = [compute_control_variability(z,L) for (z,L) in
            zip(results["z_store_$(method)"],results["L_$(method)_all"])]
        statistics["V_u_$(method)"] = getproperty.(transmission,:V_u)
        statistics["V_u_relative_$(method)"] = normalize_control_variability(statistics["V_u_$(method)"],sigma_eta_list)
        statistics["variance_transmission_ratio_$(method)"] = getproperty.(transmission,:variance_transmission_ratio)
        statistics["variance_transmission_relative_$(method)"] =
            normalize_control_variability(statistics["variance_transmission_ratio_$(method)"],sigma_eta_list)
    end
    x = results["x_store_LCO"][idx]
    for (name,component) in [("angle",2),("velocity",5)]
        statistics["$(name)_mean"] = rad2deg.(vec(mean(x[component,:,:],dims=1)))
        statistics["$(name)_SD"] = rad2deg.(vec(std(x[component,:,:],dims=1)))
    end
    statistics["cost_improvement_pct"] = 100 .* (results["analytical_cost_TOD"] .-
        results["analytical_cost_LCO"]) ./ results["analytical_cost_TOD"]
    return statistics
end
