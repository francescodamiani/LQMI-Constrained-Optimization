# Arm analysis definitions; traces use dimension,time,trial.
using MultivariateStats, StatsBase, Statistics, LinearAlgebra

# Helper functions 

# --- Log pseudodeterminant ---
log_pseudodet(A; tol=1e-12) = sum(log, svd(A).S[svd(A).S .> tol]; init=0.0)

function mean_log_pdet_over_time(tensor)
    T = size(tensor,3)
    return mean(log_pseudodet.(eachslice(tensor[:,:,2:T-2], dims=3)))
end

# --- EV1 of z ---
function compute_EV1(z_store)
    Nnoise = length(z_store)
    ev1 = zeros(Nnoise)

    for i in 1:Nnoise
        Zmat = reshape(z_store[i], size(z_store[i],1), :)
        pca_model = fit(PCA, Zmat)
        eigvals = principalvars(pca_model)
        ev1[i] = eigvals[1] / sum(eigvals)
    end
    return ev1
end

# --- Fraction relevant variance ---
function compute_relevant_variance_fraction(z_store; dimension_of_state=6)

    Nnoise = length(z_store)
    frac_rel = zeros(Nnoise)

    P_rel = zeros(dimension_of_state, dimension_of_state)
    P_rel[1:3,1:3] .= I(3)

    for i in 1:Nnoise
        Zmat = reshape(z_store[i], size(z_store[i],1), :)
        Zmat = Zmat .- mean(Zmat, dims=2)
        Σ = (Zmat * Zmat') / size(Zmat,2)

        total_var = tr(Σ)
        rel_var   = tr(P_rel * Σ)

        frac_rel[i] = rel_var / total_var
    end

    return frac_rel
end

# --- Peak velocity ---
function compute_peak_velocity_stats(x_store)
    Nnoise = length(x_store)
    mean_peak = zeros(Nnoise)
    sem_peak  = zeros(Nnoise)

    for i in 1:Nnoise
        x_vec = x_store[i]
        vel = x_vec[4:6,:,:]
        Ntrials = size(vel,3)
        trial_peaks = zeros(Ntrials)

        for trial in 1:Ntrials
            v_trial = vel[:,:,trial]
            v_norm = vec(norm.(eachcol(v_trial)))
            trial_peaks[trial] = maximum(v_norm)
        end

        mean_peak[i] = mean(trial_peaks)
        sem_peak[i]  = std(trial_peaks) / sqrt(Ntrials)
    end

    return mean_peak, sem_peak
end


function compute_global_pc1(Z::Array{Float64,3}; t_inds=nothing)
    # Z has shape: latent_dim × T × Ntrials
    latent_dim, Tfull, Ntrials = size(Z)

    if t_inds === nothing
        t_inds = collect(2:Tfull-2)
    end

    # pool trials and selected time points
    Zmat = reshape(Z[:, t_inds, :], latent_dim, :)
    pca_model = fit(PCA, Zmat)

    pc1 = projection(pca_model)[:, 1]
    pc1 ./= norm(pc1)

    return pc1
end

function compute_global_pc1_raw_readout_strength(z_store, L_all; t_inds=nothing)
    
    Nnoise = length(z_store)
    Tfull = size(z_store[1], 2)

    if t_inds === nothing
        t_inds = collect(2:Tfull-2)
    end

    strength_mean = zeros(Nnoise)
    strength_sem  = zeros(Nnoise)

    # optional: save the global PC1 for each noise level
    pc1_list = Vector{Vector{Float64}}(undef, Nnoise)

    for i in 1:Nnoise
        # one global PC1 for this internal-noise level
        v1 = compute_global_pc1(z_store[i]; t_inds=t_inds)
        pc1_list[i] = copy(v1)

        vals = Float64[]

        for t in t_inds
            Lt = L_all[i][:, :, t]
            r_t = dot(v1, (Lt' * Lt) * v1)  # = ||Lt * v1||^2
            push!(vals, r_t)
        end

        strength_mean[i] = mean(vals)
        strength_sem[i]  = std(vals) / sqrt(length(vals))
    end

    return strength_mean, strength_sem, pc1_list
end

