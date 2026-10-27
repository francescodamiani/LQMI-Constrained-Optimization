# Redundant-arm internal-noise sweep and PCA plots.
using Plots, LaTeXStrings, MultivariateStats, StatsBase, Statistics, LinearAlgebra

function plot_internal_noise_sweep(results, plot_dir)
    (; sigma_eta_list, pc_var_list, pc_var_list_TOD, pc_dirs_list, pc_dirs_list_TOD,
       x_store, z_store, u_store, x_store_TOD, z_store_TOD, u_store_TOD,
       mean_cost_MM, sem_cost_MM, mean_cost_TOD, sem_cost_TOD,
       K_TOD_all, L_TOD_all, K_Full_Solver_all, L_Full_Solver_all,
       best_cost_list_all) = results
    dimension_of_state = results.problem_data.dimension_of_state
    dimension_of_control = results.problem_data.dimension_of_control
    T, Δt, S = results.problem_data.T, results.problem_data.Δt, results.problem_data.S
default(
    size = (550, 550),
    dpi = 500,
    thickness_scaling = 2.2,
    legend_background_color = :transparent,
    legend_foreground_color = :transparent, markerstrokewidth=0, box = false, grid = false
)

###################### (a) COST ############################

plt_cost = plot(
    sigma_eta_list, mean_cost_TOD,
    ribbon=sem_cost_TOD, lw=3, color = cgrad(:Blues_5, max(6, length(sigma_eta_list)))[6],
    label=L"$\mathrm{TOD}$",
    xlabel=L"$\sigma_{\eta}$", ylabel=L"$C$",
    yscale=:log10, grid=false
)

plot!(plt_cost,
    sigma_eta_list, mean_cost_MM,
    ribbon=sem_cost_MM, lw=3, color = cgrad(:Reds_5, max(6, length(sigma_eta_list)))[6],
    label=L"$\mathrm{LCO}$" #stands for Lagrangina Constrained Optimization
)

savefig(plt_cost, joinpath(plot_dir,"A_cost.pdf"))

################ (b) log |K|_p #############################

gain_K_LCO = mean_log_pdet_over_time.(K_Full_Solver_all)
gain_K_TOD = mean_log_pdet_over_time.(K_TOD_all)

plt_K = plot(
    sigma_eta_list, gain_K_TOD,
    lw=3, marker=:circle, markerstrokecolor=:transparent,
    color = cgrad(:Blues_5, max(6, length(sigma_eta_list)))[6], label=L"$\mathrm{TOD}$",
    xlabel=L"$\sigma_{\eta}$", ylabel=L"$\log(|K|_p)$",
    grid=false
)

plot!(plt_K,
    sigma_eta_list, gain_K_LCO,
    lw=3, marker=:circle, markerstrokecolor=:transparent,
    color = cgrad(:Reds_5, max(6, length(sigma_eta_list)))[6], label=L"$\mathrm{LCO}$"
)

savefig(plt_K, joinpath(plot_dir, "B_logK.pdf"))

################ (c) log |L|_p #############################

gain_L_LCO = mean_log_pdet_over_time.(L_Full_Solver_all)
gain_L_TOD = mean_log_pdet_over_time.(L_TOD_all)

plt_L = plot(
    sigma_eta_list, gain_L_TOD,
    lw=3, marker=:circle, markerstrokecolor=:transparent,
    color = cgrad(:Blues_5, max(6, length(sigma_eta_list)))[6], label=L"$\mathrm{TOD}$",
    xlabel=L"$\sigma_{\eta}$", ylabel=L"$\log(|L|_p)$",
    grid=false
)

# Display convention; raw gains and numerical data remain untouched.
gain_L_LCO[1] = gain_L_TOD[1]

plot!(plt_L,
    sigma_eta_list, gain_L_LCO,
    lw=3, marker=:circle, markerstrokecolor=:transparent,
    color= cgrad(:Reds_5, max(6, length(sigma_eta_list)))[6], label=L"$\mathrm{LCO}$"
)

savefig(plt_L, joinpath(plot_dir, "C_logL.pdf"))

################ (d) EV1 ###################################

ev1_LCO = compute_EV1(z_store)
ev1_TOD = compute_EV1(z_store_TOD)

plt_ev1 = plot(
    sigma_eta_list, ev1_TOD,
    lw=3, marker=:circle, markerstrokecolor=:transparent,
    color=cgrad(:Blues_5, max(6, length(sigma_eta_list)))[6], label=L"$\mathrm{TOD}$",
    xlabel=L"$\sigma_{\eta}$", ylabel=L"\mathrm{EV} - 1\mathrm{st PC\ of\ } z_t",
    grid=false
)

plot!(plt_ev1,
    sigma_eta_list, ev1_LCO,
    lw=3, marker=:circle, markerstrokecolor=:transparent,
    color=cgrad(:Reds_5, max(6, length(sigma_eta_list)))[6], label=L"$\mathrm{LCO}$"
)

# highlight first and last LCO points
# scatter!(plt_ev1,
#     [sigma_eta_list[1]], [ev1_LCO[1]],
#     markersize=7, marker=:circle,
#     color=:gray55, markerstrokecolor=:transparent,
#     label=""
# )

scatter!(plt_ev1,
    [sigma_eta_list[end]], [ev1_LCO[end]],
    markersize=7, marker=:circle,
    color=:darkorange3, markerstrokecolor=:transparent,
    label="", xlims = (-0.05,0.45), ylims = (0.25,1.05)
)
savefig(plt_ev1, joinpath(plot_dir, "D_EV1.pdf"))

xtick_labels = [
    L"$z_{\theta_1}$",
    L"$z_{\theta_2}$",
    L"$z_{\theta_3}$",
    L"$z_{\dot{\theta_1}}$",
    L"$z_{\dot{\theta_2}}$",
    L"$z_{\dot{\theta_3}}$"
]

################ (e.2) PC1 BAR — LCO HIGH NOISE #############

idx_high = length(sigma_eta_list)
Zmat = reshape(z_store[idx_high], dimension_of_state, :)
pca_model = fit(PCA, Zmat)
pc1 = projection(pca_model)[:,1]
pc1 ./= norm(pc1)

p_bar = bar(
    1:dimension_of_state, abs.(pc1),
    xticks=(1:dimension_of_state, xtick_labels),
    xlabel=L"$\mathrm{Component\ of}\ z_t$",
    ylabel=L"$1\mathrm{st PC\ \ weight}$",
    grid=false,
    legend=false,
    linecolor=:transparent,
    color=:darkorange3
)

savefig(p_bar, joinpath(plot_dir, "E_PC1_LCO_high_noise.pdf"))

# ################ (e.1) PC1 BAR — LCO LOW NOISE #############

# idx_low = 2
# Zmat = reshape(z_store[idx_low], dimension_of_state, :)
# pca_model = fit(PCA, Zmat)
# pc1 = projection(pca_model)[:,1]
# pc1 ./= norm(pc1)

# p_bar = bar(
#     1:dimension_of_state, abs.(pc1),
#     xticks=(1:dimension_of_state, xtick_labels),
#     xlabel=L"$\mathrm{Component\ of}\ z_t$",
#     ylabel=L"$1\mathrm{st PC\ \ weight}$",
#     grid=false,
#     legend=false,
#     linecolor=:transparent,
#     color=:gray55
# )

# savefig(p_bar, joinpath(plot_dir, "E_PC1_LCO_low_noise.pdf"))

# ################ (f) Relevant Variance Fraction ############

# frac_rel_LCO = compute_relevant_variance_fraction(z_store)
# frac_rel_TOD = compute_relevant_variance_fraction(z_store_TOD)

# plt_frac = plot(
#     sigma_eta_list, frac_rel_TOD,
#     lw=3, marker=:circle, markerstrokecolor=:transparent,
#     color=cgrad(:Blues_5, max(6, length(sigma_eta_list)))[6], label=L"$\mathrm{TOD}$",
#     xlabel=L"$\sigma_{\eta}$",
#     ylabel=L"$\mathrm{Tr}(P_{\mathrm{rel}}\Sigma_z) / \mathrm{Tr}(\Sigma_z)$",
#     grid=false
# )

# plot!(plt_frac,
#     sigma_eta_list, frac_rel_LCO,
#     lw=3, marker=:circle, markerstrokecolor=:transparent,
#     color=cgrad(:Reds_5, max(6, length(sigma_eta_list)))[6], label=L"$\mathrm{LCO}$", legend = :left
# )

# savefig(plt_frac, joinpath(plot_dir, "F_relevant_fraction.pdf"))

################ (g) Peak Velocity #########################

mean_peak_LCO, sem_peak_LCO = compute_peak_velocity_stats(x_store)
mean_peak_TOD, sem_peak_TOD = compute_peak_velocity_stats(x_store_TOD)

plt_vel = plot(
    sigma_eta_list, mean_peak_TOD,
    ribbon=sem_peak_TOD, lw=3, color=cgrad(:Blues_5, max(6, length(sigma_eta_list)))[6],
    label=L"$\mathrm{TOD}$",
    xlabel=L"$\sigma_{\eta}$",
    ylabel=L"$\mathrm{Peak\ Velocity}$",
    grid=false
)

plot!(plt_vel,
    sigma_eta_list, mean_peak_LCO,
    ribbon=sem_peak_LCO, lw=3, color=cgrad(:Reds_5, max(6, length(sigma_eta_list)))[6],
    label=L"$\mathrm{LCO}$"
)

savefig(plt_vel, joinpath(plot_dir, "G_peak_velocity.pdf"))

################ GLOBAL PC1(z) READOUT STRENGTH ACROSS INTERNAL NOISE (LCO ONLY) ################

# Here we compute, for each internal noise level:
#   1) one single global PC1 of z_t by pooling trials and time points in a movement window
#   2) for each time t in that same window, the readout strength of that PC1:
#          r(t) = ||L_t * v1_z||^2
#   3) the average of r(t) across time

# Same movement window used for the global PCA analysis
t_inds = 2:T-2

pc1_raw_readout_mean_LCO, pc1_raw_readout_sem_LCO, pc1_list_LCO =
    compute_global_pc1_raw_readout_strength(
        z_store,
        L_Full_Solver_all;
        t_inds = t_inds
    )

################ PLOT ################

plt_pc1_raw_readout_LCO = plot(
    sigma_eta_list,
    pc1_raw_readout_mean_LCO,
    ribbon = pc1_raw_readout_sem_LCO,
    lw = 3,
    marker = :circle,
    markerstrokecolor = :transparent,
    color = cgrad(:Reds_5, max(6, length(sigma_eta_list)))[6],
    label = "",
    xlabel = L"$\sigma_{\eta}$",
    ylabel = L"\langle \|L_t v_1^z\|^2 \rangle",
    grid = false
)


scatter!(plt_pc1_raw_readout_LCO,
    [sigma_eta_list[end]], [pc1_raw_readout_mean_LCO[end]],
    markersize=7, marker=:circle,
    color=:darkorange3, markerstrokecolor=:transparent,
    label="", xlims = (-0.03,0.45), ylims = (-0.1, 2.2)
)

savefig(plt_pc1_raw_readout_LCO, joinpath(plot_dir, "F_PC1_raw_readout_strength_LCO_globalPC.pdf"))

################ APPENDIX: TOD PCs #########################

# Zero internal noise
Zmat_TOD = reshape(z_store_TOD[1], dimension_of_state, :)
pca_TOD = fit(PCA, Zmat_TOD)
eigvals_TOD = principalvars(pca_TOD)
pcs_TOD = projection(pca_TOD)

for k in 1:length(eigvals_TOD)
    if eigvals_TOD[k] / sum(eigvals_TOD) > 1e-4
        pc_k = pcs_TOD[:,k]
        pc_k ./= norm(pc_k)

        p_bar_k = bar(
            1:dimension_of_state, abs.(pc_k),
            xticks=(1:dimension_of_state, xtick_labels),
            xlabel=L"$\mathrm{Component\ of}\ z_t$",
            ylabel=L"$\mathrm{PC%$k \ weight}$",
            grid=false, legend=false,
            linecolor=:transparent, color = :gray55
        )

        savefig(p_bar_k, joinpath(plot_dir, "Appendix_TOD_zero_noise_PC$(k).pdf"))
    end
end

# Build covariance
Zmat_TOD = reshape(z_store_TOD[1], dimension_of_state, :)
pca_TOD = fit(PCA, Zmat_TOD)

eigvals_TOD = principalvars(pca_TOD)
EV_TOD = eigvals_TOD ./ sum(eigvals_TOD)
pc_indices = 1:length(EV_TOD)

plt_EV_spectrum = plot(
    pc_indices,
    EV_TOD,
    lw = 2,
    marker = :circle,
    markerstrokecolor = :transparent,
    color = :gray55,
    xlabel = L"\mathrm{PC\ index}",
    ylabel = L"\mathrm{EV}",
    grid = false, size = (1100, 550), label = L"$\sigma_{\eta} = 0.0$"
)

savefig(
    plt_EV_spectrum,
    joinpath(plot_dir, "Appendix_TOD_zero_noise_EV_spectrum.pdf")
)

# High internal noise

# Highest internal noise index
idx_high = length(sigma_eta_list)

Zmat_TOD = reshape(z_store_TOD[idx_high], dimension_of_state, :)
pca_TOD = fit(PCA, Zmat_TOD)
eigvals_TOD = principalvars(pca_TOD)
pcs_TOD = projection(pca_TOD)

for k in 1:length(eigvals_TOD)
    if eigvals_TOD[k] / sum(eigvals_TOD) > 1e-4
        pc_k = pcs_TOD[:,k]
        pc_k ./= norm(pc_k)

        p_bar_k = bar(
            1:dimension_of_state, abs.(pc_k),
            xticks=(1:dimension_of_state, xtick_labels),
            xlabel=L"$\mathrm{Component\ of}\ z_t$",
            ylabel=L"$\mathrm{PC%$k \ weight}$",
            grid=false, legend=false,
            linecolor=:transparent, color = cgrad(:Blues_5, max(6, length(sigma_eta_list)))[6]
        )

        savefig(p_bar_k, joinpath(plot_dir, "Appendix_TOD_high_noise_PC$(k).pdf"))
    end
end

# Highest internal noise index
idx_high = length(sigma_eta_list)

# Build covariance
Zmat_TOD = reshape(z_store_TOD[idx_high], dimension_of_state, :)
pca_TOD = fit(PCA, Zmat_TOD)

eigvals_TOD = principalvars(pca_TOD)
EV_TOD = eigvals_TOD ./ sum(eigvals_TOD)
pc_indices = 1:length(EV_TOD)

plt_EV_spectrum = plot(
    pc_indices,
    EV_TOD,
    lw = 2,
    marker = :circle,
    markerstrokecolor = :transparent,
    color = cgrad(:Blues_5, max(6, length(sigma_eta_list)))[6],
    xlabel = L"\mathrm{PC\ index}",
    ylabel = L"\mathrm{EV}",
    grid = false, size = (1100, 550), label = L"$\sigma_{\eta} = 0.4$"
)

savefig(
    plt_EV_spectrum,
    joinpath(plot_dir, "Appendix_TOD_high_noise_EV_spectrum.pdf")
)

return nothing
end
