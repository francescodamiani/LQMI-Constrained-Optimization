# Figure 1b–g: statistical definitions and PDF plotting.
using DataFrames, Statistics, LinearAlgebra, JLD2, Plots, LaTeXStrings, Printf

const FIG1_DIMENSIONS = "Fig1_ab_data_10000_GD.jld2"
const FIG1_NOISE = "Fig_internal_noise_sweep_full_solutions_IC_at_LCO_solutions_10000_GD.jld2"
const FIG1_CONVERGENCE = "Fig1d_LCO_cost_trajectories.jld2"

function load_fig1(data_dir, panels=collect("bcdefg"))
    files = String[]
    any(p -> p in "bc", panels) && push!(files, FIG1_DIMENSIONS)
    'd' in panels && push!(files, FIG1_CONVERGENCE)
    any(p -> p in "efg", panels) && push!(files, FIG1_NOISE)
    missing = filter(f -> !isfile(joinpath(data_dir, f)), files)
    isempty(missing) || error("Missing numerical data in $data_dir: $(join(missing, ", ")). " *
        "Restore the saved datasets or explicitly use --simulate to generate new numerical data.")
    return Dict(f => JLD2.load(joinpath(data_dir, f)) for f in files)
end

function agg_mean_sem(df::DataFrame, col::Symbol)
    gdf = combine(groupby(df, [:n, :method]),
        col => mean => Symbol(col, :_mean),
        col => (x -> (length(x) > 1 ? std(x) / sqrt(length(x)) : 0.0)) => Symbol(col, :_sem),
        col => length => :N)
    sort!(gdf, [:method, :n])
    return gdf
end

function gain_distances(solutions_all)
    Tm1, Nσ = size(solutions_all["L_GD"])[3:4]
    distances = Dict{String,Matrix{Float64}}()
    # norm(matrix) is the Frobenius norm. Average the norms over ALL T-1 slices.
    for gain in ["K", "L"], method in ["TOD", "LCO"]
        dist = zeros(Tm1, Nσ)
        for iσ in 1:Nσ, t in 1:Tm1
            dist[t, iσ] = norm(solutions_all["$(gain)_$(method)"][:, :, t, iσ] -
                               solutions_all["$(gain)_GD"][:, :, t, iσ])
        end
        distances["$(gain)_$(method)"] = dist
    end
    return distances
end

function plot_fig1(datasets, plot_dir; panels=collect("bcdefg"))
    any(d -> get(get(d, "settings", Dict()), "smoke", false), values(datasets)) &&
        error("Smoke datasets are for numerical validation only; no figures are generated from them.")
    default(size=(550, 550), dpi=500, thickness_scaling=2.2,
        legend_background_color=:transparent, legend_foreground_color=:transparent,
        markerstrokewidth=0, box=false, grid=false)
    colors = Dict("TOD"=>cgrad(:Blues_5, 20)[6], "LCO"=>cgrad(:Reds_5, 20)[6],
                  "GD"=>cgrad(:BuGn_5, 6)[5])
    plots = Dict{Char,Any}()
    if any(p -> p in "bc", panels)
        df_min = datasets[FIG1_DIMENSIONS]["df_min"]
        if 'b' in panels
            gdf_cost = agg_mean_sem(df_min, :cost)
            p_cost = plot(xaxis=(L"n_{\mathrm{task}}", :log10), yaxis=(L"C", :log10), legend=:topleft)
            for method in ["TOD", "LCO", "GD"]
                sub = gdf_cost[gdf_cost.method .== method, :]
                plot!(p_cost, sub.n, sub.cost_mean; ribbon=sub.cost_sem, lw=3,
                    marker=:circle, markerstrokecolor=:transparent, color=colors[method], label=latexstring("\\mathrm{" * method * "}"))
            end
            plots['b'] = p_cost
        end
        if 'c' in panels
            # Runtime-only exclusion: n=1,seed=1 for ALL methods.
            # The full archive and panel-b statistics retain those three rows.
            df_time = df_min[.!((df_min.n .== 1) .& (df_min.seed .== 1)), :]
            gdf_time = agg_mean_sem(df_time, :time_sec)
            p_time = plot(xaxis=(L"n_{\mathrm{task}}", :log10),
                yaxis=(L"\mathrm{runtime\ (s)}", :identity),
                yticks=([1.0, 1e4, 3e4], ["10^{0}", "10^{4}", "3 · 10^{4}"]),
                legend=:topright, ylim=(-2000, 3.3e4), xlim=(0.7, 200))
            for (method, markersize) in zip(["TOD", "LCO", "GD"], [7, 5, 6])
                sub = gdf_time[gdf_time.method .== method, :]
                plot!(p_time, sub.n, sub.time_sec_mean; lw=3, marker=:circle,
                    markersize, markerstrokecolor=:transparent, color=colors[method], label=latexstring("\\mathrm{" * method * "}"))
            end
            # Runtime ribbons are not plotted.
            plots['c'] = p_time
        end
    end
    if 'd' in panels
        cost_trajectories = datasets[FIG1_CONVERGENCE]["cost_trajectories"]
        sigma_eta_list = cost_trajectories["sigma_eta_list"]
        colors_sigma = cgrad(:Reds_7, 12)
        p = plot(xaxis=L"\mathrm{Iteration}", yaxis=L"C", legend=:topright, xlims=(0.8, 10))
        for (iσ, σ_η) in enumerate(sigma_eta_list)
            tr = cost_trajectories["traj"][iσ]
            # Initial cost at index 1, then one entry per K/L update; preserve all entries.
            plot!(p, eachindex(tr), tr; lw=3, color=colors_sigma[iσ],
                label=latexstring("\\sigma_{\\eta} = " * @sprintf("%.2f", σ_η)))
        end
        plots['d'] = p
    end
    if any(p -> p in "efg", panels)
        sweep = datasets[FIG1_NOISE]
        df_sweep, solutions_all, sigma_eta_list = sweep["df_sweep"], sweep["solutions_all"], sweep["sigma_eta_list"]
        if 'e' in panels
            p_cost = plot(xaxis=L"\sigma_{\eta}", yaxis=L"C", legend=:topleft,
                xticks=[0.0, 0.35, 0.7], xlims=(-0.05, 0.75), ylims=(4, 20))
            for (method, markersize) in zip(["TOD", "LCO", "GD"], [6, 5, 4])
                sub = sort(df_sweep[df_sweep.method .== method, :], :sigma_eta)
                plot!(p_cost, sub.sigma_eta, sub.mc_cost_mean; ribbon=sub.mc_cost_sem,
                    lw=3, marker=:circle, markersize, markerstrokecolor=:transparent,
                    color=colors[method], label=latexstring("\\mathrm{" * method * "}"), alpha=(method == "GD" ? 0.8 : 1.0))
            end
            plots['e'] = p_cost
        end
        if any(p -> p in "fg", panels)
            dist = gain_distances(solutions_all)
            for (panel, gain, limits, legend) in [('f', "K", (-5, 50), :topleft), ('g', "L", (-0.1, 0.9), :right)]
                panel in panels || continue
                ylabel = latexstring("\\langle\\|" * gain * "_t-" * gain * "_{t,\\mathrm{GD}}\\|_{\\mathrm{F}}\\rangle")
                p = plot(xaxis=L"\sigma_{\eta}", yaxis=ylabel, legend=legend,
                    ylims=limits, xlims=(-0.05, 0.75), guidefontsize=9, tickfontsize=7)
                for method in ["TOD", "LCO"]
                    plot!(p, sigma_eta_list, vec(mean(dist["$(gain)_$(method)"], dims=1));
                        color=colors[method], marker=:circle, markersize=6, lw=3, label=latexstring("\\mathrm{" * method * "}"))
                end
                plots[panel] = p
            end
        end
    end
    for panel in panels
        # Relative annotations add manuscript letters without changing axis limits.
        plot!(plots[panel]; top_margin=6Plots.mm)
        annotate!(plots[panel], (0.0, 1.10), Plots.text("($panel)", 8, :left))
        savefig(plots[panel], joinpath(plot_dir, "fig1$(panel).pdf"))
    end
    if panels == collect("bcdefg")
        # Manuscript quantitative-panel order; panel a artwork is outside this task.
        composite = plot([plots[p] for p in "bcdefg"]...; layout=(2, 3), size=(1650, 1100),
            dpi=300, thickness_scaling=1.2, left_margin=8Plots.mm, top_margin=12Plots.mm)
        savefig(composite, joinpath(plot_dir, "fig1_b-g.pdf"))
    end
    return plots
end
