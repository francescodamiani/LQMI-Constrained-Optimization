# LCO relative-perturbation violins and median overlays.
using StatsPlots, LaTeXStrings, JLD2, DataFrames, CodecZlib

const PERTURBATION_DATA_FILENAME = "LCO_local_perturbation_scan_sigma_eta_0p200.jld2"

function load_perturbations(data_dir)
    file = joinpath(data_dir, PERTURBATION_DATA_FILENAME)
    isfile(file) && return JLD2.load(file)
    # The published archive is gzip-wrapped without changing any original bytes.
    isfile(file * ".gz") || error("Missing perturbation dataset: $file[.gz]. Use --simulate to recompute it explicitly.")
    return mktemp() do path, output
        open(file * ".gz") do compressed
            input = GzipDecompressorStream(compressed)
            try
                while !eof(input)
                    write(output, read(input, 1024^2))
                end
            finally
                close(input)
            end
        end
        close(output)
        JLD2.load(path)
    end
end

function plot_perturbations(results, folder_out)
    if haskey(results, "settings") && results["settings"].validation_only
        error("Reduced numerical validation data cannot produce paper figures.")
    end
    df_scan = results["df_scan"]
    df_summary = results["df_summary"]
    epsilon_list = results["epsilon_list"]
    all(isfinite, df_scan.delta_cost) && all(>(0.0), df_scan.delta_cost) ||
        error("Log-scale violins require finite positive delta_cost values; inspect the numerical results before plotting.")
# PLOT DEFAULTS --------------------------------------------------------------------------------------------------------
default(
    size = (550, 550),
    dpi = 500,
    thickness_scaling = 2.2,
    legend_background_color = :transparent,
    legend_foreground_color = :transparent,
    markerstrokewidth = 0,
    box = false,
    grid = false
)

epsilon_labels = [epsilon_to_label(ε) for ε in epsilon_list]

# optional: common y-limits across the three plots, so they match nicely in Inkscape
all_vals = df_scan.delta_cost[isfinite.(df_scan.delta_cost)]
ymin = minimum(all_vals)
ymax = maximum(all_vals)
ypad = 0.05 * (ymax - ymin + 1e-12)
common_ylims = (ymin - ypad, ymax + ypad)

color_L     = cgrad(:Blues_5, 6)[5]
color_K     = cgrad(:Greens_5, 6)[5]
color_joint = cgrad(:Purples_5, 6)[5]

color_L     = cgrad(:Reds_5, 6)[5]
color_K     = cgrad(:Reds_5, 6)[5]
color_joint     = cgrad(:Reds_5, 6)[5]


# ==================================================================================================
# (1) L only
# ==================================================================================================
sub_L = df_scan[df_scan.perturbation_type .== "L_only", :]
sub_sum_L = sort!(df_summary[df_summary.perturbation_type .== "L_only", :], :epsilon)

p_L = plot(
    xaxis = (L"$\mathrm{Pert \ amplitude}$", :identity),
    yaxis = (L"$\Delta C$", :log10),
    title = "",
    legend = :topleft,
    xticks = (1:length(epsilon_list), epsilon_labels),
    xrotation = 45,
    #ylims = common_ylims, yticks = [0,300,600]
)

for (i, ε) in enumerate(epsilon_list)
    vals = sub_L.delta_cost[sub_L.epsilon .== ε]

    violin!(
        p_L,
        fill(i, length(vals)),
        vals,
        fillcolor = color_L,
        fillalpha = 0.65,
        linecolor = :gray35,
        linealpha = 0.6,
        linewidth = 1.4, 
        label = false
    )
end

plot!(
    p_L,
    1:length(epsilon_list),
    sub_sum_L.median_delta_cost,
    lw = 2.5,
    marker = :circle,
    linecolor = :gray35,
    linealpha = 0.6,
    markersize = 4.5,
    markercolor = :gray35,
    markeralpha = 0.6,
    label = L"$L \ \mathrm{only}$"
)

#hline!(p_L, [0.0], lw = 2, ls = :dash, color = :gray, alpha = 0.8)

savefig(p_L, joinpath(folder_out, "Fig_LCO_local_perturbation_violin_deltaC_L_only.pdf"))

# ==================================================================================================
# (2) K only
# ==================================================================================================
sub_K = df_scan[df_scan.perturbation_type .== "K_only", :]
sub_sum_K = sort!(df_summary[df_summary.perturbation_type .== "K_only", :], :epsilon)

p_K = plot(
    xaxis = (L"$\mathrm{Pert \ amplitude}$", :identity),
    yaxis = (L"$\Delta C$", :log10),
    title = "",
    legend = :topleft,
    xticks = (1:length(epsilon_list), epsilon_labels),
    xrotation = 45,
    #ylims = common_ylims, yticks = [0,300,600]
)

for (i, ε) in enumerate(epsilon_list)
    vals = sub_K.delta_cost[sub_K.epsilon .== ε]

    violin!(
        p_K,
        fill(i, length(vals)),
        vals,
        fillcolor = color_K,
        fillalpha = 0.65,
        linecolor = :gray35,
        linealpha = 0.6,
        linewidth = 1.4,
        label = false
    )
end

plot!(
    p_K,
    1:length(epsilon_list),
    sub_sum_K.median_delta_cost,
    lw = 2.5,
    marker = :circle,
    linecolor = :gray35,
    linealpha = 0.6,
    markersize = 4.5,
    markercolor = :gray35,
    markeralpha = 0.6,
    label = L"$K \ \mathrm{only}$"

)

#hline!(p_K, [0.0], lw = 2, ls = :dash, color = :gray, alpha = 0.8)

savefig(p_K, joinpath(folder_out, "Fig_LCO_local_perturbation_violin_deltaC_K_only.pdf"))

# ==================================================================================================
# (3) Joint (L,K)
# ==================================================================================================
sub_joint = df_scan[df_scan.perturbation_type .== "joint", :]
sub_sum_joint = sort!(df_summary[df_summary.perturbation_type .== "joint", :], :epsilon)

p_joint = plot(
    xaxis = (L"$\mathrm{Pert \ amplitude}$", :identity),
    yaxis = (L"$\Delta C$", :log10),
    title = "",
    legend = :topleft,
    xticks = (1:length(epsilon_list), epsilon_labels),
    xrotation = 45,
    #ylims = common_ylims, yticks = [0,300,600]
)

for (i, ε) in enumerate(epsilon_list)
    vals = sub_joint.delta_cost[sub_joint.epsilon .== ε]

    violin!(
        p_joint,
        fill(i, length(vals)),
        vals,
        fillcolor = color_joint,
        fillalpha = 0.65,
        linecolor = :gray35,
        linealpha = 0.6,
        linewidth = 1.4,
        label = false
    )
end

plot!(
    p_joint,
    1:length(epsilon_list),
    sub_sum_joint.median_delta_cost,
    lw = 2.5,
    marker = :circle,
    linecolor = :gray35,
    linealpha = 0.6,
    markersize = 4.5,
    markercolor = :gray35,
    markeralpha = 0.6,
    label = L"$\mathrm{Joint} \ (L,K)$"

)

#hline!(p_joint, [0.0], lw = 2, ls = :dash, color = :gray, alpha = 0.8)

savefig(p_joint, joinpath(folder_out, "Fig_LCO_local_perturbation_violin_deltaC_joint.pdf"))


return (p_L, p_K, p_joint)
end
