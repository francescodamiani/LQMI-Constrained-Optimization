# Eight individual Figure 2 panels. Numerical definitions are in fig2_statistics.jl.
using Plots, LaTeXStrings
using Plots.PlotMeasures

const FIG2_PDFS = ["A_elbow_angle.pdf","B_endpoint_variability.pdf","C_elbow_velocity.pdf",
    "D_expected_cost.pdf","E_logK.pdf","F_logL.pdf","G_control_variability.pdf","H_peak_velocity.pdf"]

function plot_fig2(results,plot_dir; representative_sigma_eta=0.019,panels=collect("ABCDEFGH"))
    !isempty(panels) && all(p -> p in "ABCDEFGH",panels) || error("Select Figure 2 panels A–H")
    outputs = [FIG2_PDFS[findfirst(==(panel),"ABCDEFGH")] for panel in panels]
    get(results,"validation_only",false) && error("Reduced numerical validation data cannot produce paper figures.")
    s = fig2_statistics(results;representative_sigma_eta)
    sigma_eta_list = results["sigma_eta_list"]
    time_vector = results["time_vector"]
    # Explicit fonts keep the compact panels readable when placed four across,
    # avoiding GR's excessive padding with global thickness_scaling=2.
    default(size=(400,400),dpi=500,thickness_scaling=1.0,
        guidefontsize=26,tickfontsize=20,legendfontsize=20,markersize=8,
        legend_background_color=:transparent,legend_foreground_color=:transparent,
        markerstrokewidth=0,box=false,grid=false,
        left_margin=4mm,bottom_margin=4mm,right_margin=7mm,top_margin=2mm)
    color_TOD = cgrad(:Blues_5,length(sigma_eta_list))[6]
    color_LCO = cgrad(:Reds_5,length(sigma_eta_list))[6]
    colors = Dict("TOD"=>color_TOD,"LCO"=>color_LCO)
    plots = Dict{String,Any}()

    # A/C: mean ± across-trial sample SD, in degrees and deg/s, respectively.
    for (name,output,ylabel) in [
        ("angle",FIG2_PDFS[1],L"\theta_{\mathrm{elbow}}\ (\mathrm{deg})"),
        ("velocity",FIG2_PDFS[3],L"v_{\mathrm{elbow}}\ (\mathrm{deg/s})")]
        p = plot(time_vector,s["$(name)_mean"];ribbon=s["$(name)_SD"],fillalpha=0.18,
            lw=6,color=color_LCO,label=false,xlabel=L"t\ (\mathrm{s})",
            ylabel=ylabel,xlims=(0,0.4),xticks=0:0.1:0.4,
            legend=false)
        name == "angle" && hline!(p,[0.0];ls=:dash,lw=4,color=:black,label=false)
        plots[output] = p
    end

    # B: endpoint sample SD as a percentage of the 20-degree movement amplitude.
    plots[FIG2_PDFS[2]] = plot(sigma_eta_list,s["endpoint_SD_pct_LCO"];lw=6,
        marker=:circle,color=color_LCO,label=false,xlabel=L"\sigma_\eta",
        ylabel=L"\frac{\mathrm{SD}(\theta_{\mathrm{elbow},T})}{20^\circ}\ (\%)",legend=false,
        xlims=(0,0.04),xticks=0:0.02:0.04)

    # D/E/F/G/H: source colors, marker/line conventions, and method comparison.
    definitions = [
        (FIG2_PDFS[4],"analytical_cost",L"C",:log10),
        (FIG2_PDFS[5],"gain_K",L"\log(|K|_p)",:identity),
        (FIG2_PDFS[6],"gain_L",L"\log(|L|_p)",:identity),
        (FIG2_PDFS[7],"variance_transmission_relative",L"\mathcal{V}(\sigma_\eta)/\mathcal{V}(0)",:identity),
        (FIG2_PDFS[8],"mean_peak",L"v_{\mathrm{elbow}}^{\mathrm{peak}}\ (\mathrm{deg/s})",:identity)]
    for (output,quantity,ylabel,scale) in definitions
        p = plot(xlabel=L"\sigma_\eta",ylabel=ylabel,yscale=scale,
            xlims=(0,0.04),xticks=0:0.02:0.04,
            legend=(quantity == "gain_L" ? :topright : :topleft))
        if quantity == "variance_transmission_relative"
            hline!(p,[1.0];ls=:dash,lw=2,color=:gray60,label=false)
        end
        for method in ["TOD","LCO"]
            values = quantity == "analytical_cost" ? results["$(quantity)_$(method)"] : s["$(quantity)_$(method)"]
            if quantity == "mean_peak"
                plot!(p,sigma_eta_list,values;ribbon=s["sem_peak_$(method)"],fillalpha=0.18,
                    lw=6,marker=:circle,color=colors[method],label=latexstring("\\mathrm{",method,"}"))
            else
                plot!(p,sigma_eta_list,values;lw=6,marker=:circle,
                    color=colors[method],label=latexstring("\\mathrm{",method,"}"))
            end
        end
        plots[output] = p
    end
    for output in outputs
        savefig(plots[output],joinpath(plot_dir,output))
    end
    return plots
end
