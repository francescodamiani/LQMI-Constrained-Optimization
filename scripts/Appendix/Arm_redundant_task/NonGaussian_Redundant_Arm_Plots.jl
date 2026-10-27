# Original appendix Figure 6; preserve the source's mean ± sample SEM and axes.
using Plots, LaTeXStrings

function plot_nongaussian_arm(results,plot_dir)
    default(dpi=500,thickness_scaling=2.2,size=(1300,600),
        legend_background_color=:transparent,legend_foreground_color=:transparent,
        markerstrokewidth=0,box=false,grid=false)
    noise_names = [L"\mathrm{Gaussian}",L"\mathrm{Student\!-t}\ (\nu=5)",L"\beta(2,5)"]
    labels = ["Gaussian","Student-t (nu=5)","Beta(2,5)"]
    mean_costs = [results.mean_cost_dict[label] for label in labels]
    sem_costs = [results.sem_cost_dict[label] for label in labels]
    # Preserve the paper window when it contains the saved error bars.
    # Other generated runs use automatic bounds rather than clipping data.
    limits = all(mean_costs .- sem_costs .>= 358.5) && all(mean_costs .+ sem_costs .<= 360.5) ?
        (358.5,360.5) : :auto
    plt_mean = scatter(1:3,mean_costs;yerror=sem_costs,xticks=(1:3,noise_names),
        markersize=5,markercolor=:black,markerstrokecolor=:black,linewidth=2,color=:black,
        xlabel=L"\mathrm{Noise\ distribution}",ylabel=L"C",legend=false,grid=false,
        xtickfontsize=12,ytickfontsize=12,xguidefontsize=18,yguidefontsize=18,
        ylims=limits)
    savefig(plt_mean,joinpath(plot_dir,"Noise_Distr_redundant_arm.pdf"))
    return nothing
end
