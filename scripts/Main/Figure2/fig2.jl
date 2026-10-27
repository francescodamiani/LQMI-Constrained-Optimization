#!/usr/bin/env julia
# Main Figure 2: calibrated reaching data, eight PDFs, no default simulations.
include(joinpath(@__DIR__,"..","..","..","src","project.jl"))
ENV["GKSwstype"] = get(ENV,"GKSwstype","100")
include(joinpath(@__DIR__,"fig2_data.jl"))
include(joinpath(@__DIR__,"fig2_statistics.jl"))
include(joinpath(@__DIR__,"fig2_plots.jl"))

function fig2_options(args)
    simulate,force = "--simulate" in args,"--force" in args
    data,tag = "",""
    panels = collect("ABCDEFGH")
    for arg in args
        if startswith(arg,"--data=")
            data = split(arg,"=";limit=2)[2]
            isempty(data) && error("--data requires a folder")
        elseif startswith(arg,"--tag=")
            tag = split(arg,"=";limit=2)[2]
            occursin(r"^[A-Za-z0-9][A-Za-z0-9_-]*$",tag) || error("Invalid run tag")
        elseif startswith(arg,"--panels=")
            panels = collect(replace(split(arg,"=";limit=2)[2],","=>""))
        elseif !(arg in ["--simulate","--force","--help"])
            error("Unknown argument: $arg")
        end
    end
    simulate && !isempty(data) && error("--data and --simulate cannot be combined")
    !isempty(panels) && all(p -> p in "ABCDEFGH",panels) || error("Select Figure 2 panels A–H")
    length(unique(panels)) == length(panels) || error("Repeated panel")
    panels = [p for p in "ABCDEFGH" if p in panels]
    isempty(tag) && (tag = "run_" * Dates.format(now(),"yyyymmdd_HHMMSSsss"))
    data_dir = isempty(data) ? joinpath(PROJECT_ROOT,"data","published","Main","Figure2") :
        (isabspath(data) ? data : joinpath(PROJECT_ROOT,data))
    simulate && (data_dir = joinpath(PROJECT_ROOT,"data","generated","Main","Figure2",tag))
    figure_dir = joinpath(PROJECT_ROOT,"figures","generated","Main","Figure2")
    return (;simulate,force,data_dir=normpath(data_dir),figure_dir,panels)
end

function fig2_main(args=ARGS)
    if "--help" in args
        println("""
        Figure 2: calibrated elbow-reaching task (Julia 1.8.5)
          julia scripts/Main/Figure2/fig2.jl   eight PDFs from saved data
          --simulate    recompute 20 noise levels, TOD/LCO 500 iterations, 2000 MC trials
          --tag=NAME    name the new numerical run
          --force       replace generated numerical results for this tag
          --data=DIR    plot a compatible saved run without simulations
          --panels=G    regenerate only panel G (or select letters A–H)
        Data: data/published/Main/Figure2/ (complete saved calibrated sweep).
        New runs: data/generated/Main/Figure2/TAG/.
        PDFs: figures/generated/Main/Figure2/. No composite or smoke figures.
        Representative elbow trajectories: sigma_eta=0.019, mean ± SD, no legends.
        Peak velocity: trial mean of peak absolute elbow velocity, with SEM.
        Variance transmission: sum Tr[L Cov(z) L'] / sum Tr[Cov(z)], times 1:T-2.
        Panel G: V(sigma_eta)/V(0), separately for each method, baseline=1.
        figures/published/ is never written.
        """)
        return
    end
    options = fig2_options(args)
    if options.simulate
        prepare_generated_directory(options.data_dir,"data";force=options.force)
        include(joinpath(@__DIR__,"..","..","..","src","MyFunctionsOptimalRegularizedSolutions.jl"))
        include(joinpath(@__DIR__,"reaching_task.jl"))
        include(joinpath(@__DIR__,"fig2_simulation.jl"))
        results = Base.invokelatest(simulate_fig2,options.data_dir)
    else
        results = load_fig2(options.data_dir)
    end
    prepare_generated_directory(options.figure_dir,"figures";force=true)
    plot_fig2(results,options.figure_dir;panels=options.panels)
    println("Saved Figure 2 panels ",join(options.panels,",")," to ",options.figure_dir)
end

abspath(PROGRAM_FILE) == (@__FILE__) && fig2_main()
