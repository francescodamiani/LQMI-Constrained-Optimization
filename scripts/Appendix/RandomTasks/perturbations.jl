#!/usr/bin/env julia
# LCO local perturbations: saved data by default; simulation is explicitly opt-in.
include(joinpath(@__DIR__, "..", "..", "..", "src", "project.jl"))
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
include(joinpath(@__DIR__, "perturbation_helpers.jl"))
include(joinpath(@__DIR__, "perturbation_plots.jl"))

function perturbation_options(args)
    simulate, force = "--simulate" in args, "--force" in args
    data, tag = "", ""
    for arg in args
        if startswith(arg, "--data=")
            data = split(arg, "="; limit=2)[2]
            isempty(data) && error("--data requires a folder")
        elseif startswith(arg, "--tag=")
            tag = split(arg, "="; limit=2)[2]
            occursin(r"^[A-Za-z0-9][A-Za-z0-9_-]*$", tag) || error("Invalid run tag")
        elseif !(arg in ["--simulate", "--force", "--help"])
            error("Unknown argument: $arg")
        end
    end
    simulate && !isempty(data) && error("--data and --simulate cannot be combined")
    isempty(tag) && (tag = "run_" * Dates.format(now(), "yyyymmdd_HHMMSSsss"))
    data_dir = isempty(data) ? joinpath(PROJECT_ROOT, "data", "published", "Appendix", "RandomTasks") :
        (isabspath(data) ? data : joinpath(PROJECT_ROOT, data))
    simulate && (data_dir = joinpath(PROJECT_ROOT, "data", "generated", "Appendix", "RandomTasks", tag))
    figure_dir = joinpath(PROJECT_ROOT, "figures", "generated", "Appendix", "RandomTasks")
    return (;simulate, force, data_dir=normpath(data_dir), figure_dir)
end

function perturbation_main(args=ARGS)
    if "--help" in args
        println("""
        Local perturbations of LCO (Julia 1.8.5)
          julia scripts/Appendix/RandomTasks/perturbations.jl
            Generate all three PDF panels from the included numerical dataset.
          --simulate   explicitly recompute the full scan (2.4 million costs)
          --tag=NAME   name a newly simulated numerical run
          --force      replace generated numerical data for that tag
          --data=DIR   plot a compatible saved run without simulations
        PDFs: figures/generated/Appendix/RandomTasks/.
        New data: data/generated/Appendix/RandomTasks/TAG/.
        The default command runs no simulations. figures/published/ is never written.
        """)
        return
    end
    options = perturbation_options(args)
    if options.simulate
        prepare_generated_directory(options.data_dir, "data"; force=options.force)
        include(joinpath(@__DIR__, "..", "..", "..", "src", "MyFunctionsOptimalRegularizedSolutions.jl"))
        include(joinpath(@__DIR__, "Test_LCO_perurbations_for_local_minimum.jl"))
        results = Base.invokelatest(simulate_perturbations, options.data_dir)
    else
        results = load_perturbations(options.data_dir)
    end
    prepare_generated_directory(options.figure_dir, "figures"; force=true)
    plot_perturbations(results, options.figure_dir)
    println("Saved LCO perturbation PDFs to ", options.figure_dir)
end

abspath(PROGRAM_FILE) == (@__FILE__) && perturbation_main()
