#!/usr/bin/env julia
# Redundant-arm stabilization experiments: run simulations and write PDFs.
# Full trajectory datasets for these arm figures are not bundled.
include(joinpath(@__DIR__,"..","..","..","src","project.jl"))
ENV["GKSwstype"] = get(ENV,"GKSwstype","100")
using JLD2
include(joinpath(@__DIR__,"arm_statistics.jl"))
include(joinpath(@__DIR__,"Muscoskeletal_Internal_Noise_Sweep_Plots.jl"))
include(joinpath(@__DIR__,"NonGaussian_Redundant_Arm_Plots.jl"))

function arm_options(args)
    experiment, data, tag = "all", "", ""
    save_data, force = "--save-data" in args, "--force" in args
    for arg in args
        if startswith(arg,"--experiment=")
            experiment = split(arg,"=";limit=2)[2]
        elseif startswith(arg,"--data=")
            data = split(arg,"=";limit=2)[2]
        elseif startswith(arg,"--tag=")
            tag = split(arg,"=";limit=2)[2]
        elseif !(arg in ["--simulate","--save-data","--force","--help"])
            error("Unknown argument: $arg")
        end
    end
    experiment in ["all","sweep","non_gaussian"] || error("Unknown arm experiment: $experiment")
    !isempty(data) && ("--simulate" in args || save_data) && error("--data cannot be combined with --simulate or --save-data")
    isempty(tag) || occursin(r"^[A-Za-z0-9][A-Za-z0-9_-]*$",tag) || error("Invalid run tag")
    isempty(tag) && (tag = "run_"*Dates.format(now(),"yyyymmdd_HHMMSSsss"))
    data_dir = isempty(data) ? joinpath(PROJECT_ROOT,"data","generated","Appendix","Arm_redundant_task",tag) :
        (isabspath(data) ? data : joinpath(PROJECT_ROOT,data))
    figure_dir = joinpath(PROJECT_ROOT,"figures","generated","Appendix","Arm_redundant_task")
    return (;experiment,data_dir=normpath(data_dir),figure_dir,save_data,force,simulate=isempty(data))
end

function arm_main(args=ARGS)
    if "--help" in args
        println("""
        Redundant-arm stabilization task (Julia 1.8.5)
          julia scripts/Appendix/Arm_redundant_task/arm_redundant_task.jl
            Recompute the noise sweep and non-Gaussian test, then plot PDFs.
            This performs full simulations; figure trajectory datasets are not bundled.
          --simulate                 explicit alias for the default computation
          --experiment=sweep         internal-noise sweep and PCA diagnostics
          --experiment=non_gaussian  non-Gaussian noise comparison
          --save-data --tag=NAME      optionally cache numerical results locally
          --data=FOLDER              plot a compatible cache without simulations
          --force                    replace an explicitly named numerical cache
        PDFs: figures/generated/Appendix/Arm_redundant_task/.
        Optional cache: data/generated/Appendix/Arm_redundant_task/TAG/.
        A complete sweep cache is about 1 GB; keep it local and Git-ignored.
        figures/published/ is reserved for final manuscript figures.
        """)
        return
    end
    options = arm_options(args)
    experiments = options.experiment == "all" ? ["sweep","non_gaussian"] : [options.experiment]
    filenames = Dict("sweep"=>"Muscoskeletal_Internal_Noise_Sweep.jld2",
        "non_gaussian"=>"NonGaussian_Redundant_Arm.jld2")
    if !options.simulate
        for experiment in experiments
            isfile(joinpath(options.data_dir,filenames[experiment])) ||
                error("Missing arm cache: $(joinpath(options.data_dir,filenames[experiment])). Run the experiment to generate plots; figure trajectory datasets are not bundled.")
        end
    end
    prepare_generated_directory(options.figure_dir,"figures";force=true)
    if options.simulate
        options.save_data && prepare_generated_directory(options.data_dir,"data";force=options.force)
        println("Running redundant-arm simulations with full publication settings.")
        include(joinpath(@__DIR__,"..","..","..","src","MyFunctionsOptimalRegularizedSolutions.jl"))
        isdefined(@__MODULE__,:ArmNoiseRNG) || include(joinpath(@__DIR__,"arm_task.jl"))
        include(joinpath(@__DIR__,"Muscoskeletal_Internal_Noise_Sweep.jl"))
        include(joinpath(@__DIR__,"NonGaussian_Redundant_Arm.jl"))
    end
    for experiment in experiments
        if options.simulate
            results = experiment == "sweep" ? Base.invokelatest(simulate_internal_noise_sweep) :
                Base.invokelatest(simulate_nongaussian_arm)
            options.save_data && JLD2.jldsave(joinpath(options.data_dir,filenames[experiment]);results)
        else
            results = JLD2.load(joinpath(options.data_dir,filenames[experiment]),"results")
        end
        if experiment == "sweep"
            plot_internal_noise_sweep(results,options.figure_dir)
        else
            plot_nongaussian_arm(results,options.figure_dir)
        end
        results = nothing
        GC.gc()
    end
    println("Saved arm PDFs to $(options.figure_dir)")
end

abspath(PROGRAM_FILE) == (@__FILE__) && arm_main()
