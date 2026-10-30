#!/usr/bin/env julia
# Numerical tables from stored data by default. No plotting dependency.
include(joinpath(@__DIR__,"..","..","..","..","src","project.jl"))
include(joinpath(@__DIR__,"..","..","benchmark_tables.jl"))
include(joinpath(@__DIR__,"noise_statistics.jl"))

function noise_options(args;experiments=["behavioral","original_ablation","high_noise"])
    reoptimize = "--reoptimize" in args
    common_args = String[]
    for arg in args
        if startswith(arg,"--experiment=")
            value = split(arg,"=";limit=2)[2]
            experiments = value == "all" ? ["behavioral","original_ablation","high_noise"] : split(value,",")
        elseif arg != "--reoptimize"
            push!(common_args,arg)
        end
    end
    experiments = [name == "baseline_ablation" ? "original_ablation" : name for name in experiments]
    !isempty(experiments) && all(in(["behavioral","original_ablation","high_noise"]),experiments) || error("Unknown noise experiment")
    length(unique(experiments)) == length(experiments) || error("Repeated noise experiment")
    options = benchmark_options(common_args,"Arm_redundant_task/Noise_ablation")
    reoptimize && !options.simulate && error("--reoptimize requires --simulate")
    return merge(options,(;experiments,reoptimize))
end

function noise_ablation_main(args=ARGS;experiments=["behavioral","original_ablation","high_noise"])
    if "--help" in args
        println("""
        Redundant-arm stabilization: behavioral variability and fixed-policy noise ablations.
          julia scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl
          Default: tables from the included numerical archives; no simulation or plots.
          behavioral_variability.jl       physical-state variability across internal-noise levels
          baseline_noise_ablation.jl      fixed-gain ablations in the baseline noise regime
          high_sensory_motor_noise.jl     TOD/LCO comparison and ablations with increased sensory/motor noise
          --experiment=behavioral,baseline_ablation,high_noise   choose experiments
          Default analyses: $(join(replace.(experiments,"original_ablation"=>"baseline_ablation"),","))
          --simulate      explicitly rerun selected analyses and save numerical results
          --reoptimize    also recompute baseline-regime full-model policies; never ablated policies
          --validate      with --simulate, use a reduced, clearly marked numerical configuration
          --tag=NAME      generated-run directory name
          --data=DIR      print an existing numerical run without recomputation
          --force         permit replacing numerical files inside an existing generated tag
        Baseline analyses reuse archived full-model gains unless --reoptimize is requested.
        High-noise runs recompute TOD and initialize LCO from TOD at each full-model noise level.
        Data: data/published/Appendix/Arm_redundant_task/Noise_ablation/.
        New runs: data/generated/Appendix/Arm_redundant_task/Noise_ablation/TAG/.
        """)
        return
    end
    options = noise_options(args;experiments)
    if options.simulate
        prepare_generated_directory(options.data_dir,"data";force=options.force)
        include(joinpath(PROJECT_ROOT,"src","MyFunctionsOptimalRegularizedSolutions.jl"))
        include(joinpath(@__DIR__,"noise_experiments.jl"))
        results = Base.invokelatest(simulate_noise_experiments,options.data_dir;
            experiments=options.experiments,validation_only=options.validation_only,reoptimize=options.reoptimize)
    else
        results = load_noise_experiments(options.data_dir,options.experiments)
    end
    print_noise_experiments(results)
    return results
end

abspath(PROGRAM_FILE) == (@__FILE__) && noise_ablation_main()
