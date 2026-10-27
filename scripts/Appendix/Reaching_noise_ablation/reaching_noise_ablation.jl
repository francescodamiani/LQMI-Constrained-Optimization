#!/usr/bin/env julia
# Archived numerical tables by default; no figures, optimization or MC.
include(joinpath(@__DIR__,"..","..","..","src","project.jl"))
include(joinpath(@__DIR__,"..","benchmark_tables.jl"))
include(joinpath(@__DIR__,"ablation_tables.jl"))
using CodecZlib

function reaching_ablation_options(args)
    "--validate" in args && error("The complete analytical ablation is inexpensive; use --simulate without --validate")
    return benchmark_options(args,"Reaching_noise_ablation")
end

function reaching_noise_ablation_main(args=ARGS)
    if "--help" in args
        println("""
        Calibrated Figure 2 reaching: fixed-gain leave-one-noise-source-out tables.
          julia scripts/Appendix/Reaching_noise_ablation/reaching_noise_ablation.jl
          --simulate --tag=NAME   recompute analytical variances with existing Figure 2 LCO gains
          --data=DIR             print a saved generated run
          --force                replace numerical results within an existing generated tag
        Default: saved numerical tables, no computation or output files.
        Published results: data/published/Appendix/Reaching_noise_ablation/.
        New results: data/generated/Appendix/Reaching_noise_ablation/NAME/.
        No gains are reoptimized; no Monte Carlo or figures are generated.
        """)
        return
    end
    options=reaching_ablation_options(args)
    if options.simulate
        prepare_generated_directory(options.data_dir,"data";force=options.force)
        include(joinpath(PROJECT_ROOT,"src","MyFunctionsOptimalRegularizedSolutions.jl"))
        include(joinpath(@__DIR__,"fixed_gain_ablation.jl"))
        Base.invokelatest(recompute_reaching_ablation,options.data_dir)
    end
    data=load_reaching_ablation(options.data_dir)
    print_reaching_ablation(data)
    return data
end

abspath(PROGRAM_FILE)==(@__FILE__) && reaching_noise_ablation_main()
