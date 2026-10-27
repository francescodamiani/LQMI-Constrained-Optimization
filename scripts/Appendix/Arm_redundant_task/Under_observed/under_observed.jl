#!/usr/bin/env julia
# Numerical tables by default; recomputation is explicitly opt-in. No plots.
include(joinpath(@__DIR__,"..","..","..","..","src","project.jl"))
include(joinpath(@__DIR__,"..","..","benchmark_tables.jl"))

function benchmark_main(args=ARGS)
    if "--help" in args
        println("""
        Historical stabilization arm: m=6,p=9,k=1,T=300; six noise levels; TOD/LCO 500 passes. No GD.
          julia scripts/Appendix/Arm_redundant_task/Under_observed/under_observed.jl
          Default: print Markdown tables from the complete original JLD2 archive.
          --simulate   rerun the original publication settings (can be expensive)
          --validate   with --simulate, run a small explicitly marked numerical test
          --tag=NAME   save new numerical results under this generated-run tag
          --data=DIR   print tables from an existing compatible run, without recomputation
          --force      allow an existing generated run directory (benchmarks resume saved rows)
        No figures or duplicate table files are generated. Historical data are never overwritten.
        """)
        return
    end
    options = benchmark_options(args,"Arm_redundant_task/Under_observed")
    if options.simulate
        prepare_generated_directory(options.data_dir,"data";force=options.force)
        include(joinpath(PROJECT_ROOT,"src","MyFunctionsOptimalRegularizedSolutions.jl"))
        include(joinpath(@__DIR__,"run_arm_redundant_task_severely_underobserved_benchmark.jl"))
        data = Base.invokelatest(simulate_underobserved_arm,options.data_dir;validation_only=options.validation_only)
    else
        data = load_benchmark(options.data_dir,"arm_redundant_task_severely_underobserved_results.jld2")
    end
    print_arm_benchmark(data)
    return data
end

abspath(PROGRAM_FILE) == (@__FILE__) && benchmark_main()
