#!/usr/bin/env julia
# Numerical tables by default; recomputation is explicitly opt-in. No plots.
include(joinpath(@__DIR__,"..","..","..","..","src","project.jl"))
include(joinpath(@__DIR__,"..","..","benchmark_tables.jl"))

function benchmark_main(args=ARGS)
    if "--help" in args
        println("""
        Random tasks: m=p=6,8,10; k=1; T=10; TOD/LCO 500 passes, GD 10000 iterations.
          julia scripts/Appendix/RandomTasks/Under_observed/under_observed.jl
          Default: print Markdown tables from the complete saved JLD2 archive.
          --simulate   rerun the publication settings (can be expensive)
          --validate   with --simulate, run a small explicitly marked numerical test
          --tag=NAME   save new numerical results under this generated-run tag
          --data=DIR   print tables from an existing compatible run, without recomputation
          --force      allow an existing generated run directory (benchmarks resume saved rows)
        No figures or duplicate table files are generated. Published data are never overwritten.
        """)
        return
    end
    options = benchmark_options(args,"RandomTasks/Under_observed")
    if options.simulate
        prepare_generated_directory(options.data_dir,"data";force=options.force)
        include(joinpath(PROJECT_ROOT,"src","MyFunctionsOptimalRegularizedSolutions.jl"))
        include(joinpath(@__DIR__,"run_random_task_severely_underobserved_benchmark.jl"))
        data = Base.invokelatest(simulate_underobserved_random,options.data_dir;validation_only=options.validation_only)
    else
        data = load_benchmark(options.data_dir,"random_task_severely_underobserved_results.jld2")
    end
    print_random_benchmark(data)
    return data
end

abspath(PROGRAM_FILE) == (@__FILE__) && benchmark_main()
