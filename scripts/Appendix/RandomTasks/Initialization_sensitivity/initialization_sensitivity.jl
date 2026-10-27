#!/usr/bin/env julia
# Numerical tables by default; recomputation is explicitly opt-in. No plots.
include(joinpath(@__DIR__,"..","..","..","..","src","project.jl"))
include(joinpath(@__DIR__,"..","..","benchmark_tables.jl"))

function benchmark_main(args=ARGS)
    if "--help" in args
        println("""
        Task/direction seeds 4321/4322; m=p=k=5,T=10; 202 initializations, at most 500 complete LCO passes.
          julia scripts/Appendix/RandomTasks/Initialization_sensitivity/initialization_sensitivity.jl
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
    options = benchmark_options(args,"RandomTasks/Initialization_sensitivity")
    if options.simulate
        prepare_generated_directory(options.data_dir,"data";force=options.force)
        include(joinpath(PROJECT_ROOT,"src","MyFunctionsOptimalRegularizedSolutions.jl"))
        include(joinpath(@__DIR__,"run_lco_initialization_study.jl"))
        data = Base.invokelatest(simulate_initialization_study,options.data_dir;validation_only=options.validation_only)
    else
        data = load_benchmark(options.data_dir,"lco_initialization_study_results.jld2")
    end
    print_initialization_benchmark(data)
    return data
end

abspath(PROGRAM_FILE) == (@__FILE__) && benchmark_main()
