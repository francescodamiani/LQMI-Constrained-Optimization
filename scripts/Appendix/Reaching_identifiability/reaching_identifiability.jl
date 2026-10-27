#!/usr/bin/env julia
# Saved-data tables by default; no plotting or simulation dependencies loaded.
include(joinpath(@__DIR__,"..","..","..","src","project.jl"))
include(joinpath(@__DIR__,"..","benchmark_tables.jl"))
include(joinpath(@__DIR__,"identifiability_statistics.jl"))
using CodecZlib # Load before reading compressed JLD2 data on Julia 1.8.5.

const IDENTIFIABILITY_FILES = Dict("overestimation"=>"sensory_overestimation.jld2",
    "noise_sources"=>"matched_noise_sources.jld2","underestimation"=>"sensory_underestimation.jld2")

function identifiability_options(args)
    experiments=["overestimation","noise_sources","underestimation"]
    common=String[]
    for arg in args
        if startswith(arg,"--experiment=")
            value=split(arg,"=";limit=2)[2]
            experiments=value=="all" ? experiments : String.(split(value,","))
        elseif arg!="--full"
            push!(common,arg)
        end
    end
    !isempty(experiments) && all(haskey(IDENTIFIABILITY_FILES,e) for e in experiments) || error("Unknown experiment")
    length(unique(experiments))==length(experiments) || error("Repeated experiment")
    return merge(benchmark_options(common,"Reaching_identifiability"),(;experiments,full="--full" in args))
end

function reaching_identifiability_main(args=ARGS)
    if "--help" in args
        println("""
        Internal-noise identifiability comparisons: numerical tables, no figures.
          julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl
          --full                    every A/B condition (C always includes both complete sweeps)
          --experiment=overestimation,noise_sources,underestimation   select analyses
          --simulate --tag=NAME     explicitly rerun selected experiments with publication settings
          --validate                with --simulate: 8 samples, 64 trials, 3 iterations, reduced grids
          --data=DIR                tables from an existing saved run
          --force                   replace results within an existing generated tag
        Published data are never overwritten. Generated data: data/generated/Appendix/Reaching_identifiability/NAME/.
        A/B reuse the saved internal-noise reference, including its zero-noise policy.
        """)
        return
    end
    options=identifiability_options(args)
    if options.simulate
        prepare_generated_directory(options.data_dir,"data";force=options.force)
        include(joinpath(@__DIR__,"identifiability_simulation.jl"))
        published=joinpath(PROJECT_ROOT,"data","published","Appendix","Reaching_identifiability")
        Base.invokelatest(IdentifiabilitySimulation.simulate_identifiability,
            options.data_dir,published,options.experiments;validation_only=options.validation_only)
    end
    data=Dict(e=>load_benchmark(options.data_dir,IDENTIFIABILITY_FILES[e]) for e in options.experiments)
    if "noise_sources" in options.experiments
        data["reference"]=load_benchmark(joinpath(options.data_dir,"reference"),"sweep.jld2")
    end
    print_identifiability(data;full=options.full)
    return data
end

abspath(PROGRAM_FILE)==(@__FILE__) && reaching_identifiability_main()
