#!/usr/bin/env julia
# Figure 1b–g: saved data by default; simulation is explicitly opt-in.
include(joinpath(@__DIR__, "..", "..", "..", "src", "project.jl"))
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100") # file rendering on machines without a display
using DataFrames, JLD2, Plots, LaTeXStrings, Statistics, LinearAlgebra, Random
include(joinpath(@__DIR__, "fig1_plots.jl"))

function main(args=ARGS)
    if "--help" in args
        println("""
        Figure 1b–g (Julia 1.8.5)
          julia scripts/Main/Figure1/fig1.jl                       all six panels from saved data
          julia scripts/Main/Figure1/fig1.jl --simulate            full publication settings (GD: days)
          julia scripts/Main/Figure1/fig1.jl --simulate --smoke    reduced numerical test; no figures
          --tag=NAME    isolate newly simulated numerical data
          --force       replace existing generated numerical data for this tag
          --data=FOLDER load a compatible dataset folder (relative to repository root)
          --panels=b,c,d,e,f,g select panels and the simulations needed for them
        PDFs: figures/generated/Main/Figure1/fig1b.pdf through fig1g.pdf, plus fig1_b-g.pdf.
        Plots replace these derived PDFs; figures/published is reserved for final manuscript figures.
        New numerical data: data/generated/Main/Figure1/TAG.
        """)
        return
    end
    options = fig1_options(args)
    if options.simulate
        prepare_generated_directory(options.data_dir, "data"; force=options.force)
        options.smoke || prepare_generated_directory(options.figure_dir, "figures"; force=true)
        options.smoke && println("SMOKE TEST: reduced settings; outputs are NOT publication results.")
        include(joinpath(@__DIR__, "..", "..", "..", "src", "MyFunctionsOptimalRegularizedSolutions.jl"))
        include(joinpath(@__DIR__, "fig1_tasks.jl"))
        include(joinpath(@__DIR__, "fig1_dimensions.jl"))
        include(joinpath(@__DIR__, "fig1_convergence.jl"))
        include(joinpath(@__DIR__, "fig1_noise.jl"))
        # Invokelatest allows these optional helpers to be included only in simulation mode.
        if any(p -> p in "bc", options.panels)
            Base.invokelatest(simulate_dimensions, options.data_dir; smoke=options.smoke)
        end
        if 'd' in options.panels
            Base.invokelatest(simulate_convergence, options.data_dir; smoke=options.smoke)
        end
        if any(p -> p in "efg", options.panels)
            Base.invokelatest(simulate_internal_noise, options.data_dir; smoke=options.smoke)
        end
        if options.smoke
            println("Saved smoke numerical data to $(options.data_dir). No figures generated.")
            return
        end
    end
    datasets = load_fig1(options.data_dir, options.panels)
    # Derived PDFs have stable names; saved numerical datasets remain protected.
    # The plotting function rejects reduced smoke data before writing any files.
    prepare_generated_directory(options.figure_dir, "figures"; force=true)
    plot_fig1(datasets, options.figure_dir; panels=options.panels)
    println("Saved panels $(join(options.panels, ',')) to $(options.figure_dir)")
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
