#!/usr/bin/env julia
# Fixed-gain noise ablations in the baseline stabilization model; saved data by default.
include(joinpath(@__DIR__,"noise_ablation.jl"))

abspath(PROGRAM_FILE) == (@__FILE__) &&
    noise_ablation_main(;experiments=["original_ablation"])
