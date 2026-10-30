#!/usr/bin/env julia
# Increased sensory/motor noise: TOD/LCO costs and fixed-gain ablations; saved data by default.
include(joinpath(@__DIR__,"noise_ablation.jl"))

abspath(PROGRAM_FILE) == (@__FILE__) &&
    noise_ablation_main(;experiments=["high_noise"])
