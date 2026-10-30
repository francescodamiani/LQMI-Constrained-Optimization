#!/usr/bin/env julia
# Physical-state variability across the internal-noise sweep; saved data by default.
include(joinpath(@__DIR__,"noise_ablation.jl"))

abspath(PROGRAM_FILE) == (@__FILE__) &&
    noise_ablation_main(;experiments=["behavioral"])
