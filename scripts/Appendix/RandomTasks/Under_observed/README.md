# Under-observed random-task benchmark

This experiment compares TOD, LCO and direct GD with one scalar
observation as state dimension increases. Run from the repository root using
Julia **1.8.5** and the existing pinned environment:

```bash
# Print Markdown tables from the saved numerical archive; no optimization:
julia scripts/Appendix/RandomTasks/Under_observed/under_observed.jl
# Explicit full rerun, including the expensive numerical GD:
julia scripts/Appendix/RandomTasks/Under_observed/under_observed.jl --simulate --tag=new_random
# Fast reduced numerical check, clearly marked as validation:
julia scripts/Appendix/RandomTasks/Under_observed/under_observed.jl --simulate --validate --tag=check_random
# Inspect a generated run without recomputation:
julia scripts/Appendix/RandomTasks/Under_observed/under_observed.jl --data=data/generated/Appendix/RandomTasks/Under_observed/new_random
```

`under_observed.jl` selects archived tables or explicit simulation.
`run_random_task_severely_underobserved_benchmark.jl` provides the
task builder, optimization body, histories, residuals and status diagnostics.
It reuses the public Figure 1 stable-matrix and joint-gain packing helpers.
The small shared [benchmark_tables.jl](../../benchmark_tables.jl) handles
command-line paths, observability diagnostics and Markdown output. No plotting
packages are loaded and no figures or duplicate table exports are created.

Publication settings are `m=p ∈ {6,8,10}`, `k=1`, `T=10`,
`H=[1,0,…,0]`, spectral radius 0.9 for independent Gaussian A/B draws,
`C=0.2B`, and D equal to 0.2 times the first row of a third rescaled matrix
after taking its absolute value. All external amplitudes are 0.2 and internal
SD is 0.5. Initial means are ones, centered covariances zero, Q/R identity
except `Q₁=0`, `R_T=0`. Independent `MersenneTwister` seeds are **61235,
81235, 101235**, from `1234+10000m+1`.

TOD uses 500 iterations; LCO uses 500 L/K passes from TOD, with `rho=[0]`,
fallbacks `[1e-10,1e-8,1e-6]`, one update each, and both internal GD routes
disabled. The separate comparison is Optim `GradientDescent()` from TOD,
10,000 iterations with its default finite-difference objective and
line search. All costs use the general analytical routines.
`--validate` uses m=2, T=4, three TOD/LCO passes and two GD iterations;
these results are not publication data.

The archive reproduces LCO improvements of **34.9613%, 65.7148%, 60.1891%**.
All GD runs reached the iteration limit without declaring convergence; the
LCO/GD cost gaps are below 0.3%. The m=8,10 LCO histories have small positive
increments and are flagged as non-monotone. These flags remain visible.

The saved random-task archive contains summaries and complete LCO histories,
but **no final gain arrays**. Saved GD objectives cannot therefore
be independently reevaluated from policies without its expensive rerun. New
runs additionally save all six TOD/LCO/GD gain arrays per task. Saved
timings are reported directly. The first LCO timing
is un-warmed; compilation overhead is plausible but not separately measured.

New results go only to `data/generated/Appendix/RandomTasks/Under_observed/TAG/`.
Existing tags are protected; `--force` permits resuming saved rows. Validation
and full runs cannot be mixed. Default commands never write numerical data.
