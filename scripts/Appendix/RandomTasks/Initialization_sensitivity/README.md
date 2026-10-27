# LCO sensitivity to initialization

This experiment tests empirical robustness to initial gains on one random
task; it does not establish global optimality. With Julia **1.8.5**, run from
the repository root:

```bash
# Print tables from all 202 archived runs, without optimization:
julia scripts/Appendix/RandomTasks/Initialization_sensitivity/initialization_sensitivity.jl
# Explicitly repeat all 202 runs with publication parameters:
julia scripts/Appendix/RandomTasks/Initialization_sensitivity/initialization_sensitivity.jl --simulate --tag=new_initializations
# Fast reduced check: one direction, both signs/scales, three complete passes:
julia scripts/Appendix/RandomTasks/Initialization_sensitivity/initialization_sensitivity.jl --simulate --validate --tag=check_initializations
# Inspect a generated run without recomputing:
julia scripts/Appendix/RandomTasks/Initialization_sensitivity/initialization_sensitivity.jl --data=data/generated/Appendix/RandomTasks/Initialization_sensitivity/new_initializations
```

`initialization_sensitivity.jl` selects cached tables or explicit simulation.
`run_lco_initialization_study.jl` provides the random-task builder,
joint-vector construction, one-pass LCO loop and regularization groups.
[benchmark_tables.jl](../../benchmark_tables.jl) recomputes numerical summaries
from every stored row and prints them as Markdown. No plots are generated.

Settings: m=p=k=5, T=10, H=I, external
amplitudes 0.2, internal SD 0.5, spectral radius 0.9, initial means ones and
covariances zero, Q/R identity except Q₁=0 and R_T=0, TOD 500 iterations.
Task RNG is `MersenneTwister(4321)`; direction RNG is
`MersenneTwister(4322)`.

The joint vector is **`theta=vcat(vec(L),vec(K))`**, with length 450 and its
Euclidean norm (equivalently, combined Frobenius norm). Fifty Gaussian vectors
are individually normalized to unit length. For each direction, the order is
negative/unit, negative/large, positive/unit, positive/large. Random starts
are `sign*scale*norm(theta_TOD)*direction`, with scales 1 and 10: these are
gains from the origin, not displacements added to TOD. TOD and zero are the
first two starts. Every direction and final policy is archived.

The simulation calls `Optimal_Control_Estimation_with_Regularized_Inverse` for one
complete L/K pass at a time. This primitive is part of the validated core;
substituting `FullProcedure` per pass would add different acceptance rules.
The whole-run regularization groups are `[0]`, then
`[1e-10,1e-8,1e-6]` only if needed. All saved runs selected zero.
The workflow performs **500 complete passes** using
the `1e-10` relative-change tolerance as a diagnostic; all 202 archived
runs used 500 passes without reaching that tolerance. The one-pass best-checkpoint
return can plateau, so reruns need not reproduce every saved final value
or convergence flag. For example, a zero-start validation rerun gave 6.342619525945239,
versus archived 6.338369494611399 (about 0.0671% higher). The archived result
is retained; identical final values on rerun are not guaranteed.

The archive gives TOD **9.102512125005155**, LCO from TOD
**6.32081865520458**, from zero **6.338369494611399**, and random finals
**6.32013073152–6.33468549686**. Random runs within 0.1% of
TOD-initialized LCO: **183/200 = 91.5%**. Within 0.2%: **199/200 = 99.5%**,
not all; the largest gap is **0.2193836338%**. Every final cost is finite
and below TOD. No solver exceptions occurred. These statistics are calculated from the saved results.

New data go to `data/generated/Appendix/RandomTasks/Initialization_sensitivity/TAG/`.
`--force` permits replacing an existing generated run; this study does not
resume individual starts. Reduced results are clearly labeled and retained
separately. Published files are never overwritten; tables go to stdout only.
