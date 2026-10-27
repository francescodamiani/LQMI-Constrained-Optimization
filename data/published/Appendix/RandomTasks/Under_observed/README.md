# Original under-observed random-task results

`random_task_severely_underobserved_results.jld2` is copied byte-for-byte from
the rebuttal's `Reviews/results/random_task_severely_underobserved_benchmark/`.
Its original generating script is
`Reviews/scripts/run_random_task_severely_underobserved_benchmark.jl`.
`SHA256SUMS` verifies the complete original file (about 39 KiB).

All top-level fields are retained: the three-row, 31-column result table,
three complete 1,001-value LCO update histories, seeds/dimensions, horizon,
noise amplitudes, rescaling and original iteration/diagnostic settings.
No final TOD/LCO/GD gain arrays were saved historically. GD objectives can be
read at full precision and their comparison statistics verified, but cannot
be reevaluated from absent policies without an expensive rerun. The duplicate
original CSV and old plots are not needed to reproduce the numerical table.

Historical wall-clock timings include compilation and machine effects.
All GD rows hit 10,000 iterations without declaring convergence; LCO histories
at m=8,10 have small positive increments. The source diagnostics are preserved.

```bash
julia scripts/Appendix/RandomTasks/Under_observed/under_observed.jl
```
