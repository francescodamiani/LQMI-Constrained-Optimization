# Under-observed random-task results

`random_task_severely_underobserved_results.jld2` contains the complete
under-observed random-task benchmark (about 39 KiB). `SHA256SUMS` verifies it.
The simulation script is
`scripts/Appendix/RandomTasks/Under_observed/run_random_task_severely_underobserved_benchmark.jl`.

All top-level fields are retained: the three-row, 31-column result table,
three complete 1,001-value LCO update histories, seeds/dimensions, horizon,
noise amplitudes, rescaling and iteration/diagnostic settings.
Final TOD/LCO/GD gain arrays are not included in the archive. GD objectives can be
read at full precision and their comparison statistics verified, but cannot
be reevaluated from absent policies without an expensive rerun. The numerical table is reproduced directly from the JLD2 archive.

Saved wall-clock timings include compilation and machine effects.
All GD rows hit 10,000 iterations without declaring convergence; LCO histories
at m=8,10 have small positive increments. The source diagnostics are preserved.

```bash
julia scripts/Appendix/RandomTasks/Under_observed/under_observed.jl
```
