# Under-observed random-task data

`random_task_severely_underobserved_results.jld2` contains a three-row,
31-column result table, three full 1,001-value LCO histories, task/noise
settings, seeds and solver diagnostics.

Final TOD/LCO/GD gains are not included, so objectives cannot be independently
reevaluated from saved policies. GD reaches 10,000 iterations without declaring
convergence; dimension-8/10 LCO histories have small positive increments.
Saved timings include compilation and machine effects.

Run `julia scripts/Appendix/RandomTasks/Under_observed/under_observed.jl` from
the repository root to print tables. The
[experiment guide](../../../../../scripts/Appendix/RandomTasks/Under_observed/README.md)
gives settings and rerun commands. Generated runs are saved under
`data/generated/Appendix/RandomTasks/Under_observed/TAG/` and read with `--data`;
no figures are produced.
