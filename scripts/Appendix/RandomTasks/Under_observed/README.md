# Under-observed random-task benchmark

Compare TOD, LCO and direct GD as state dimension increases with only one
scalar observation.

## Run

With Julia **1.8.5**, run from the repository root:

```bash
julia scripts/Appendix/RandomTasks/Under_observed/under_observed.jl
julia scripts/Appendix/RandomTasks/Under_observed/under_observed.jl --simulate --tag=my_run
julia scripts/Appendix/RandomTasks/Under_observed/under_observed.jl --data=data/generated/Appendix/RandomTasks/Under_observed/my_run
```

The default prints tables from `data/published/Appendix/RandomTasks/Under_observed/`.
It runs no simulation and produces no figures. Explicit reruns save to
`data/generated/Appendix/RandomTasks/Under_observed/TAG/`; GD can be expensive.
Use `--help` for options, including reduced numerical checks with `--validate`.

## Settings and limitations

State/control dimensions are 6, 8 and 10; observation dimension is 1, T=10,
`H=[1,0,…,0]`, external-noise amplitudes 0.2 and internal SD 0.5.
Task seeds are 61235, 81235 and 101235. TOD uses 500 iterations; LCO uses
500 passes from TOD with regularization `[0]`, fallback `[1e-10,1e-8,1e-6]`.
Direct GD starts from TOD and runs for up to 10,000 iterations.

Saved LCO improvements over TOD are 34.9613%, 65.7148% and 60.1891%.
All GD runs reached their iteration limit; dimension-8/10 LCO histories have
small positive increments. The dataset contains costs, timings and full LCO
histories but no final gains, so policies cannot be reevaluated without a rerun.
Timings include machine and compilation effects.

`under_observed.jl` is the entry point;
`run_random_task_severely_underobserved_benchmark.jl` runs the experiment.
See the [dataset guide](../../../../data/published/Appendix/RandomTasks/Under_observed/README.md).
