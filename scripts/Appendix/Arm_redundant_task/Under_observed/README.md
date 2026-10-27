# Under-observed redundant-arm stabilization

Compare TOD and LCO when only elbow angle is observed in the three-joint
stabilization task.

## Run

With Julia **1.8.5**, run from the repository root:

```bash
julia scripts/Appendix/Arm_redundant_task/Under_observed/under_observed.jl
julia scripts/Appendix/Arm_redundant_task/Under_observed/under_observed.jl --simulate --tag=my_run
julia scripts/Appendix/Arm_redundant_task/Under_observed/under_observed.jl --data=data/generated/Appendix/Arm_redundant_task/Under_observed/my_run
```

The default prints saved tables from
`data/published/Appendix/Arm_redundant_task/Under_observed/`, with no simulation
or figures. Reruns save to `data/generated/Appendix/Arm_redundant_task/Under_observed/TAG/`.
Use `--help` for options, including reduced numerical checks with `--validate`.

## Settings and results

The [arm model](../README.md) has 6 states, 9 controls, T=300 and Δt=0.01 s.
Here `H=[0,1,0,0,0,0]`, with observability rank 2/6. All external-noise amplitudes
are 0.1; internal SDs are `[.05,.10,.15,.20,.30,.40]`. TOD and LCO use
500 iterations/passes, with LCO initialized from TOD and regularization `[0]`,
fallback `[1e-10,1e-8,1e-6]`. There is no GD or Monte Carlo step.

LCO cost improvements range from 15.0056% to 91.7578%. All twelve policies and
six full LCO histories are included, allowing analytical objectives to be
reevaluated. Saved timings are machine-dependent.

`under_observed.jl` is the entry point;
`run_arm_redundant_task_severely_underobserved_benchmark.jl` runs the experiment.
See the [dataset guide](../../../../data/published/Appendix/Arm_redundant_task/Under_observed/README.md).
