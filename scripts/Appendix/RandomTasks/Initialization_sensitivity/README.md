# LCO sensitivity to initialization

Test how LCO solutions depend on initial control and estimation gains for one
random task. The study measures empirical robustness, not global optimality.

## Run

With Julia **1.8.5**, run from the repository root:

```bash
julia scripts/Appendix/RandomTasks/Initialization_sensitivity/initialization_sensitivity.jl
julia scripts/Appendix/RandomTasks/Initialization_sensitivity/initialization_sensitivity.jl --simulate --tag=my_run
julia scripts/Appendix/RandomTasks/Initialization_sensitivity/initialization_sensitivity.jl --data=data/generated/Appendix/RandomTasks/Initialization_sensitivity/my_run
```

The default prints tables from
`data/published/Appendix/RandomTasks/Initialization_sensitivity/` without optimization
or figures. Reruns save to `data/generated/Appendix/RandomTasks/Initialization_sensitivity/TAG/`.
Use `--help` for options, including reduced numerical checks with `--validate`.

## Settings and results

The task has state/control/observation dimensions 5/5/5, T=10, external-noise
amplitudes 0.2 and internal SD 0.5. Task/direction seeds are 4321/4322.
After 500 TOD iterations, LCO is run from TOD, zero and 200 random starts:
50 normalized joint-gain directions, both signs and norm scales 1 and 10.
Random starts are measured from the origin, not added to TOD.

Each run uses 500 complete L/K passes. Regularization groups are `[0]`, then
`[1e-10,1e-8,1e-6]` if needed; all saved runs selected zero. The `1e-10`
relative-change tolerance is reported as a diagnostic, not an early-stop rule.
All 202 costs are finite and below TOD. Of the random starts, 183/200 finish
within 0.1% of TOD-initialized LCO and 199/200 within 0.2%; the largest gap is
0.2193836338%. Reruns can return different best checkpoints, so exact saved
final values and convergence flags are not guaranteed.

`initialization_sensitivity.jl` prints results;
`run_lco_initialization_study.jl` runs the study. See the
[dataset guide](../../../../data/published/Appendix/RandomTasks/Initialization_sensitivity/README.md).
