# Redundant-arm task

Study stabilization of a three-joint arm: internal-noise effects on costs,
gains, PCA structure and peak velocity, plus a non-Gaussian noise comparison.

## Run

With Julia **1.8.5**, run from the repository root:

```bash
julia scripts/Appendix/Arm_redundant_task/arm_redundant_task.jl
julia scripts/Appendix/Arm_redundant_task/arm_redundant_task.jl --experiment=sweep
julia scripts/Appendix/Arm_redundant_task/arm_redundant_task.jl --experiment=non_gaussian
```

**The default runs full simulations:** figure trajectory datasets are not
bundled. All outputs are PDFs in `figures/generated/Appendix/Arm_redundant_task/`.
The sweep produces `A_cost.pdf`–`G_peak_velocity.pdf` and the
`Appendix_TOD_*` PCA plots; the distribution comparison produces
`Noise_Distr_redundant_arm.pdf`.

## Settings

The model has 6 states, 9 controls and 3 angle observations, T=300, Δt=0.01 s,
zero initial means/covariances, state weights `(1,1,1,.001,.001,.001)` and
control weight 0.01. External-noise amplitudes are 0.1. The sweep uses
20 internal-noise levels in `[0,0.4]`, 500 TOD iterations, 500 LCO cycles and
500 Monte Carlo trials per method/level. LCO starts from TOD; regularization
is `[0]` with fallback `[1e-10,1e-8,1e-6]`.

The non-Gaussian comparison fixes internal noise at 0.1 and uses 100,000 trials
each for matched-variance Gaussian, Student-t(5) and Beta(2,5) noise. Costs
are mean ± SEM. PCA uses default 99% variance retention and plots returned
PCs with explained variance above `1e-4`; fewer than six PCs may appear.
The zero-noise control-gain volumes are equated for display only.

## Save and replot simulations

```bash
julia scripts/Appendix/Arm_redundant_task/arm_redundant_task.jl --save-data --tag=my_run
julia scripts/Appendix/Arm_redundant_task/arm_redundant_task.jl --data=data/generated/Appendix/Arm_redundant_task/my_run
```

`--save-data` stores results in `data/generated/Appendix/Arm_redundant_task/TAG/`;
a full sweep cache is about 1 GB. Use the same `--experiment` when replotting
a partial cache. `--data` never triggers simulations; existing cache names
require `--force`.

`arm_redundant_task.jl` is the entry point. `arm_task.jl` and `arm_statistics.jl`
define the model and statistics; `Muscoskeletal_Internal_Noise_Sweep.jl` and
`NonGaussian_Redundant_Arm.jl` run experiments, with matching `_Plots.jl` helpers.
The [behavioral variability and noise-source analyses](Noise_ablation/README.md)
follow the paper's discussion in three parts: physical-state variability,
baseline fixed-gain ablations, and the higher sensory/motor-noise regime.
These and the [under-observed arm](Under_observed/README.md) analysis have saved
numerical data and print tables by default. See `--help` for options.
