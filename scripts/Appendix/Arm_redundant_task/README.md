# Supplementary figures: redundant-arm task

This folder contains the redundant-arm stabilization task, its internal-noise
sweep, PCA diagnostics and non-Gaussian noise comparison.

**Full trajectory datasets for the k=3 arm figures are not bundled.** The
figure entry point continues to run full simulations before plotting. It is
separate from the calibrated reaching task in main Figure 2.

The [under-observed arm benchmark](Under_observed/README.md) is a
separate k=1 experiment with complete saved numerical data. Its entry point
prints tables from saved results without simulation or figures.

The separate [noise-ablation analyses](Noise_ablation/README.md) now include
the k=3 full-model policy cache, 5,000-trial component-variance
summaries, fixed-policy noise removals and tenfold sensory/motor-noise results.
Their default command prints numerical tables from archives and creates no figures.

## Run the experiments and plots

Run these commands from the repository root with **Julia 1.8.5**:

```bash
# Run both full experiments and generate their PDF plots:
julia scripts/Appendix/Arm_redundant_task/arm_redundant_task.jl
# Or select one experiment:
julia scripts/Appendix/Arm_redundant_task/arm_redundant_task.jl --experiment=sweep
julia scripts/Appendix/Arm_redundant_task/arm_redundant_task.jl --experiment=non_gaussian
```

These are substantial computations. The sweep uses T=300, 20 internal-noise
levels, 500 TOD iterations, 500 LCO cycles and 500 Monte Carlo trials per
method/noise level. The non-Gaussian experiment uses internal noise 0.1 and
100,000 trials for each distribution.

Only PDF figures are written, under
`figures/generated/Appendix/Arm_redundant_task/`. Repeated runs replace the
corresponding PDFs. Numerical results are kept in memory unless caching is
explicitly requested. Nothing is written to `figures/published/`.

## What each script does

Use `arm_redundant_task.jl` as the command-line entry point. The other files
provide functions loaded by this entry point, rather than separate commands.

| Script | Purpose |
|---|---|
| [arm_redundant_task.jl](arm_redundant_task.jl) | Activates the project, selects experiments, runs simulations by default, optionally saves or loads numerical caches, and dispatches PDF plotting. |
| [arm_task.jl](arm_task.jl) | Defines the arm model, TOD/LCO solve, analytical cost and calls to the core Monte Carlo simulator. Provides standardized Gaussian, Student-t and Beta noise sampling through the simulator's RNG interface. |
| [arm_statistics.jl](arm_statistics.jl) | Supplies gain log-volume, PCA, variance-fraction, global-PC1 readout and peak-velocity statistics, including the source time windows and SEM definitions. |
| [Muscoskeletal_Internal_Noise_Sweep.jl](Muscoskeletal_Internal_Noise_Sweep.jl) | Runs TOD and LCO across the 20 internal-noise levels; retains gains, trajectories, controls, convergence histories, Monte Carlo cost summaries and PCA results for the sweep and PCA plots. |
| [Muscoskeletal_Internal_Noise_Sweep_Plots.jl](Muscoskeletal_Internal_Noise_Sweep_Plots.jl) | Plots cost, gains, PCA, readout strength and peak velocity from the sweep results. |
| [NonGaussian_Redundant_Arm.jl](NonGaussian_Redundant_Arm.jl) | Solves the fixed-noise LCO task once, then evaluates its gains with matched-variance Gaussian, standardized Student-t(5) and standardized Beta(2,5) noise in all five noise channels. Runs Monte Carlo in batches and retains all trial costs. |
| [NonGaussian_Redundant_Arm_Plots.jl](NonGaussian_Redundant_Arm_Plots.jl) | Plots the mean cost ± SEM for the three noise distributions. |

## PDF outputs

All filenames below are relative to
`figures/generated/Appendix/Arm_redundant_task/`.

| Content | PDF output |
|---|---|
| Monte Carlo mean cost ± SEM versus internal noise | `A_cost.pdf` |
| Filter-gain log volume | `B_logK.pdf` |
| Control-gain log volume | `C_logL.pdf` |
| First-PC explained variance | `D_EV1.pdf` |
| LCO PC1 components at high internal noise | `E_PC1_LCO_high_noise.pdf` |
| LCO global-PC1 readout strength | `F_PC1_raw_readout_strength_LCO_globalPC.pdf` |
| Peak joint velocity | `G_peak_velocity.pdf` |
| TOD PCA spectrum and PC weights at zero internal noise | `Appendix_TOD_zero_noise_EV_spectrum.pdf`, `Appendix_TOD_zero_noise_PCk.pdf` |
| TOD PCA spectrum and PC weights at high internal noise | `Appendix_TOD_high_noise_EV_spectrum.pdf`, `Appendix_TOD_high_noise_PCk.pdf` |
| Cost under Gaussian, Student-t and Beta noise | `Noise_Distr_redundant_arm.pdf` |

In the PC-weight filenames, `k` is replaced by the PC number. PCA uses default
99% variance retention and PC bars only for returned explained variance above
1e−4. Consequently, fewer than six PC bars may be
generated. The plots equate zero-noise TOD/LCO control-gain volumes for display only;
raw gains are unchanged.

## Optional numerical caching and plotting without simulations

```bash
# Compute both experiments, save their results locally, and plot:
julia scripts/Appendix/Arm_redundant_task/arm_redundant_task.jl --save-data --tag=arm_run
# Later, regenerate plots from that cache without simulations:
julia scripts/Appendix/Arm_redundant_task/arm_redundant_task.jl --data=data/generated/Appendix/Arm_redundant_task/arm_run
```

Caches are saved under `data/generated/Appendix/Arm_redundant_task/TAG/` as
`Muscoskeletal_Internal_Noise_Sweep.jld2` and
`NonGaussian_Redundant_Arm.jld2`, each containing a `results` object. When caching
only one experiment, use the same `--experiment` selection when replotting.
A missing cache produces an error; it does not trigger simulations.

The complete sweep cache is about 1 GB and must remain local and Git-ignored.
Use a new tag or explicit `--force` to replace an existing cache. `--simulate`
is an alias for the default computation; `--data` cannot be combined with
`--simulate` or `--save-data`. Relative `--data` paths resolve from the repository
root. Run the entry point with `--help` for its options.

See the [data guide](../../../data/README.md) for archive provenance and availability.

Generated caches with fewer noise levels are supported. Plot palettes retain
the full-sweep colors for the 20-level grid; non-Gaussian costs outside
the paper display window use automatic vertical limits.
