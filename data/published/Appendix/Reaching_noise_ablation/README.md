# Calibrated reaching noise-ablation data

`calibrated_fixed_gain_noise_ablation_subset.csv` is a **byte-for-byte copy**
of the complete original ablation-specific numerical output, **4,918 bytes**,
35 rows × 8 columns. Despite “subset” in its historical filename, it contains
all seven original conditions at all five selected internal-noise levels.
Every available numerical value is retained at its original full precision.
No JLD conversion or extra dependency is needed. Default tables display 12
significant digits; CSV parsing restores the original Float64 values.

Original source paths below are relative to `reference_private/Reviews/`
and document provenance only; public code never accesses them:

| Historical file | Public use |
|---|---|
| `scripts/analyze_elbow_reaching_calibrated_fixed_gain_noise_ablation.jl` | Migrated under `scripts/Appendix/Reaching_noise_ablation/`; original loops, noise removals, tolerance and SD metrics retained using existing public analysis helpers. |
| `results/elbow_reaching_internal_noise_sweep_baseline_configuration_calibrated/calibrated_fixed_gain_noise_ablation_subset.csv` | Complete original numerical CSV included here without modification. |
| Same directory, `elbow_reaching_internal_noise_sweep_baseline_configuration_calibrated.jld` | **Already migrated** as the complete public Figure 2 archive. Reuse `data/published/Main/Figure2/sweep.jld2` directly for model/gains; no duplicate large archive. |
| Same directory, `calibrated_sigma_xi_0_01_leave_one_noise_out_results_latex.txt` and `short_terminal_elbow_noise_ablation_summary.txt` | Read to verify the reported reductions, combined sensory definition, baseline exclusions and nonadditive interpretation. Duplicate text exports are not included. |

SHA-256 provenance:

```text
0829a0cd68ab1c8c6f9c722e16f3c19c8e901aa8d37eb23bf9d332d63049ea9e  original and public ablation CSV
1ed6c2a5784404544484201f92567fbe2485de608445956871a8c0d53b895e23  original ablation script
5a02bcfcb9c2f9f5ac3440b0304713c27da68719a85408051ee704b3b9d83a0d  original full calibrated JLD
55527e57d4b9e2faace3cf8ad2413892229da80daab7ac473dd9549045f8118f  reused public Figure 2 sweep.jld2
```

[`SHA256SUMS`](SHA256SUMS) verifies the sole new published numerical file.
The historical full archive has 75 fields; all 69 non-trajectory fields of
its public metadata/gain file were checked with `isequal` against the
original. The existing Figure 2 migration preserves all six trajectory fields
separately. Nothing from that approximately 500 MB archive is duplicated here.

| CSV column | Definition |
|---|---|
| `sigma_eta` | Actual selected amplitudes 0,.011,.019,.029,.040, from nominal nearest-grid targets 0,.01,.02,.03,.04 |
| `condition` | Full model; five individual source removals; genuine combined additive+multipl. sensory removal |
| `terminal_elbow_SD_pct` | Analytical centered terminal elbow SD, as percent of the 20° movement |
| `mean_angle_SD` | Mean instantaneous SD over angles 1:3 and samples 2:41; radians |
| `mean_velocity_SD` | Mean instantaneous SD over angular velocities 4:6 and samples 2:41; rad/s |
| `terminal_elbow_SD_reduction_pct` | 100*(1-terminal SD after removal / full terminal SD) |
| `mean_angle_SD_reduction_pct` | Same relative reduction for the mean-angle SD metric |
| `mean_velocity_SD_reduction_pct` | Same relative reduction for the mean angular-velocity SD metric |

All evaluations freeze the full-model optimized LCO gains. The combined
sensory removal zeros Ωω and D simultaneously; it is not a sum of individual
percentages. Mean SD is not sqrt(mean variance). The original physical-variance
negativity threshold is −1e-10 before clipping roundoff to zero.

The rounded reviewer ranges refer to **nonzero** internal noise. Analytical
full-model endpoint SD is 8.63464–10.54364%; the approximately 8.55–10.46%
Figure 2 values are archived Monte Carlo sample SDs. These distinctions are
preserved rather than changing the original analytical output.

```bash
julia scripts/Appendix/Reaching_noise_ablation/reaching_noise_ablation.jl
julia scripts/Appendix/Reaching_noise_ablation/reaching_noise_ablation.jl --simulate --tag=new_run
```

Default tables read only this small CSV and write no files. Explicit reruns
load existing Figure 2 metadata/gains and save **one**
`calibrated_reaching_ablation.jld2` under
`data/generated/Appendix/Reaching_noise_ablation/TAG/`. It stores the complete
35-row table, all physical variances (6×41×7×5), selected/nominal levels,
indices, condition names, window, normalization, tolerance, fixed-policy/no-
optimization/no-MC flags and the source Figure 2 metadata hash. It saves no
second copy of gains, model or trajectory data. No plots or Monte Carlo runs
are produced; no existing experiment is rerun during ablation.
