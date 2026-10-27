# Calibrated reaching noise-ablation data

`calibrated_fixed_gain_noise_ablation_subset.csv` contains the complete
ablation results: **4,918 bytes**, 35 rows × 8 columns, covering all seven
conditions at all five selected internal-noise levels. All numerical values
retain full precision.
No JLD conversion or extra dependency is needed. Default tables display 12
significant digits; CSV parsing restores the stored Float64 values.

The plotting-independent analysis is in
`scripts/Appendix/Reaching_noise_ablation/`. Reruns reuse
`data/published/Main/Figure2/sweep.jld2` for the model and gains; no trajectory
archive is duplicated.

SHA-256 provenance:

```text
0829a0cd68ab1c8c6f9c722e16f3c19c8e901aa8d37eb23bf9d332d63049ea9e  published ablation CSV
1ed6c2a5784404544484201f92567fbe2485de608445956871a8c0d53b895e23  source ablation script
5a02bcfcb9c2f9f5ac3440b0304713c27da68719a85408051ee704b3b9d83a0d  source calibrated JLD
55527e57d4b9e2faace3cf8ad2413892229da80daab7ac473dd9549045f8118f  reused public Figure 2 sweep.jld2
```

[`SHA256SUMS`](SHA256SUMS) verifies the published CSV.
The Figure 2 dataset has 75 fields. Its metadata/gain file retains all
69 non-trajectory fields; the six trajectory fields are stored separately.
All fields were checked with `isequal` against the source container. Nothing from that approximately 500 MB archive is duplicated here.

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
percentages. Mean SD is not sqrt(mean variance). The physical-variance
negativity threshold is −1e-10 before clipping roundoff to zero.

The reported reduction ranges refer to **nonzero** internal noise. Analytical
full-model endpoint SD is 8.63464–10.54364%; the approximately 8.55–10.46%
Figure 2 values are archived Monte Carlo sample SDs. Analytical and Monte Carlo results are reported separately.

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
