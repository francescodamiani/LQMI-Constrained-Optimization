# Reaching-identifiability data

| File/directory | Contents |
|---|---|
| `sensory_overestimation.jld2` | 21 assumed sensory-noise conditions, gains, true/assumed costs and behavioral summaries |
| `matched_noise_sources.jld2` | 100 comparison rows, 80 reoptimized external-noise policies, matched amplitudes and extended statistics |
| `candidate_noise_grids.jld2` | All 125 candidate policies, costs, variability and completion flags |
| `sensory_underestimation.jld2` | Ten assumed-noise and thirteen internal-noise conditions, gains, costs and behavioral summaries |
| `reference/` | Shared model/gains/settings in `sweep.jld2`, plus full TOD/LCO trajectories by internal-noise level |

The reference has 20 internal-noise levels and 2,000 trials per method/level.
Trajectory axes are `(dimension,trial,time)`. A/B reuse this reference; C uses
1,000 trials. The datasets total about 485 MiB.

Baseline external amplitudes are `σ_ξ=σ_ω=σ_ϵ=σ_ρ=0.01`; C sets true sensory
amplitude to 0.05 while other external amplitudes stay at .01. These settings
are specific to this experiment; main Figure 2 uses its calibrated amplitudes.

Run `julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl`
from the repository root; `--full` prints all comparisons. Default tables use
summary files, without loading large trajectories. See the
[experiment guide](../../../../scripts/Appendix/Reaching_identifiability/README.md)
for metrics, limitations and rerun commands. Generated results go to
`data/generated/Appendix/Reaching_identifiability/TAG/` and are read with `--data`;
no figures are produced.

Multiplicative grids do not bracket all endpoint targets; nearest-grid
fallbacks are not successful matches. Condition-level mismatch/matched
trajectories and optimization histories are not included; shared-reference
trials and histories are included.
