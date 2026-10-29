# Reaching: internal noise and alternative explanations

Compare internal-noise signatures with sensory-model mismatch and alternative
noise sources. Results apply to the tested models and ranges; they do not
establish general identifiability of internal noise from behavior alone.

## Run

With Julia **1.8.5**, run from the repository root:

```bash
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl --full
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl --simulate --experiment=overestimation,underestimation --tag=my_run
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl --experiment=overestimation,underestimation --data=data/generated/Appendix/Reaching_identifiability/my_run
```

The default prints saved tables from `data/published/Appendix/Reaching_identifiability/`;
`--full` includes every comparison. Reruns save to
`data/generated/Appendix/Reaching_identifiability/TAG/` and may be expensive.
No figures are generated. `--experiment=NAME` selects `overestimation`,
`noise_sources` or `underestimation`. Use `--help` for options.

Available JLD2 datasets:

| File/directory | Contents |
|---|---|
| `sensory_overestimation.jld2` | A: all 21 assumed-noise conditions, gains, costs and behavioral summaries |
| `matched_noise_sources.jld2` | B: 100 comparisons, 80 reoptimized external-noise policies and extended statistics |
| `candidate_noise_grids.jld2` | All 125 candidate policies, costs, variability and completion flags |
| `sensory_underestimation.jld2` | C: ten assumed-noise and twenty internal-noise conditions, gains, costs and summaries |
| `reference/` | B reference model/gains/settings in `sweep.jld2`, plus full TOD/LCO trajectories by noise level |

A/C use the [main Figure 2](../../Main/Figure2/README.md) model and saved
internal-noise sweep from `data/published/Main/Figure2/`. B alone uses the
separate `reference/` dataset. Both references have 20 internal-noise levels
and 2,000 trials per method/level, with `(dimension,trial,time)` trajectory
axes. Default tables read summary files without loading trajectories.

## Settings

A/C obtain dynamics, initial conditions and costs directly from
`Main/Figure2/reaching_task.jl`: 6 states, 9 controls, 3 observations,
41 stored states at Δt=0.01 s, and a 20° reach over 0.40 s. External amplitudes
are `σ_ξ=0.01`, `σ_ω=deg2rad(0.9)`, `σ_ϵ=0.10`, `σ_ρ=0.01`.
Both use the same 20 internal levels in `[0,0.04]` as Figure 2. Saved A/C files
include a `model` dictionary recording the complete calibrated definition.

Each optimized condition uses 500 TOD iterations and 500 LCO cycles from TOD,
with regularization `[0]`, fallback `[1e-10,1e-8,1e-6]` and no internal GD.
Monte Carlo seed is 1; A uses 2,000 trials and C retains 1,000. Internal-noise
policies and trials are reused from Figure 2; C uses the first 1,000 trials.
The correctly specified, zero-internal-noise gains are shared by A/C.

- **A — sensory-model mismatch:** true sensory amplitude `deg2rad(0.9)`,
  internal noise zero. The existing 21-point logarithmic grid spans
  `0.1–10 ×` the true amplitude (approximately `0.001571–0.157080` rad).
  Default tables show the baseline and overestimation conditions; `--full`
  also includes the lower half of this same grid. Gains are optimized using
  assumed sensory noise and then frozen for evaluation with true noise.
  Internal comparisons interpolate saved statistics only within overlapping
  analytical endpoint-SD ranges, without extrapolation or additional fitting.
- **B — alternative noise sources:** this separate experiment retains its
  own model, results and reference, with all four baseline external amplitudes
  at `.01`. It independently optimizes each external-noise source at zero
  internal noise, then reoptimizes at amplitudes matched to 20 internal-reference
  endpoint-SD targets. Candidate grids have 31 process, 32 additive-sensory
  and 31 levels per multiplicative source, spanning `[0,.1]`, `[.0001,.1]`
  (plus .01), and `[0,.3]`, respectively. It uses 2,000 Monte Carlo trials.
- **C — sensory underestimation:** the same true calibrated system as A;
  ten linearly spaced assumed amplitudes decrease from `1 ×` to `0.32 ×`
  the true value (approximately `0.015708–0.005027` rad), with internal noise
  zero. Gains are frozen for true-system evaluation. The comparison uses
  Figure 2's correctly specified internal-noise policies at all 20 levels.
  These are raw sweeps without endpoint matching.

## Interpretation

Endpoint variability is centered terminal elbow SD as a percentage of 20°;
peak velocity is the mean per-trial maximum absolute elbow velocity, with SEM.
A/C use Figure 2's gain log-pseudodeterminants: sum natural logs of singular
values above `1e-6`, then average indices `2:38`. Frobenius summaries are also
retained. B's gain threshold remains `1e-6`. Costs are evaluated analytically
in the true system; percentage changes use the correctly specified,
zero-internal-noise baseline. Assumed-model costs are stored separately.
B's PCA uses centered trials over `2:39` and default 99% variance retention.

The multiplicative grids each fail to bracket 13/20 targets; 12/20 achieved
comparisons per grid exceed 1% endpoint-SD error. Unbracketed nearest-grid
results are not successful matches. Some gain/velocity trends have small
reversals. The A–C summary files omit per-condition trajectories and
optimization histories; their referenced internal-sweep trials and histories
are included in the respective reference datasets.
Monte Carlo reruns need not match saved draws exactly.

## Saved calibrated results

Overestimation up to `10 ×` true sensory noise increases analytical endpoint SD
from 8.6346% to 12.4012% of the reach and true cost by 10.0371%; mean peak
velocity decreases overall from 76.7221 to 74.4990 deg/s. Eight of the eleven
baseline/overestimation conditions overlap the internal-reference endpoint-SD
range (8.6346–10.5436%); the three highest amplitudes remain unmatched.

Underestimation down to `0.32 ×` true sensory noise decreases endpoint SD to
8.2445% and increases true cost by 1.1624%; mean peak velocity increases overall
from 76.9706 to 77.5626 deg/s. Its nonbaseline endpoint SDs lie below the
internal-reference range. Both mismatch sweeps have small sampled peak-velocity
reversals near baseline. The filter-gain metric decreases with overestimation
and increases with underestimation; the control-gain metric changes little.
Internal noise instead produces a non-monotonic peak-velocity curve and a much
larger decrease in the control-gain metric. These conclusions concern the
sampled conditions, not general identifiability.

`reaching_identifiability.jl` is the entry point; `identifiability_statistics.jl`
prints comparisons, and `identifiability_simulation.jl` dispatches
`sensory_overestimation.jl`, `matched_noise_sources.jl` and
`sensory_underestimation.jl`.
