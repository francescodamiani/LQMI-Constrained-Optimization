# Reaching: internal noise and alternative explanations

Compare internal-noise signatures with sensory-model mismatch and alternative
noise sources. Results apply to the tested models and ranges; they do not
establish general identifiability of internal noise from behavior alone.

## Run

With Julia **1.8.5**, run from the repository root:

```bash
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl --full
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl --simulate --tag=my_run
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl --data=data/generated/Appendix/Reaching_identifiability/my_run
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
| `sensory_underestimation.jld2` | C: ten assumed-noise and thirteen internal-noise conditions, gains, costs and summaries |
| `reference/` | Shared model/gains/settings in `sweep.jld2`, plus full TOD/LCO trajectories by noise level |

The reference has 20 internal-noise levels and 2,000 trials per method/level;
trajectory axes are `(dimension,trial,time)`. A/B reuse this reference. Default
tables read summary files without loading the large trajectories.

## Settings

The reaching model has 6 states, 9 controls, 3 observations, T=41, Δt=0.01 s
and a 20° reach. Costs and initial conditions are stored in `reference/sweep.jld2`.
Baseline external-noise amplitudes are
**`σ_ξ=σ_ω=σ_ϵ=σ_ρ=0.01`**. These experiments use different noise settings from
[main Figure 2](../../Main/Figure2/README.md), which uses `σ_ξ=0.01`,
`σ_ω=deg2rad(0.9)`, `σ_ϵ=0.10`, `σ_ρ=0.01`.

Each optimized condition uses 500 TOD iterations and 500 LCO cycles from TOD,
with regularization `[0]`, fallback `[1e-10,1e-8,1e-6]` and no internal GD.
Monte Carlo seed is 1; A/B use 2,000 trials and C uses 1,000.

- **A — sensory-model mismatch:** true sensory amplitude 0.01, internal noise
  zero, 21 assumed amplitudes from 0.001 to 0.1. Gains are optimized under
  assumed noise and evaluated under true noise. Comparisons interpolate the
  internal-noise reference only within overlapping endpoint-SD ranges.
- **B — alternative noise sources:** independently optimize each external-noise
  source at zero internal noise, then reoptimize at amplitudes matched to
  20 internal-reference endpoint-SD targets. Candidate grids have 31 process,
  32 additive-sensory and 31 levels per multiplicative source, spanning
  `[0,.1]`, `[.0001,.1]` (plus .01), and `[0,.3]`, respectively.
- **C — sensory underestimation:** true sensory amplitude 0.05, ten assumed
  values from .05 to .016 at zero internal noise, compared with correct-model
  internal noise `0:.0025:.03`. Other external amplitudes stay at .01.
  These are raw sweeps without endpoint matching.

## Interpretation

Endpoint variability is centered terminal elbow SD as a percentage of 20°;
peak velocity is the mean per-trial maximum absolute elbow velocity, with SEM.
Gain metrics average indices `2:38`; log-volume singular-value thresholds are
`1e-6` for B and `1e-4` for C. Cost changes use each experiment's
correct-model baseline. B's PCA uses centered trials over `2:39` and default
99% variance retention.

The multiplicative grids each fail to bracket 13/20 targets; 12/20 achieved
comparisons per grid exceed 1% endpoint-SD error. Unbracketed nearest-grid
results are not successful matches. Some gain/velocity trends have small
reversals. The A–C summary files omit per-condition trajectories and
optimization histories; the shared internal-reference trials and histories
are included.
Monte Carlo reruns need not match saved draws exactly.

`reaching_identifiability.jl` is the entry point; `identifiability_statistics.jl`
prints comparisons, and `identifiability_simulation.jl` dispatches
`sensory_overestimation.jl`, `matched_noise_sources.jl` and
`sensory_underestimation.jl`.
