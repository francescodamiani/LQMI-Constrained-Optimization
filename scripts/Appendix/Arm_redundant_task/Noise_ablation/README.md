# Stabilization variability and noise-source ablations

Measure behavioral variability, remove noise sources with fixed full-model
gains, and compare TOD/LCO when sensory and motor noise are increased.

## Run

With Julia **1.8.5**, run from the repository root:

```bash
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl --simulate --tag=my_run
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl --data=data/generated/Appendix/Arm_redundant_task/Noise_ablation/my_run
```

The default prints tables from `data/published/Appendix/Arm_redundant_task/Noise_ablation/`
without simulation or figures. Reruns save to
`data/generated/Appendix/Arm_redundant_task/Noise_ablation/TAG/`.
`--experiment=behavioral`, `--experiment=original_ablation` and
`--experiment=high_noise` select individual analyses; use the same selection
when reading a partial run. `--reoptimize` with `--simulate` also recomputes
baseline full-model policies. Use `--help` for all options.

Published JLD2 datasets:

| File | Contents |
|---|---|
| `behavioral_variability.jld2` | 20-level component variances, analytical/Monte Carlo summaries and relative changes |
| `original_fixed_policy_ablation.jld2` | Full/ablated x/z variances and relative-variance reductions |
| `original_full_model_policies.jld2` | All 20 TOD/LCO K/L sequences, grid and completion flags |
| `external_noise_robustness.jld2` | Multipliers 1, 5 and 10: gains, costs, variances, ablations and histories; default tables select 10 |

## Settings and definitions

The [arm model](../README.md) uses 6 states, 9 controls, 3 angle observations,
T=300, Δt=0.01 s and 500 TOD/LCO iterations/passes. LCO starts from TOD with
regularization `[0]`, fallback `[1e-10,1e-8,1e-6]`.

- **Behavioral variability:** external amplitudes 0.1, 20 internal-noise levels
  in `[0,0.4]`, 5,000 Monte Carlo trials and seed 1234. Component SD is
  `sqrt(mean_t Var_trials(x_i,t))` over `1:298`, with unbiased trial variance.
  Angles and velocities are reported separately in rad and rad/s.
- **Baseline ablation:** exact centered moments with K/L frozen. Combined
  sensory removal zeros both additive and multiplicative sensory noise.
  Component SD reduction is `100*(1-SD_ablated/SD_full)`; group reductions
  compare mean SD across the three angles or three velocities.
- **High-noise comparison:** sensory/motor amplitudes 1.0, process amplitude
  0.1 and internal SD `[0,.1,.2,.3,.4]`. Full-model gains are optimized per
  condition; 2,000 Monte Carlo trials use seed 1234. Analytic ablations remove
  the two sensory channels separately, without reoptimizing gains.

Ablation percentages are not additive variance shares; separate sensory
reductions cannot be summed into a measured combined effect. Raw saved
percentage fields refer to variance; printed SD tables use square roots.
High-noise LCO cost improvements are 0%, 8.7661%, 25.2846%, 40.2810% and 52.1711%.
Raw Monte Carlo trajectories are not included; seeded reruns need not reproduce
saved sample statistics exactly. Analytical ablations can be recomputed
from the included full-model gains.

`noise_ablation.jl` is the entry point, `noise_statistics.jl` loads data and
prints tables, `noise_experiments.jl` runs experiments, and
`component_variances.jl` propagates exact moments.
