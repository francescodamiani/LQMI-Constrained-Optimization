# Behavioral variability and noise-source contributions

Three stabilization analyses, ordered as in the paper's discussion. The
[arm model](../README.md) has 6 states, 9 controls, 3 angle observations,
T=300, Δt=0.01 s and zero initial means/covariances. Optimization uses
500 TOD iterations and 500 LCO cycles, initialized from TOD, with
regularization `[0]` and fallback `[1e-10,1e-8,1e-6]`.

## Reproduce saved results

With Julia **1.8.5**, run from the repository root:

```bash
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/behavioral_variability.jl
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/baseline_noise_ablation.jl
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/high_sensory_motor_noise.jl
```

Each command prints its numerical tables from
`data/published/Appendix/Arm_redundant_task/Noise_ablation/`, without simulations
or figures. `noise_ablation.jl` runs all three in the same order.

## 1. Behavioral variability

`behavioral_variability.jl` reports physical-state SD under LCO across 20
internal-noise levels in `[0,0.4]`, with external amplitudes 0.1 and 5,000
Monte Carlo trials per condition (seed 1234). Component SD is
`sqrt(mean_t Var_trials(x_i,t))` over `1:298`, using unbiased trial variance.
Angles and velocities are reported separately in rad and rad/s; components
1–3 are shoulder, elbow and wrist. Their SD increases from zero to maximum
internal noise are 15.50%, 12.07% and 10.95%, respectively.

## 2. Baseline noise-source ablations

`baseline_noise_ablation.jl` uses exact centered moments at the same 20 noise
levels. Full-model LCO gains remain fixed while removing process, internal
or multiplicative motor noise individually, or both sensory channels together.
It reports component SD-reduction ranges and reductions in mean component SD
for angles and velocities separately.

Component reduction is `100*(1-SD_ablated/SD_full)`; group reduction compares
the mean of the three component SDs before and after removal. All use the
same `1:298` window. These effects are not additive fractions of total variance.
Process-noise removal reduces component SD by about 81–98%; maximum internal
effects are 1.83% for angles and 1.70% for velocities.

## 3. Higher sensory/motor-noise regime

`high_sensory_motor_noise.jl` selects tenfold increases in additive sensory,
multiplicative sensory and multiplicative motor amplitudes (0.1 to 1.0),
with process amplitude fixed at 0.1 and internal SD `[0,.1,.2,.3,.4]`.
TOD is recomputed at each level; LCO starts from that TOD solution. Analytical
costs give LCO improvements of 0%, 8.77%, 25.28%, 40.28% and 52.17%.

Fixed-gain analytical ablations remove additive and multiplicative sensory
noise **separately**, along with the other noise sources. The mean joint-angle
SD reductions at nonzero internal noise are 1.27–1.58% (additive sensory),
0.57–0.63% (multiplicative sensory) and 0.43–0.56% (internal).
Separate sensory reductions cannot be summed into a combined-removal effect.
The dataset also includes 2,000-trial Monte Carlo summaries (seed 1234).

## Explicit simulation reruns

Every entry script accepts `--simulate`, `--tag`, `--data` and `--help`. For example:

```bash
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/behavioral_variability.jl --simulate --tag=my_run
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/behavioral_variability.jl --data=data/generated/Appendix/Arm_redundant_task/Noise_ablation/my_run
```

Reruns save numerical results under
`data/generated/Appendix/Arm_redundant_task/Noise_ablation/TAG/`. Baseline
analyses reuse saved full-model gains; adding `--reoptimize` also recomputes
those gains. Higher-noise reruns always recompute TOD/LCO. Read a partial run
with the same entry script that generated it. The combined runner also accepts
`--experiment=behavioral`, `--experiment=baseline_ablation` or `--experiment=high_noise`;
existing commands remain supported.

## Data and shared scripts

Published JLD2 files are loaded with `JLD2.load(path)`:

| File | Contents |
|---|---|
| `behavioral_variability.jld2` | Component variances, analytical/Monte Carlo summaries and relative changes |
| `original_fixed_policy_ablation.jld2` | Full/ablated x/z variances and relative-variance reductions |
| `original_full_model_policies.jld2` | All 20 TOD/LCO K/L sequences, noise grid and completion flags |
| `external_noise_robustness.jld2` | Multipliers 1, 5 and 10: gains, costs, variances, ablations and histories; this analysis selects 10 |

Raw percentage fields refer to **variance**; the printed SD tables take square
roots first. Raw Monte Carlo trajectories are not included; seeded reruns need
not reproduce saved sample statistics exactly. Analytical ablations can be
recomputed from saved gains.

The three entry scripts select analyses through `noise_ablation.jl`.
`noise_statistics.jl` loads data and prints tables; `noise_experiments.jl`
runs simulations; `component_variances.jl` propagates exact moments.
