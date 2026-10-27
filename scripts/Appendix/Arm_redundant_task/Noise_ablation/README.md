# Stabilization variability and fixed-policy noise-source ablations

These numerical analyses concern the zero-mean
three-joint stabilization task. Their purpose is to compare the behavioral
effects of process, sensory, motor and internal noise, and verify the LCO
advantage when sensory/motor amplitudes are increased. They are separate
from both the under-observed arm and calibrated reaching task.

## Print the saved numerical results

With Julia **1.8.5**, run from the repository root:

```bash
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl
# Print an individual experiment:
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl --experiment=behavioral
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl --experiment=original_ablation
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl --experiment=high_noise
```

Default commands load complete saved numerical data and print Markdown
tables. They run no optimization, moment propagation or Monte Carlo, load no
plotting packages and write no files. No figures or duplicate table exports
are generated. Results are under
`data/published/Appendix/Arm_redundant_task/Noise_ablation/`.

| Script | Responsibility |
|---|---|
| [noise_ablation.jl](noise_ablation.jl) | Entry point, experiment selection, archived-data default and explicit generated runs. |
| [noise_statistics.jl](noise_statistics.jl) | Loads numerical archives and computes SDs, changes, component ranges, mean-SD reductions and cost tables. |
| [noise_experiments.jl](noise_experiments.jl) | Baseline-regime analysis reruns, high-noise full-model optimization, fixed-gain ablations, validated-core Monte Carlo and numerical storage. |
| [component_variances.jl](component_variances.jl) | Exact first/second-moment analysis helper. It is not a separate stochastic simulator. |

The existing [arm_task.jl](../arm_task.jl) supplies the plant, TOD/LCO and
general analytical cost; [benchmark_tables.jl](../../benchmark_tables.jl)
supplies shared paths and Markdown printing.

## Parameters and statistical definitions

All experiments use m=6, p=9, k=3, T=300, Δt=0.01, the arm inertia,
damping, nine-actuator mapping and Euler dynamics. H observes all three angles.
Initial physical/internal means and centered covariances are zero. State cost
is `diag(1,1,1,0.001,0.001,0.001)` throughout; control cost is 0.01I with zero
terminal slice. TOD uses 500 iterations. LCO starts from the matching TOD
policy and uses 500 passes, one L/K update each, rho=[0], fallbacks
[1e-10,1e-8,1e-6], and both internal GD routes disabled.

The baseline regime uses all external amplitudes **0.1** and the
20-point `round.(range(0,0.4,length=20),digits=2)` internal-noise grid.
Behavioral variability uses **5,000 trials per method/condition**, seed 1234
reset for each comparison. Each component's across-trial **sample variance**
(N−1 denominator) is computed separately at each time, then averaged over
**1:298**. Reported SD is the square root of that average:

```text
V_i = mean_t Var_trials(x_i,t)
SD_i = sqrt(V_i)
SD change (%) = 100 * (SD_i(sigma_eta)/SD_i(0) - 1)
```

It is not the average of instantaneous SDs or a variance pooled across time.
The window includes the zero-variance initial state and excludes the last
two samples. Angles (components 1:3) and velocities (4:6) remain separate,
in rad and rad/s. The full LCO sweep gives endpoint changes
**15.50%, 12.07%, 10.95% for joints 1,2,3**, respectively.

The leave-one-out analysis uses **exact centered moment propagation**,
not additional Monte Carlo. At each internal-noise level, the complete-model
LCO gains are frozen for all ablations. Combined sensory removal zeros both
additive sensory covariance and multiplicative sensory matrix D. Motor removal
zeros C; process/internal removal zeros only its own additive covariance.
H, K and L are unchanged. Raw archived percentage fields are **variance**
reductions; the displayed SD reductions are derived as:

```text
Component reduction (%) = 100 * (sqrt(V_full) - sqrt(V_ablated)) / sqrt(V_full)
Group reduction (%) = 100 * (mean_i SD_full,i - mean_i SD_ablated,i) / mean_i SD_full,i
```

Group means are over the three angles or three velocities separately. They
are not means of component percentages or square roots of component-averaged
variances. Sweep changes include policy adaptation; fixed-policy ablations
measure a different counterfactual and need not match them.

The alternative regime uses sensory/motor amplitudes **1.0**, while process
amplitude stays **0.1**, at internal SD [0,0.1,0.2,0.3,0.4]. Additive sensory
covariance therefore increases 100-fold. TOD and LCO are independently
recomputed at every full-model noise level, then their policies are frozen
for each ablation. The archived analysis removes **additive sensory and
multiplicative sensory channels separately**, alongside process, motor and
internal noise; it does not contain a combined sensory removal. These separate
percentages must not be summed and described as a measured combined effect.
Full-model MC summaries use 2,000 trials and seed 1234; ablations remain analytic.

LCO cost improvements are **0%, 8.7661%, 25.2846%, 40.2810%, 52.1711%**.
Across nonzero internal noise, mean-angle SD reductions are 79.9968–80.9654%
for process, 1.2687–1.5797% for additive sensory, 0.5656–0.6264% for
multiplicative sensory, and **0.4304–0.5614%** for internal noise. Sensory
channels each exceed the internal effect in these conditions. Ablation effects
are closed-loop SD changes, not additive fractions of variance.

## Optional numerical reruns

```bash
# Rerun selected analyses using their publication settings:
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl --simulate --tag=new_noise
# Rerun only the analytic fixed-policy ablation (no MC or optimization):
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl --simulate --experiment=original_ablation --tag=new_ablation
# Also recompute baseline-regime full-model policies when desired:
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl --simulate --reoptimize --experiment=original_ablation --tag=new_policies
# Small numerical check; no figures or publication claims:
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl --simulate --validate --tag=check_noise
# Inspect generated results without rerunning:
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl --data=data/generated/Appendix/Arm_redundant_task/Noise_ablation/new_noise
```

The baseline analyses reuse the full-model policy cache
unless `--reoptimize` is explicit. High-noise reruns always recompute the
matching full-model TOD/LCO solutions independently. No ablated policy is
ever reoptimized. Newly used baseline-regime policies are stored once in the
generated run so both analyses share the exact cache.

`--validate` uses T=12, internal SD [0,0.1], three TOD/LCO passes and 64 MC
trials; its archives and printed tables are explicitly marked reduced.
Generated results are confined to `data/generated/Appendix/Arm_redundant_task/Noise_ablation/TAG/`.
Existing tags require `--force`; no published file is overwritten. When
loading a run containing only one analysis, select that same experiment.

New trajectories use only `monte_carlo_simulation`. The seed value and reset
convention are preserved. The draw order and RNG used for the saved trials
can differ from the core's seeded draw order, so identical MC realizations
are not guaranteed. Raw trajectories/noise draws were not
archived, but all numerical summaries are preserved losslessly.
