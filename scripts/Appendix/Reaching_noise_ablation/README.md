# Calibrated reaching: fixed-gain noise ablation

This appendix analysis supports main Figure 2: internal noise can affect the
optimal policy while its direct effect on endpoint variability remains smaller
than process disturbances. All ablations freeze the full-model K/L.

With Julia 1.8.5 and the root pinned environment:

```bash
# Saved numerical tables; no computation or file output:
julia scripts/Appendix/Reaching_noise_ablation/reaching_noise_ablation.jl
# Recompute all 35 analytical cases with existing Figure 2 gains:
julia scripts/Appendix/Reaching_noise_ablation/reaching_noise_ablation.jl --simulate --tag=my_run
# Read that run without recomputing:
julia scripts/Appendix/Reaching_noise_ablation/reaching_noise_ablation.jl --data=data/generated/Appendix/Reaching_noise_ablation/my_run
```

`--simulate` performs analytical propagation only: no optimization, Monte Carlo
or figures. New results are one JLD2 under `data/generated/Appendix/Reaching_noise_ablation/TAG/`;
existing tags require `--force`. The complete calculation is inexpensive, so
there is no reduced mode.

## Model and source removals

The unchanged [Figure 2 model](../../Main/Figure2/reaching_task.jl) and
`data/published/Main/Figure2/sweep.jld2` supply the calibrated task and gains.
No large trajectory archive is loaded or duplicated. Model matrices and each
full-model objective must agree with the saved record before recomputation.

The 20° elbow reach starts at −π/9, with 6 states, 9 controls, 3 angle observations,
Δt=.01 s and T=41. Noise amplitudes are ξ=.01, ω=deg2rad(.9), ε=.10, ρ=.01:
Ωξ=σξ²I₆, Ωω=σω²I₃, Ωη=ση²I₆, C=σεB, D=σρH. The same costs are used for all ablations.
Nominal targets `[0,.01,.02,.03,.04]` select actual levels
**[0,.011,.019,.029,.040]**, indices **[1,6,10,15,20]**; no gains are interpolated.

| Condition | Removed terms |
|---|---|
| `full_noise` | None |
| `no_process_noise` | Ωξ |
| `no_additive_sensory_noise` | Ωω |
| `no_all_sensory_noise` | **Ωω and D together** |
| `no_motor_multiplicative_noise` | C |
| `no_sensory_multiplicative_noise` | D |
| `no_internal_noise` | Ωη |

All other terms, H and the deterministic estimator A+BL−KH stay fixed.
The existing exact `component_variances.jl` helper is reused with singleton
windows to recover instantaneous centered variances without another moment
algorithm or a core change. Physical variance below −1e-10 is rejected;
accepted roundoff negatives are clipped before taking square roots.

## Metrics and qualifications

- Terminal elbow SD: `100*sqrt(Var(x[2,T]))/(π/9)`, percent of the 20° reach.
- Mean-angle SD: `mean(sqrt.(variance_x[1:3,2:41]))`, **rad**.
- Mean angular-velocity SD: `mean(sqrt.(variance_x[4:6,2:41]))`, **rad/s**.
- Each reduction: `100*(1-SD_ablated/SD_full)` for that same policy and metric.

The means average instantaneous SD across three joints and 40 samples, not
sqrt(mean variance). Combined sensory removal is evaluated simultaneously;
its percentage is not a sum of the individual effects. Counterfactual SD
changes are not additive variance shares.

For nonzero η, terminal-SD reductions range from **65.3551–78.2432%** for
process, **1.1623–5.1861%** for internal, **.3839–.7804%** for combined sensory,
and **.6256–.7310%** for motor noise. At η=0 they are **83.7861%, 0%, .9207%,
.5133%**, respectively; the narrower reported ranges exclude this baseline.
Analytical full-model endpoint SD is **8.63464–10.54364%**. Figure 2's
**8.55236–10.45508%** uses archived Monte Carlo sample SD, explaining the
finite-sample difference without a calibration change.

`reaching_noise_ablation.jl` is the entry point; `ablation_tables.jl` loads
CSV/JLD2 and prints all six metric tables; `fixed_gain_ablation.jl` checks the
model, evaluates fixed policies and saves results. See the
[dataset README](../../../data/published/Appendix/Reaching_noise_ablation/README.md)
for the complete 35-row CSV, units, hashes and provenance.
