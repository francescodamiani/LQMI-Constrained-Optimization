# Arm variability and noise-ablation data

Four JLD2 files contain the complete numerical results:

| File | Contents |
|---|---|
| `behavioral_variability.jld2` | 20-level TOD/LCO component variances, analytical/Monte Carlo summaries and relative-variance changes |
| `original_fixed_policy_ablation.jld2` | Full and ablated x/z variances for both methods, conditions and relative-variance reductions |
| `original_full_model_policies.jld2` | All 20 TOD/LCO K/L gain sequences, grid and completion flags |
| `external_noise_robustness.jld2` | Sensory/motor multipliers 1, 5 and 10: gains, costs, variances, ablations, histories and settings; default tables select 10 |

Run `julia scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl`
from the repository root to print saved tables. The
[experiment guide](../../../../../scripts/Appendix/Arm_redundant_task/Noise_ablation/README.md)
gives selections and rerun commands. Generated results go to
`data/generated/Appendix/Arm_redundant_task/Noise_ablation/TAG/` and are read
with `--data`; no figures are produced.

Behavioral SD is `sqrt(mean_t Var_trials(x_i,t))` over `1:298`, with unbiased
trial variance. Saved relative fields describe variance; SD tables apply
square roots first. Ablations freeze full-model gains. Baseline sensory
removal combines both sensory channels; the high-noise analysis removes them
separately. Percentages are not additive variance shares.

Raw Monte Carlo trajectories/draws are not included, so exact saved sample
statistics cannot be reconstructed without rerunning. Analytical ablations
can be recomputed from the included gains.
