# Stabilization variability and noise-ablation data

These four archives retain **all 82 fields** from the relevant
source JLD files. They were converted to JLD2 solely to use the existing pinned
Julia 1.8.5 environment without adding JLD/HDF5 dependencies. The conversion
verified identical key sets and `isequal(source[k],stored[k])` for every
field, including complete arrays, statuses, histories and metadata. The complete numerical records are included.
`SHA256SUMS` verifies the public JLD2 files.

| Public JLD2 file | Source container | Contents |
|---|---|---|
| `behavioral_variability.jld2` | `Data/Arm_Redundant_Task/Internal_Noise_Component_Variances_all_time_steps.jld` | 19 fields: 20-level TOD/LCO physical/internal MC and analytical component variances, all saved relative-variance fields, window and 5,000-trial count |
| `original_fixed_policy_ablation.jld2` | `Data/Arm_Redundant_Task/Leave_One_Noise_Source_Out/Leave_One_Noise_Source_Out_Component_Variances.jld` | 10 fields: full and four ablated x/z variance arrays for both methods, names, grid/window and raw relative-variance reductions |
| `original_full_model_policies.jld2` | `Data/Arm_Redundant_Task/Internal_Noise_Optimal_Solutions.jld` | 6 fields: all 20 TOD/LCO K/L sequences, grid and completion flags |
| `external_noise_robustness.jld2` | `results/arm_redundant_task_external_noise_robustness_NEW_TEST/arm_redundant_task_external_noise_robustness_results.jld` | 47 fields: all s=1,5,10 conditions, policies, costs, MC/analytic variances, individual ablations, full histories, statuses and settings; default tables select s=10 |

Together these are about **11.7 MiB**; the largest file is 6.6 MiB. The archives supply all default tables without duplicate CSV summaries or plots.

Source-container SHA-256 checksums, in the table's order:

```text
b5012034d9593ac3c0906e5aeb6b9e3d9bb078d2e3752193563329636793b0bb
ed22f8ee2b62946a08a742884db43dfcd01575226686236ae093194e63342d8a
a41f4a0eaa9c6b429a0b9bba0aec64e9098e82c0d42b6236b5b1c46f6b2a68d4
8b8159faa1d91b85da15ab576e3af1a16950b031733e2d51163f24ca3a198962
```

```bash
julia scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl
```

This prints saved-data tables and writes nothing. The behavioral SD is
`sqrt(mean_t Var_trials(x_i,t))`, using unbiased across-trial variances and
window **1:298**. The stored relative fields describe variance,
while reported SD changes/reductions are derived after square roots.
The baseline sensory ablation removes both sensory channels; s=10 stores
them separately. K/L remain the complete-model gains during ablations.

Raw MC trajectories and raw random draws are not present. Exact
5,000-trial sample statistics cannot be independently rebuilt from trajectories
without new draws, but every archived variance and derived table is retained.
Analytic ablations can be independently recomputed from the complete policy
cache. New simulations save only under the corresponding
`data/generated/Appendix/Arm_redundant_task/Noise_ablation/TAG/` directory and
are clearly marked generated/reduced. Their RNG realizations may differ.
