# Numerical datasets

`published/` contains the numerical datasets accompanying the paper.
`generated/` holds user simulation results and is Git-ignored.

| Dataset guide | Contents |
|---|---|
| [Figure 1](published/Main/Figure1/README.md) | Dimension benchmark, convergence histories and internal-noise sweep |
| [Figure 2](published/Main/Figure2/README.md) | Calibrated reaching model, gains, costs and full Monte Carlo trajectories |
| [LCO perturbations](../scripts/Appendix/RandomTasks/README.md) | 2.4 million local-perturbation cost evaluations |
| [Under-observed random tasks](published/Appendix/RandomTasks/Under_observed/README.md) | Costs, timings and LCO histories |
| [Initialization sensitivity](published/Appendix/RandomTasks/Initialization_sensitivity/README.md) | All 202 starts, final gains and cost histories |
| [Redundant-arm task](published/Appendix/Arm_redundant_task/README.md) | Data availability and simulation instructions |
| [Under-observed arm](published/Appendix/Arm_redundant_task/Under_observed/README.md) | Policies, objectives and histories |
| [Arm noise ablation](published/Appendix/Arm_redundant_task/Noise_ablation/README.md) | Behavioral variability, fixed-policy ablations and external-noise comparisons |
| [Reaching identifiability](published/Appendix/Reaching_identifiability/README.md) | Sensory-model mismatch and alternative noise-source comparisons |
| [Reaching noise ablation](published/Appendix/Reaching_noise_ablation/README.md) | Fixed-gain ablations using the Figure 2 policies |

Main-figure commands load saved data without simulation. Arm figure trajectories
are not bundled; that workflow runs simulations and can optionally save a cache.
The complete numerical download is about 1 GB.

Experiment READMEs give the `--simulate`, `--tag` and `--data` commands.
Relative `--data` paths resolve from the repository root. Plotting commands
write PDFs under `figures/generated/`; numerical table commands print to the
terminal. Published data are never overwritten. Seeded Monte Carlo reruns can
use a different draw order from the saved trials and need not match them exactly.
