# LCO local-perturbation experiments

Test changes in expected cost around an LCO solution under perturbations of
control gains L, estimation gains K, and their joint vector.

## Run

With Julia **1.8.5**, run from the repository root:

```bash
julia scripts/Appendix/RandomTasks/perturbations.jl
julia scripts/Appendix/RandomTasks/perturbations.jl --simulate --tag=my_run
julia scripts/Appendix/RandomTasks/perturbations.jl --data=data/generated/Appendix/RandomTasks/my_run
```

The default reads the complete saved scan in `data/published/Appendix/RandomTasks/`
without optimization. The compressed JLD2 file is loaded automatically; allow
about 307 MB temporary space. Explicit simulations save to
`data/generated/Appendix/RandomTasks/TAG/` and evaluate **2.4 million costs**.

Three PDFs are written to `figures/generated/Appendix/RandomTasks/`:
`Fig_LCO_local_perturbation_violin_deltaC_L_only.pdf`,
`Fig_LCO_local_perturbation_violin_deltaC_K_only.pdf`, and
`Fig_LCO_local_perturbation_violin_deltaC_joint.pdf`.
They show red violins and median overlays on a logarithmic cost-increase axis.

## Settings and interpretation

The task has n=5, T=10 and internal noise 0.2, with 500 TOD iterations and
500 LCO cycles. Regularization is `[0]`, with fallback `[1e-10,1e-8,1e-6]`;
internal GD routes are disabled. Each perturbation family uses 100,000 normalized
Gaussian directions, both signs and relative amplitudes `[1e-4,1e-3,1e-2,1e-1]`.
Displacements scale by the baseline L, K or joint-vector norm. Costs are
analytical. All saved cost increments are positive: evidence from sampled
directions, not a proof of local optimality.

`perturbations.jl` is the entry point; `perturbation_helpers.jl` defines tasks
and statistics, `Test_LCO_perurbations_for_local_minimum.jl` runs the scan,
and `perturbation_plots.jl` loads data and plots. Use `--help` for options.
Related table-only experiments have their own guides:
[under-observed tasks](Under_observed/README.md) and
[initialization sensitivity](Initialization_sensitivity/README.md).
