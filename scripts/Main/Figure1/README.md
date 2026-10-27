# Figure 1: random LQMI tasks

Compare TOD, LCO and gradient descent (GD) across task dimensions and internal
noise, and inspect LCO convergence. Panel a is manually prepared artwork;
this workflow produces panels b–g.

## Run

With Julia **1.8.5**, run from the repository root:

```bash
julia scripts/Main/Figure1/fig1.jl
julia scripts/Main/Figure1/fig1.jl --panels=d
```

The default loads `data/published/Main/Figure1/` without simulation and writes
`fig1b.pdf`–`fig1g.pdf` plus `fig1_b-g.pdf` to `figures/generated/Main/Figure1/`.

## Settings and panels

- **b,c:** analytical cost and runtime across 10 tasks per dimension. TOD/LCO
  cover dimensions 1–100; GD covers 1, 2, 5 and 10. Cost bars are mean ± SEM.
  Runtime means exclude dimension 1, seed 1 for every method; saved rows remain complete.
- **d:** seven full LCO cost histories for a five-dimensional task, T=10,
  internal noise `[0,.05,.1,.2,.3,.5,.7]`, with 500 TOD iterations and 500 LCO
  cycles (1,001 recorded costs per history).
- **e–g:** the same fixed task and noise grid; 1,000 Monte Carlo trials per
  method/level, followed by time-averaged Frobenius gain distances to GD.
  The final filter slice is retained although its value is not determined by this objective.

LCO starts from TOD; numerical GD uses 10,000 iterations. Regularization is
`[0]`, with fallback `[1e-10,1e-8,1e-6]`. Details and dataset contents are in the
[data README](../../../data/published/Main/Figure1/README.md).

## Rerun and replot

```bash
julia scripts/Main/Figure1/fig1.jl --simulate --panels=d --tag=my_run
julia scripts/Main/Figure1/fig1.jl --data=data/generated/Main/Figure1/my_run --panels=d
```

Omit `--panels=d` to recompute all panels, including expensive GD runs that can
take days. Results go to `data/generated/Main/Figure1/TAG/`. Select the same
panels when replotting a partial run. Existing run names require `--force`;
`--smoke` is available for reduced numerical checks, not paper plots.

`fig1.jl` is the entry point; `fig1_tasks.jl` defines tasks,
`fig1_dimensions.jl`, `fig1_convergence.jl` and `fig1_noise.jl` run experiments,
and `fig1_plots.jl` computes plot statistics and renders PDFs. Use `--help`
for all options.
