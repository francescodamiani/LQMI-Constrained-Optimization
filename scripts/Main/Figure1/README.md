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

LCO starts from TOD. GD uses 10,000 iterations, starting from TOD for b/c
and LCO for e–g. Regularization is `[0]`, with fallback `[1e-10,1e-8,1e-6]`.
Tasks use H=I, independent A/B scaled to spectral radius 0.9, C=0.2B and
D=0.2 times the absolute entries of a separately drawn stable matrix. Additive
process/sensory amplitudes are 0.2; initial means are ones and covariances zero.
Q/R are identity except Q₁=0 and R_T=0. Problem seeds are `1234+10000*n+seed`
(seed 51235 for d–g).

## Saved data

`load_fig1` in `fig1_plots.jl` loads the required JLD2 files from the data folder:

| File | Contents | Panels |
|---|---|---|
| `Fig1_ab_data_10000_GD.jld2` | `df_min`: 240 rows of costs, runtimes, dimensions, seeds and methods | b,c |
| `Fig1d_LCO_cost_trajectories.jld2` | Full histories, summaries, task matrices and settings | d |
| `Fig_internal_noise_sweep_full_solutions_IC_at_LCO_solutions_10000_GD.jld2` | Fixed task, gains, objectives and Monte Carlo summaries | e,f,g |

Fixed-task gain arrays have shape `(5,5,9,7)`: row, column, time, noise level.
Gain distances average all nine time slices. The dimension archive has no
final gains; the noise sweep has no per-trial trajectories or costs and does
not store the GD iteration limit. Panel d displays only the early portion
of the full saved histories.

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
