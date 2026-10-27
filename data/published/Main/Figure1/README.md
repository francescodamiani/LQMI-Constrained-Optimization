# Figure 1 datasets

All numerical data for panels b–g are included here.

| File | Contents | Panels |
|---|---|---|
| `Fig1_ab_data_10000_GD.jld2` | `df_min`: 240 rows of dimension, seed, method, runtime and analytical cost | b,c |
| `Fig1d_LCO_cost_trajectories.jld2` | Seven full 1,001-value histories, summaries, task matrices and solver settings | d |
| `Fig_internal_noise_sweep_full_solutions_IC_at_LCO_solutions_10000_GD.jld2` | Fixed task, TOD/LCO/GD gains, objectives and Monte Carlo summaries | e,f,g |

Run `julia scripts/Main/Figure1/fig1.jl` from the repository root to plot these
files without simulations. Outputs go to `figures/generated/Main/Figure1/`.
See the [experiment guide](../../../../scripts/Main/Figure1/README.md) to rerun
experiments and plot results from `data/generated/Main/Figure1/TAG/`.

## Numerical conventions

Tasks are square with H=I, independently drawn A/B scaled to spectral radius
0.9, C=0.2B, and D=0.2 times the absolute entries of a separately drawn stable
matrix. Additive process/sensory amplitudes are 0.2. Initial means are ones,
initial covariances zero, and Q/R are identity except Q₁=0 and R_T=0.
The master seed is 1234; problem seeds use `1234+10000*n+seed`.

The fixed task has n=5, T=10, seed 51235 and internal amplitudes
`[0,.05,.1,.2,.3,.5,.7]`. Gain arrays have axes `(row,column,time,noise)` and
shape `(5,5,9,7)`. Panel e uses saved 1,000-trial means/SEMs; f/g average
Frobenius distances to GD across all nine gain slices. GD starts from TOD
for b/c and LCO for e–g. The noise-sweep file does not store the GD iteration
limit; reruns use 10,000 iterations.

Dimension-benchmark data contain no final gains; noise-sweep data contain
no per-trial trajectories or costs. Convergence histories retain every gain
update, although panel d displays only their early portion.
