# Figure 2 calibrated reaching dataset

The complete dataset contains 75 fields, including model matrices, gains,
analytical and Monte Carlo costs, solver histories and all trial trajectories.

- `sweep.jld2`: model, settings, gains, trial costs, summaries and diagnostics.
- `trajectories_sigma_eta_*.jld2`: full x/z/u trajectories for TOD and LCO,
  one file per internal-noise level.

There are 20 noise levels and 2,000 trials per method/level. State arrays have
shape `(6,2000,41)`, controls `(9,2000,40)`, with axes `(dimension,trial,time)`.
K/L tensors have shapes `(6,3,40)` and `(9,6,40)`. The dataset is about 479 MiB.
Load all parts through `load_fig2` in `scripts/Main/Figure2/fig2_data.jl`.

## Use

Run `julia scripts/Main/Figure2/fig2.jl` from the repository root. It writes
eight PDFs to `figures/generated/Main/Figure2/` without simulations.
The [experiment guide](../../../../scripts/Main/Figure2/README.md) gives rerun
commands and the `--data` option for `data/generated/Main/Figure2/TAG/`.

## Settings and panel G

External amplitudes are `σ_ξ=0.01`, `σ_ω=deg2rad(0.9)`, `σ_ϵ=0.10`, `σ_ρ=0.01`;
internal amplitudes span `[0,0.04]`. Δt=0.01 s and T=41. Saved convergence
flags report finite-cost/nonincrease checks, not a proof of convergence.

For each method, panel G uses unbiased across-trial centered covariance over
`1:39` to compute
`𝒱_{z→u}=sum_t Tr[L_t Cov(z_t) L_tᵀ]/sum_t Tr[Cov(z_t)]`.
It plots `𝒱_{z→u}(σ_η)/𝒱_{z→u}(0)`, normalized separately for TOD and LCO.
The baseline is **1**, not a percentage change. No data are resampled.
