# Final calibrated elbow-reaching dataset

This directory contains **all 75 fields** from the original refined-calibration
archive, with every numerical value and full trajectory preserved. No simulation
was needed to obtain these results.

Original source:
`Reviews/results/elbow_reaching_internal_noise_sweep_baseline_configuration_calibrated/elbow_reaching_internal_noise_sweep_baseline_configuration_calibrated.jld`.
Original producer: `run_elbow_reaching_internal_noise_sweep_calibrated.jl`.
The corresponding calibrated plotting script and detailed result/rebuttal text
agree with the stored noise settings and benchmarks. Older uniform-noise and
low-process-noise archives are not used.

Original JLD bytes: **547,185,214**. Original SHA-256:

```text
5a02bcfcb9c2f9f5ac3440b0304713c27da68719a85408051ee704b3b9d83a0d
```

## Lossless public layout

The original JLD/HDF5 container is converted to JLD2 using the existing pinned
public environment. JLD/HDF5 was used only in a temporary migration environment;
it is not required by public scripts. Exact equality was checked for **every
original key and value after reconstructing the split archive**.

- `sweep.jld2`: all original metadata, model matrices, covariances, gain arrays,
  complete TOD/LCO histories, individual trial costs, analytical/MC costs and
  summaries, convergence flags and process-calibration history.
- `trajectories_sigma_eta_*.jld2`: one file per internal-noise level, containing
  all six original arrays: `x_store_TOD`, `z_store_TOD`, `u_store_TOD`,
  `x_store_LCO`, `z_store_LCO`, `u_store_LCO` for that level.
- [SHA256SUMS](SHA256SUMS): SHA-256 checksums for all 21 JLD2 files. Verify from
  this directory with `shasum -a 256 -c SHA256SUMS`.

The 21 numerical files total **502,630,037 bytes**; the largest is 25,102,218
bytes. Splitting by condition keeps every file below GitHub's ordinary file
limit while preserving the entire experiment, with no LFS or external download.
The layout changes storage only; no values, trials, precision or fields are
removed. Checksums-manifest SHA-256:

```text
b562dfa18f87698d66367689d9561c509c9f8b7aeec0d18bba3651a4dc36b476
```

## Key fields and axes

All original names and array conventions are retained. Load the reconstructed
dictionary through `load_fig2` in `scripts/Main/Figure2/fig2_data.jl`.

| Fields | Contents |
|---|---|
| `x_store_*`, `z_store_*` | 20 arrays per method, each `(6,2000,41)`, axes **dimension,trial,time** |
| `u_store_*` | 20 arrays per method, each `(9,2000,40)`; nominal `u=Lz` |
| `K_*_all`, `L_*_all` | 20 full gain tensors per method, respectively `(6,3,40)` and `(9,6,40)` |
| `trial_costs_*` | All 2,000 trial totals per method and noise level |
| `cost_history_LCO` | All original 1,001-entry update histories per noise level |
| `cost_history_TOD`, `cost_history_iterations_TOD` | Source initial/final TOD cost endpoints and indices, not unsaved intermediate TOD iterations |
| `A`, `B`, `H`, `C`, `D`, `S`, `Ω_ξ`, `Ω_ω`, initial means/covariances, Q/R fields | Complete saved task definition |
| `analytical_cost_*`, `mean_MC_cost_*`, `sem_MC_cost_*` | Exact objectives and original MC estimates/SEMs |
| Endpoint, peak, gain, improvement fields | Original stored summaries, recomputed and checked by the public validation suite |
| `calibration_*`, `calibrated_sigma_xi`, target/tolerance/iteration-limit fields | Complete original process-calibration record: first candidate 0.01 accepted |
| `TOD_converged`, `LCO_converged`, `method_LCO`, `ρ_LCO`, iteration counts | Original flags and solver selections; the source flags are finite-cost/nonincrease checks, not a proof of mathematical convergence |

The exact sensory amplitude is `deg2rad(0.9)`, not rounded 0.0157. Process
amplitude is 0.01, motor amplitude 0.10 and multiplicative sensory amplitude
0.01. Horizon is 41 stored states at 0.01 s; internal amplitudes are
`[0,.002,.004,.006,.008,.011,.013,.015,.017,.019,.021,.023,.025,.027,.029,.032,.034,.036,.038,.04]`.
The original MC seed is 1 and there are 2,000 trials per method/condition.

## Use and derived statistics

```bash
julia scripts/Main/Figure2/fig2.jl
julia scripts/Main/Figure2/fig2.jl --simulate --tag=new_run
julia scripts/Main/Figure2/fig2.jl --data=data/generated/Main/Figure2/new_run
```

The default command uses only this archive. It computes plot statistics from
stored trajectories/gains and writes exactly eight PDF panels without running
optimization, Monte Carlo, PCA or any new experiment. Panel D uses saved
general analytical expected costs, validated against the current core. A/C use
LCO at `σ_η=.019` with SD bands; H uses mean trial peak absolute elbow velocity
with SEM. G uses the original unbiased within-time centered covariance and
times **1:39** to compute
`V = sum_t Tr[L_t Cov(z_t) L_t'] / sum_t Tr[Cov(z_t)]`.
It plots `𝒱(σ_η)/𝒱(0)`, separately normalized for TOD and LCO,
on a linear axis with a reference line at 1. Both curves start at 1.
Run `julia scripts/Main/Figure2/fig2.jl --panels=G` to update only that PDF.
All quantities are calculated directly from these unchanged numerical files.
