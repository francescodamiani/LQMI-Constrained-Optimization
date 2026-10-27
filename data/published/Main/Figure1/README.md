# Figure 1 numerical data

This directory contains **all numerical datasets required for Figure 1b–g**.
`julia scripts/Main/Figure1/fig1.jl` loads these files directly, runs no simulations, and
writes the six PDF panels and their composite under `figures/generated/Main/Figure1/`.
`figures/published/` is reserved for complete figures added manually.

| File | Bytes | Contents | Panels |
|---|---:|---|---|
| `Fig1_ab_data_10000_GD.jld2` | 24,349 | `df_min`, all 240 rows × 5 columns | b,c |
| `Fig1d_LCO_cost_trajectories.jld2` | 77,360 | Seven complete 1,001-entry histories, summary, problem and settings | d |
| `Fig_internal_noise_sweep_full_solutions_IC_at_LCO_solutions_10000_GD.jld2` | 139,580 | Fixed problem, full gain sequences and all numerical summaries | e,f,g |

SHA-256 checksums:

```text
7fee4c6ba31a6e28f06c0f245bacf4838d3e01f80116e962b54dd161e5491808  Fig1_ab_data_10000_GD.jld2
0851f07105fa4b970276ca3a324da1568e952e5d3ac44720fbba20a7e132b07b  Fig1d_LCO_cost_trajectories.jld2
ca0fc697d918b8eb208f89cc00c6f02a7bfc25f28c989ef7140d777edb800e20  Fig_internal_noise_sweep_full_solutions_IC_at_LCO_solutions_10000_GD.jld2
```

The dimension and noise-sweep files were copied byte for byte from the original
NeurIPS supplementary submission, with no removed fields or reduced precision.
The convergence file is a byte-identical copy of the validated full-settings
run regenerated with the validated core. All three datasets share this directory and
the same plotting workflow. No script accesses the private historical source.
New simulations write exclusively to `generated/Main/Figure1/TAG/` (Git-ignored) and
never replace these saved reference files. Camera-ready provenance is not assumed.

## Convergence: d

Producer: `scripts/Main/Figure1/fig1_convergence.jl`, adapted from
`Monotonic_Convergence_LCO.jl`. The saved file contains `cost_trajectories`,
`df_summary`, `problem_data`, and `settings`.

`cost_trajectories["sigma_eta_list"]` stores [0,0.05,0.1,0.2,0.3,0.5,0.7];
`cost_trajectories["traj"]` stores all seven Float64 histories, each with
**1,001 entries**, including the initial objective and every control/filter
update. The plot shows the original x range 0.8–10; the dataset retains the
entire histories without downsampling. `df_summary` stores `sigma_eta`,
`n_iter_recorded`, `cost_init`, `cost_final`, `solve_time_sec`, `best_method`,
and `best_rho`. `cost_final` means the last recorded cost, which need not be
identical to a solver's best checkpoint in general.

The fixed task uses n=5, T=10 and problem seed 51235, with 500 TOD iterations
and 500 LCO cycles, matching the settings below. The initial and final costs
agree with saved TOD/LCO objectives within 1.5e−12. Settings record the master
seed, iteration counts and `smoke=false`. Problem matrices and initial moments
match the fixed noise-sweep problem. The original submission script saved only
a plot, so these histories were regenerated using the validated core; their
provenance is recorded here without a separate directory structure.

## Dimension benchmark: b,c

Original producer: `RandomTask_sweep_dimensions.jl`; original renderer:
`RandomTask_sweep_dimensions_Plots_ab.jl`. The renderer's "a" is manuscript c
(runtime), and its "b" is manuscript b (cost). Public producer:
`scripts/Main/Figure1/fig1_dimensions.jl`; renderer: `scripts/Main/Figure1/fig1_plots.jl`.

`df_min` has columns `n::Int`, `seed::Int`, `method::String`,
`time_sec::Float64`, `cost::Float64`. All 10 realizations per method and
dimension are retained. TOD/LCO: n = [1,2,5,10,15,20,40,60,80,100]; GD:
n = [1,2,5,10]. Horizon T=10; internal noise amplitude 0.5. Costs are analytical
trace-formula costs. Panel b uses the mean and sample-SD/√10 SEM over all rows.
Panel c preserves the original exclusion of n=1,seed=1 **for all three methods**
(9 observations there, 10 elsewhere). Runtime SEM is computed but was not drawn
in the source, so no runtime ribbons are added. This display exclusion does not
remove any archived row. LCO/GD times exclude the preceding TOD initialization;
original JIT overhead can be included in first-call timings.

This historical file contains no gains, iteration histories, raw Monte Carlo
trials, or intermediate optimizer states; the full original file is preserved.
New simulations additionally retain gain sequences, LCO histories/selection
diagnostics, GD iteration/convergence diagnostics, and run settings.

## Fixed noise sweep: e,f,g

Original producer: `RandomTask_sweep_internal_noise.jl`; original renderer:
`RandomTask_sweep_internal_noise_Plots.jl`. Public producer:
`scripts/Main/Figure1/fig1_noise.jl`. The archived sweep has n=5, seed index 1, T=10,
σ_eta = [0,0.05,0.1,0.2,0.3,0.5,0.7], and 1,000 MC trials per method/noise.

Top-level keys: `df_sweep`, `solutions_all`, `problem_data`, `sigma_eta_list`,
`n_fixed`, `seed_fixed`, `N_trials_MC`.

`df_sweep` (21 × 9) contains `n`, `seed`, `sigma_eta`, `method`, `solve_time_sec`,
`mc_cost_mean`, `mc_cost_sem`, `raw_objective`, `N_trials_MC`. Panel e uses the
stored MC means and SEMs, without resampling or analytical substitution.

`solutions_all` contains all nine `L_METHOD`, `K_METHOD`, `M_METHOD` Float64
arrays for METHOD = TOD,LCO,GD, each with shape **(5,5,9,7)**:
matrix rows, matrix columns, time 1:T−1, noise level. M_t=A+B L_t−K_t H.
It also contains `sigma_eta_list` and seven-entry vectors for each method's
`raw_objective`, `mc_cost_mean`, `mc_cost_sem`, and `solve_time_sec`.
Panels f,g compute each slice's Frobenius distance to the corresponding GD gain,
then the arithmetic mean of those norms over **all nine time slices**, including
the final filter slice. There is no averaging over random systems or SEM here.
That last filter slice affects z_T only and is unidentifiable from this objective
(no terminal control or internal-state penalty). It is retained for fidelity;
filter-gain distances should be interpreted with this qualification.

`problem_data` preserves A,B,H,C,D, Q/R arrays (5,5,10), Ω_ξ and Ω_ω, initial
means and covariances, x0,z0, n_fixed,seed_fixed,Tloc, ρ_max, the four fixed
noise amplitudes and random_seed_master. Ω_η = σ_eta² I is reconstructed from
each saved sweep value. This archive contains neither per-trial trajectories
nor per-trial costs nor LCO optimization histories. New sweeps additionally
save per-trial costs, MC seeds, LCO histories/selection diagnostics, GD
iteration/convergence diagnostics, and smoke/settings metadata.

## Model, parameters, initialization, and algorithms

All tasks are square: state/control/observation dimensions n; H=I.
Independent Gaussian A and B are individually scaled to spectral radius 0.9.
D=0.2 abs(A_D), where A_D is a separately drawn stable matrix; C=0.2 B.
Additive process/sensory covariance is 0.2² I; multiplicative motor/sensory
channels have independent scalar unit-variance draws, with amplitudes already
in C,D. Internal covariance is σ_eta² I. All five noise sources are retained.
Initial means x1=z1=ones(n), both initial covariances and cross covariance zero.
Q_t=R_t=I except Q_1=0 and R_T=0; terminal cost is x_T'Q_T x_T.
The controller convention is **u_t=L_t z_t** throughout; TOD already returns
this sign. There is no additional sign flip.

The master seed is 1234; each problem uses
`MersenneTwister(1234 + 10000*n + seed)` with seed indices 1:10 for b,c and
1 for d–g. The source MC evaluation seed is
`1234 + 400000 + round(Int,1e6*σ_eta)`, reset for each method at a noise level.
New simulations use fresh MersenneTwisters with those documented seeds.
The new simulator samples the same mathematical model with a different draw
order: exact original trial realizations are neither required nor claimed.
Arrays use `(dimension,time,trial)`; reported empirical cost and sample SEM
come directly from the core simulator. Dense copies of additive covariances
are passed to MC to support singular zero-noise covariance on Julia 1.8.5;
all covariance values are unchanged.

TOD: 500 iterations. LCO starts from TOD, calls the unchanged core FullProcedure
with ρ_list=[0], ρ_list_large=[1e−10,1e−8,1e−6], 500 cycles, iters_K=iters_L=1,
no analytical GD refinement/from-IC (both flags false, iteration counts 0).
Retained inactive GD settings: Adam, learning rate 1e−3, maxiter 50, tol 1e−6,
β1=0.9, β2=0.999, epsilon=1e−8. The core controls regularization, pseudoinverses,
update order, admissibility, and best checkpoints.
Numerical GD uses Optim.GradientDescent(), finite-difference gradients and
Optim defaults, iterations=10000, store_trace=show_trace=false. GD starts from
TOD in b,c and from LCO in e–g. Timing covers optimization only.

The supplied sweep source was edited after its seven-level archive was saved:
it selects only σ_eta=0.2, leaves `iterations_GD` unfinished and names an
`ONLY_ONE_NOISE` output. Publication settings above follow the archived grid,
the `_10000_GD` filename, the benchmark source and manuscript Table 1. The
archive does not itself store the optimizer iteration limit, so that limit's
provenance is documented rather than asserted as recovered metadata.
