# Reaching: internal noise, sensory-model mismatch and alternative noise sources

These experiments ask whether alternative mechanisms reproduce the
joint endpoint-variability, peak-velocity and gain signatures of internal noise
within the tested reaching model and parameter ranges. They do **not** establish
general identifiability of internal computation noise from behavior alone.

Run with Julia **1.8.5**, using the repository's pinned environment:

```bash
# Saved-data tables only: no optimization, simulation, figures or output files.
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl
# All raw A conditions, all 100 B rows, and both complete C sweeps:
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl --full
# One analysis, or a comma-separated selection:
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl --experiment=underestimation
# Explicit full-settings rerun; potentially expensive:
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl --simulate --tag=new_run
# Small numerical validation only, clearly marked as reduced:
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl --simulate --validate --tag=check
# Print saved generated results, without rerunning:
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl --data=data/generated/Appendix/Reaching_identifiability/check
```

`--experiment` accepts `overestimation`, `noise_sources`, `underestimation`, or
`all`. Published data are read from
`data/published/Appendix/Reaching_identifiability/`; reruns write only to
`data/generated/Appendix/Reaching_identifiability/TAG/`. Existing tags require
an explicit `--force` and matching checkpoint grids. Simulation loops support checkpointing and resuming. Default tables show 12 significant digits;
JLD2 arrays retain full available precision. No command creates figures.

| Script | Purpose |
|---|---|
| `reaching_identifiability.jl` | Entry point, experiment selection, archive loading and generated-directory protection. Simulations load only with `--simulate`. |
| `identifiability_statistics.jl` | Numerical tables, overlap-only interpolation for A, gain norms/log volumes, achieved-vs-target SD errors for B, own-baseline cost changes for A/C. |
| `identifiability_simulation.jl` | Loads the validated core once; adapts the core Monte Carlo array axes, writes JLD2 checkpoints and dispatches the three experiments. Reuses the Figure 2 storage helpers and matching dynamics/cost constructor for reduced validation. |
| `sensory_overestimation.jl` | Sensory-model mismatch sweep, including its underestimation branch and saved internal-reference comparisons. Optimizes under assumed noise and evaluates frozen gains in the true system. |
| `matched_noise_sources.jl` | Four independently reoptimized external-noise candidate grids, interpolation/fallback, reoptimization at matched amplitudes and extended Monte Carlo statistics. |
| `sensory_underestimation.jl` | Dense raw mismatch and correct-model internal-noise sweeps at true sensory amplitude 0.05; includes Frobenius norms and log-pseudodeterminants. |

## Common model and numerical settings

The plant, cost matrices, initial conditions and noise covariances for these
experiments are stored in `reference/sweep.jld2`.
State/control/observation dimensions are **6/9/3**, Δt=0.01 s, duration=0.40 s,
T=41 stored samples and 40 gain slices. State order is joint angles 1–3,
then angular velocities 1–3. Initial state/estimate means are
`[0,-π/9,0,0,0,0]`, with zero covariances and target zero.

The arm has inertia diag(1.2,0.8,0.5), damping 2I, and the same sparse 3×9
actuator map as main Figure 2. H observes angles. Running Q is
diag(1000,1e-5,1000,100,1e-5,100), terminal Q is
diag(1000,5000,1000,100,100,100), and R=0.005I₉ except its unused terminal
slice is zero. There is no time-step multiplier in the objective.

Baseline amplitudes ξ,ω,ε,ρ are all **0.01**, with
Ωξ=σξ²I₆, Ωω=σω²I₃, Ωη=ση²I₆, C=σεB and D=σρH.
These are the baseline noise amplitudes for this experiment.

Each optimized condition uses **500 TOD iterations**, then **500 LCO cycles**
initialized by TOD, `iters_K=iters_L=1`, regularization list `[0.0]` and fallback
`[1e-10,1e-8,1e-6]`. Both GD routes are disabled. Inactive Adam settings
are retained: η=1e-3, maxiter=50, tolerance=1e-6, β₁=.9, β₂=.999, ε=1e-8.
The Monte Carlo seed is **1**, reset per condition. Simulation reruns use the core Monte Carlo routine. Their trial draw order
can differ from that used for the saved trials, so identical stochastic
realizations are not guaranteed.

## A: assumed sensory-noise sweep, true ω=0.01

η is zero. All **21 log-spaced assumed amplitudes from 0.001 to 0.1** are
retained; the central value is exactly 0.01. Default tables show the baseline
and overestimation branch; `--full` also shows the lower branch.
The correct baseline reuses the archived η=0 LCO policy and trajectories.
Every other condition first solves TOD/LCO with the assumed Ωω, freezes K/L,
then evaluates analytical cost, terminal SD and **2,000-trial** behavior using
true Ωω=0.01²I₃. No true-model reoptimization occurs during evaluation.

The internal-noise reference has 20 amplitudes:
`round.(range(0,stop=.04,length=20),digits=3)`.
Saved internal gains, trial data and costs are reused. A interpolates internal statistics against analytical terminal SD
**only inside the overlap**, with no extrapolation. It does not reoptimize an
interpolated internal policy. Nine of 21 wrong-model conditions lie inside this
overlap; the largest two overestimation conditions do not. All remain available.

## B: independent optimization of other noise sources

This experiment tests all four external-noise sources. Each comparison varies
one source, holds the others at 0.01 and sets η=0. Candidate grids are:

- ξ: 31 points, `range(0,.1,length=31)`, including .01.
- ω: `range(.0001,.1,length=31)` plus .01, sorted and deduplicated: **32 points**.
- ε and ρ: 31 points each, `range(0,.3,length=31)`, including .01.

Each candidate has independently optimized TOD/LCO gains. For each of the 20
internal-reference endpoint-SD targets, adjacent candidate crossings linearly
interpolate the **noise amplitude**, then TOD/LCO are optimized again at that
amplitude and **2,000 trials** are evaluated. A baseline amplitude of .01 reuses
the identical η=0 reference. Multiple crossings use the minimum
endpoint-distance score; an unbracketed target uses the closest grid point and
is explicitly not matched. The 1% achieved-SD warning is retained.

**Archive limitations:** the ε and ρ grids each fail to bracket **13/20**
targets; **12/20** achieved conditions in each have SD errors above 1%, reaching
14.5444% and 14.9403%. Their highest-target amplitudes saturate at .3.
The additive sources do match the full range accurately. Additive sensory
peak velocity is nearly flat with small MC reversals; the archived values do
not support strict monotonicity of every alternative curve. No grids are
extended and no missing matched conditions are invented.

## C: dense underestimation comparison, true ω=0.05

Wrong-model branch: `range(.05,.016,length=10)`, η=0. Optimize with assumed
Ωω, freeze gains, then evaluate using true Ωω=.05²I₃. Correct-model branch:
η=`0:.0025:.03` (13 points), with Ωω=.05²I₃ during optimization and evaluation.
All other noise amplitudes remain .01. The zero-η solution is computed for
this true model and reused at the correct wrong-model baseline.
These are **raw sweeps without endpoint matching or interpolation**.

Both the saved dataset and simulation use **1,000 MC trials**.

## Statistical definitions

- Endpoint variability is centered terminal elbow SD, not RMSE or bias. Its
  analytical value propagates raw moments and subtracts the
  squared mean; MC uses unbiased sample SD. Degrees convert to percent of the
  20° reach as `100*SD_deg/20`.
- Peak elbow velocity is `mean_trial(max_t(abs(rad2deg(x[5,trial,t]))))`, over
  all 41 samples; uncertainty is sample SD/√N. Peak time uses the first argmax
  and the saved time vector. This is not the peak of the mean trajectory.
- A/B gain magnitudes are mean Frobenius norms over **2:38**. B's additional
  source CSV log volumes are reproducible from saved gains using singular
  values >1e-6. C uses singular values >**1e-4**. For each slice, sum their
  natural logs; then average over **2:38**. Tables never imply entrywise gain
  changes. C has local gain-metric reversals despite its overall directions.
- A/C cost increases are `100*(C/C_correct_baseline-1)` within each experiment;
  C's two baselines coincide. Assumed-model costs remain distinct archived
  fields. Costs use the appropriate evaluation covariances with fixed gains.
- B retains extended statistics: component SD is sqrt(mean unbiased temporal
  variances), over 1:41 for x/z and 1:40 for u. PCA centers z separately at each
  time, pools trials over **2:39**, uses default 99% retention,
  and divides the first retained eigenvalue by the retained-eigenvalue sum.
  PC1 readout is mean ‖LₜPC1‖² over 2:39, with SEM over those times. The
  secondary joint-velocity-norm peak is in rad/s; the primary elbow peak is
  in deg/s. Success is |terminal elbow angle|<2°, bias is |mean endpoint|,
  and control effort is mean trial sum of squared controls without Δt.

`--validate` uses 8 stored samples, 64 trials, 3 TOD/LCO iterations, three
reference/assumed levels and nominal three-point candidate grids with .01
inserted. It checks code paths, not the full behavioral trends.
