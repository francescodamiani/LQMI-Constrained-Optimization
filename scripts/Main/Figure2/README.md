# Camera-ready Figure 2: calibrated elbow reaching

This experiment replaces the main-text stabilization example with a 20-degree
elbow reach using the final psychophysically informed noise calibration. The
original stabilization experiments remain under
[`Appendix/Arm_redundant_task`](../../Appendix/Arm_redundant_task/README.md).

## Reproduce the eight panels

Run from the repository root with **Julia 1.8.5**:

```bash
# Install the pinned environment once:
julia --project=. -e 'using Pkg; Pkg.instantiate()'
# Plot the complete archived calibrated sweep; no optimization or simulation:
julia scripts/Main/Figure2/fig2.jl
# Regenerate only G from the same archive:
julia scripts/Main/Figure2/fig2.jl --panels=G
```

Exactly eight individual PDFs are written to `figures/generated/Main/Figure2/`.
There is no composite, panel lettering or smoke-test figure. Repeated runs
replace the corresponding PDFs. Scripts never write to `figures/published/`.

| Panel | PDF | Quantity and statistics |
|---|---|---|
| A | `A_elbow_angle.pdf` | LCO elbow angle, mean ± across-trial sample SD; dashed zero target |
| B | `B_endpoint_variability.pdf` | LCO terminal elbow sample SD divided by 20 degrees, displayed as percent |
| C | `C_elbow_velocity.pdf` | LCO elbow velocity, mean ± across-trial sample SD, deg/s |
| D | `D_expected_cost.pdf` | TOD/LCO analytical expected cost, logarithmic vertical axis; no sampling ribbon |
| E | `E_logK.pdf` | TOD/LCO mean thresholded log-pseudodeterminant of filter gains |
| F | `F_logL.pdf` | TOD/LCO corresponding control-gain measure |
| G | `G_control_variability.pdf` | Normalized internal-state variance transmission to control, `𝒱(σ_η)/𝒱(0)`, relative to each method's own zero-noise value; baseline 1 |
| H | `H_peak_velocity.pdf` | TOD/LCO mean of per-trial peak absolute elbow velocity ± sample SEM, deg/s |

**Caption details:** A and C both use the saved `σ_η=0.019` condition; shading
is SD, not SEM. A/B/C have no legends; captions should identify LCO and the
representative noise level for A/C. B is `100 SD(θ_elbow,T)/(20°)`, not an
increase from the zero-noise baseline, with no invented uncertainty ribbon.
H averages each trial's `max_t abs(x[5,trial,t])`; it differs from the maximum
of C's mean trajectory. Elbow angle is state component 2 and elbow velocity
component 5. Velocity is labeled `v_elbow`, cost is `C`, and the abbreviated
E/F axes read `log(|K|_p)` and `log(|L|_p)`; the time averaging remains unchanged.
All kinematics use seconds and degrees.

## What each script does

Use `fig2.jl` as the command-line entry point. The other files are included
helpers, not separate command-line programs.

| Script | Purpose |
|---|---|
| [fig2.jl](fig2.jl) | Activates the root environment, selects saved-data plotting or explicit simulation, and protects published/manual output paths. |
| [reaching_task.jl](reaching_task.jl) | Original reaching model, calibrated noise matrices and costs; calls the validated core for TOD, LCO and Monte Carlo. |
| [fig2_simulation.jl](fig2_simulation.jl) | Recomputes the 20-point sweep, initializing LCO independently from the matching TOD solution, and saves complete results. |
| [fig2_data.jl](fig2_data.jl) | Saves/loads all numerical fields using one metadata/gain file and one complete trajectory file per noise level. |
| [fig2_statistics.jl](fig2_statistics.jl) | Computes kinematics, endpoint SD, elbow peaks, gain measures, centered control variance, variance transmission and its ratio to the zero-noise baseline. |
| [fig2_plots.jl](fig2_plots.jl) | Renders only the eight agreed panels from numerical results; performs no optimization. |

## Final model and calibration

The state order is shoulder/elbow/wrist angles followed by their velocities.
There are six states, nine controls and three angle observations. The original
inertias, damping, muscle map and Euler dynamics are retained. Initial physical
and internal means are `(0,-π/9,0,0,0,0)`, with zero covariances and zero target.
The timestep is 0.01 s: **40 transitions, 41 states**, lasting exactly 0.40 s.

| Parameter | Final value |
|---|---|
| Running state cost | `diag(1000,1e-5,1000,100,1e-5,100)` |
| Terminal state cost | `diag(1000,5000,1000,100,100,100)` |
| Control cost | `0.005 I₉` for 40 controls; unused terminal slice zero |
| Process amplitude `σ_ξ` | **0.01**, recovered from both the calibrated script and saved results |
| Additive sensory amplitude | `deg2rad(0.9) = 0.015707963267948967`, retaining full precision |
| Multiplicative motor amplitude | 0.10, implemented as `C=0.10B` |
| Multiplicative sensory amplitude | 0.01, implemented as `D=0.01H` |
| Internal amplitude | 20 rounded values in `[0,0.04]`; `Ω_η=σ_η² I₆` |
| TOD / LCO | 500 iterations / 500 cycles; one K and one L update per cycle |
| Regularization | `ρ_list=[0]`, fallback `[1e-10,1e-8,1e-6]`; both GD routes disabled |
| Monte Carlo | 2,000 trials per method/level; seed 1 reset for each comparison |

The original preliminary calibration targeted 8.5% terminal elbow SD with
tolerance 0.2 percentage points. Its first trial at `σ_ξ=0.01` gave 8.55236%,
so no rescaling occurred. Reruns hold this final calibrated process amplitude
fixed; they do not refit it to a different simulator's random realization.

E/F retain the source's absolute singular-value threshold `s>1e-6`, summing
natural logarithms of retained singular values, then averaging gain indices
**2:38**. An empty retained spectrum contributes zero. Gains may be rectangular
or rank deficient; ordinary determinants are not used. No zero-noise gain
point is overwritten to force equality between methods.

G uses unbiased within-time centered covariance, with denominator `N_trials−1`,
and the matching reference's state/gain indices **1:39** (0–0.38 s):
`V_u = sum_t Tr[L_t Cov(z_t) L_t']`. First normalize by summed internal-state
variance, as in the supplied reference:

```text
𝒱(σ_η) = V_u(σ_η) / sum_t Tr[Cov(z_t)]
Plotted ratio = 𝒱(σ_η)/𝒱(0)
```

TOD and LCO each use their own `𝒱(0)`; both curves begin at **1**. G plots
the ratio against `σ_η` on a linear axis labeled `𝒱(σ_η)/𝒱(0)`,
with a horizontal reference line at 1. There is no time averaging or `Δt` factor.
The numerator equals summed sample control variance because `u_t=L_t z_t`;
this provides an independent check from archived controls. Raw sums and ratios
remain available in the statistics. This measures normalized transmission of
internal-state variance, not the isolated contribution of internal computation noise.

## Optional full rerun

```bash
julia scripts/Main/Figure2/fig2.jl --simulate --tag=new_reaching_run
julia scripts/Main/Figure2/fig2.jl --data=data/generated/Main/Figure2/new_reaching_run
julia scripts/Main/Figure2/fig2.jl --help
```

New results go to `data/generated/Main/Figure2/TAG/`, remain Git-ignored, and
cannot replace published data. Choose a new tag or explicitly use `--force`
to replace generated data for a tag. Relative `--data` paths resolve from the
repository root; `--data` and `--simulate` cannot be combined.

The new simulator uses the definitive core's native `(dimension,time,trial)`
arrays, converting saved arrays to the original `(dimension,trial,time)` schema.
Explicit RNG resets preserve common random numbers between TOD/LCO, but draw
order and RNG type differ from the historical simulator: new Monte Carlo
realizations are statistically reproducible, not bit-identical to the archive.

The saved dataset is about **479 MiB**, with all original trajectories retained.
Allow additional Julia/plotting memory and a similar disk budget for each new
full run.

See the [data guide](../../../data/README.md) for archive provenance and availability.
