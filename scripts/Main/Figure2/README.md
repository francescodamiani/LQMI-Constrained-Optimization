# Figure 2: calibrated elbow reaching

Study a 20° elbow reach under increasing internal noise using TOD and LCO.

## Run

With Julia **1.8.5**, run from the repository root:

```bash
julia scripts/Main/Figure2/fig2.jl
julia scripts/Main/Figure2/fig2.jl --panels=G
```

The default loads `data/published/Main/Figure2/` without simulation and writes
eight PDFs to `figures/generated/Main/Figure2/`. The second command updates only G.

| Panel | Quantity | PDF |
|---|---|---|
| A | LCO elbow angle, mean ± SD | `A_elbow_angle.pdf` |
| B | Terminal elbow SD as a percentage of the 20° reach | `B_endpoint_variability.pdf` |
| C | LCO elbow velocity, mean ± SD | `C_elbow_velocity.pdf` |
| D | Analytical expected cost C | `D_expected_cost.pdf` |
| E,F | Time-averaged gain log-pseudodeterminants | `E_logK.pdf`, `F_logL.pdf` |
| G | Variance transmission divided by its zero-internal-noise baseline | `G_control_variability.pdf` |
| H | Mean per-trial peak absolute elbow velocity ± SEM | `H_peak_velocity.pdf` |

A/C use `σ_η=0.019`. B is `100 SD(θ_elbow,T)/20°`; it is not a baseline change.
E/F sum log singular values above `1e-6`, then average gain indices `2:38`.

## Model and noise settings

There are 6 states, 9 controls and 3 angle observations, with 41 stored states
at Δt=0.01 s (0.40 s). Initial angle is −20° at the elbow; initial covariances
and the target are zero. Running/terminal state costs are
`diag(1000,1e-5,1000,100,1e-5,100)` / `diag(1000,5000,1000,100,100,100)`;
control cost is `0.005 I₉`.

External-noise amplitudes are
`σ_ξ=0.01`, `σ_ω=deg2rad(0.9)`, `σ_ϵ=0.10`, `σ_ρ=0.01`.
The sweep has 20 internal-noise levels in `[0,0.04]`, 500 TOD iterations,
500 LCO cycles and 2,000 Monte Carlo trials per method/level (seed 1).
LCO starts from TOD with regularization `[0]`, fallback `[1e-10,1e-8,1e-6]`.
The process amplitude is held fixed on reruns; calibration targeted 8.5%
terminal elbow SD and gave 8.55236% at `σ_ξ=0.01`.

## Panel G definition

For each method, using unbiased across-trial centered covariance at each time:

```math
\mathcal V_{z\to u}(\sigma_\eta)=
\frac{\sum_{t=1}^{39}\operatorname{Tr}[L_t\operatorname{Cov}(z_t)L_t^\top]}
     {\sum_{t=1}^{39}\operatorname{Tr}[\operatorname{Cov}(z_t)]}.
```

G plots `𝒱_{z→u}(σ_η)/𝒱_{z→u}(0)`, normalized separately for TOD and LCO.
Both curves start at **1**, with a horizontal reference at 1. This is a ratio,
not a percentage change or the isolated contribution of internal noise.

## Rerun and replot

```bash
julia scripts/Main/Figure2/fig2.jl --simulate --tag=my_run
julia scripts/Main/Figure2/fig2.jl --data=data/generated/Main/Figure2/my_run
```

Results go to `data/generated/Main/Figure2/TAG/`; existing names require
`--force`. Reruns preserve the model and seeds but need not reproduce the saved
Monte Carlo draws exactly.

`fig2.jl` is the entry point; `reaching_task.jl` defines the model,
`fig2_simulation.jl` runs the sweep, `fig2_data.jl` handles saved results,
`fig2_statistics.jl` defines statistics, and `fig2_plots.jl` renders PDFs.
See the [dataset guide](../../../data/published/Main/Figure2/README.md) and `--help`.
