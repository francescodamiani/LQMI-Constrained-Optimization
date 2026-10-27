# LQMI constrained optimization

Code and numerical data for **“Solving Stochastic Control under Multiplicative
and Internal Noise via Constrained Optimization” (NeurIPS 2026)**.

The Julia library solves finite-horizon linear stochastic control and
estimation problems with additive, multiplicative and internal computation
noise. It implements analytical expected costs, the TOD comparison, alternating
LCO with regularized inverses, optional gradient-based optimization, and
Monte Carlo evaluation. Controls follow `u_t=L_t z_t`.

## Install

Use **Julia 1.8.5**. From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Project.toml` and `Manifest.toml` pin the scientific/plotting environment.
The core itself needs only Julia standard libraries:

```julia
include("src/MyFunctionsOptimalRegularizedSolutions.jl")
```

Its docstrings describe the APIs, dimensions, initialization convention and
numerical options. `T` counts stored states; gain arrays have `T-1` slices.

## Repository layout

```text
src/                         # Core library and shared path/output helpers
scripts/
  Main/{Figure1,Figure2}/    # Main-text figures
  Appendix/                 # Supplementary figures and numerical tables
data/
  published/                # Complete available paper datasets; never overwritten
  generated/                # Initially empty; optional simulation results
figures/
  published/                # Complete final manuscript PDFs
  generated/                # Initially empty; PDFs produced by scripts
```

**Each experiment folder under `scripts/Main/` and `scripts/Appendix/` has its
own README**, with its parameters, script roles, commands and numerical details.
The [data guide](data/README.md) links dataset schemas and provenance.

## Reproduce figures and tables

| Experiment guide | Command (after `julia `) | Default behavior |
|---|---|---|
| [Figure 1](scripts/Main/Figure1/README.md) | `scripts/Main/Figure1/fig1.jl` | Six panels b–g and their composite from saved data |
| [Figure 2](scripts/Main/Figure2/README.md) | `scripts/Main/Figure2/fig2.jl` | Eight calibrated reaching panels from saved data |
| [Random-task perturbations](scripts/Appendix/RandomTasks/README.md) | `scripts/Appendix/RandomTasks/perturbations.jl` | Three Appendix Figure 3 panels from saved data |
| [Historical arm figures](scripts/Appendix/Arm_redundant_task/README.md) | `scripts/Appendix/Arm_redundant_task/arm_redundant_task.jl` | **Runs full simulations:** original figure caches were unavailable |
| [Under-observed random tasks](scripts/Appendix/RandomTasks/Under_observed/README.md) | `scripts/Appendix/RandomTasks/Under_observed/under_observed.jl` | Saved numerical tables |
| [Under-observed arm](scripts/Appendix/Arm_redundant_task/Under_observed/README.md) | `scripts/Appendix/Arm_redundant_task/Under_observed/under_observed.jl` | Saved numerical tables |
| [Initialization sensitivity](scripts/Appendix/RandomTasks/Initialization_sensitivity/README.md) | `scripts/Appendix/RandomTasks/Initialization_sensitivity/initialization_sensitivity.jl` | Saved numerical tables |
| [Stabilization noise ablation](scripts/Appendix/Arm_redundant_task/Noise_ablation/README.md) | `scripts/Appendix/Arm_redundant_task/Noise_ablation/noise_ablation.jl` | Saved numerical tables |
| [Reaching identifiability](scripts/Appendix/Reaching_identifiability/README.md) | `scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl` | Saved numerical tables |
| [Calibrated reaching ablation](scripts/Appendix/Reaching_noise_ablation/README.md) | `scripts/Appendix/Reaching_noise_ablation/reaching_noise_ablation.jl` | Saved numerical tables |

Numerical-only experiments do not create figures. Plotting commands write
**PDF only**, under the corresponding `figures/generated/{Main,Appendix}/`
folder. Complete manuscript PDFs in `figures/published/` are separate manual
compositions and are never overwritten.

## Rerun and plot your own data

For example, recompute Figure 1 convergence and then replot it:

```bash
julia scripts/Main/Figure1/fig1.jl --simulate --panels=d --tag=my_run
julia scripts/Main/Figure1/fig1.jl --data=data/generated/Main/Figure1/my_run --panels=d
```

Figure 2 and perturbation entry points support the same `--simulate --tag=NAME`
and `--data=PATH` pattern. Arm figures use `--save-data --tag=NAME` to cache a
simulation and `--data=PATH` to replot it. Table-only experiments use `--data`
to print a generated run. Use each entry point's `--help` for selections.

Generated data and figures are Git-ignored. Existing numerical run names are
protected unless `--force` is explicit. Full GD benchmarks can take days;
arm experiments and large Monte Carlo sweeps are also substantial computations.
Reduced validation datasets are not substitutes for published plots.

Run commands above from the repository root, or invoke a script by its absolute
path from any directory. Relative `--data` paths resolve from the repository
root; absolute paths are also accepted. No experiment requires private files.
