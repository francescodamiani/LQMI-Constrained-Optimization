# LQMI constrained optimization

A Julia library for linear-quadratic stochastic optimal control with additive,
multiplicative and internal computation noise, accompanying **“Solving Stochastic
Control under Multiplicative and Internal Noise via Constrained Optimization”
(NeurIPS 2026)**.

LCO alternates updates of the control and estimation gains to reduce expected
cost. The main solver is `Optimal_Control_Estimation_FullProcedure()` in
[src/MyFunctionsOptimalRegularizedSolutions.jl](src/MyFunctionsOptimalRegularizedSolutions.jl).
It uses pseudoinverses with regularization to handle numerical difficulties,
including singular or poorly conditioned matrices. Its arguments and
implementation are documented in the source code.

The solver can also serve as an inner solver in inverse optimal control
frameworks that require repeated stochastic optimal control solves.

## Install

Use **Julia 1.8.5**. From the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Project.toml` and `Manifest.toml` pin the environment. To load the library:

```julia
include("src/MyFunctionsOptimalRegularizedSolutions.jl")
```

## Repository guide

```text
src/                 Core library and shared helpers
scripts/Main/        Figure 1 and Figure 2 experiments
scripts/Appendix/    Supplementary experiments and numerical tables
data/published/      Numerical datasets accompanying the paper
data/generated/      Results produced by users (Git-ignored)
figures/published/   Final paper figures
figures/generated/   PDF plots produced by scripts (Git-ignored)
```

**Every experiment folder contains its own README** with commands, settings,
and data/output locations. The [data guide](data/README.md) links all datasets.

## Reproduce results

From the repository root, plot both main figures from saved data:

```bash
julia scripts/Main/Figure1/fig1.jl
julia scripts/Main/Figure2/fig2.jl
```

These commands run no simulations and write PDFs to `figures/generated/Main/`.
See the [Figure 1](scripts/Main/Figure1/README.md) and
[Figure 2](scripts/Main/Figure2/README.md) guides for panel selections.
Appendix commands are listed in their experiment READMEs. Most load saved data;
the [redundant-arm figure workflow](scripts/Appendix/Arm_redundant_task/README.md)
runs full simulations because its figure datasets are not bundled.

## Run simulations and plot the results

For example, rerun the Figure 1 convergence experiment and plot its saved result:

```bash
julia scripts/Main/Figure1/fig1.jl --simulate --panels=d --tag=my_run
julia scripts/Main/Figure1/fig1.jl --data=data/generated/Main/Figure1/my_run --panels=d
```

Figure 2 and perturbation experiments use the same `--simulate` / `--data`
pattern. Arm figures use `--save-data` to retain simulation results. Table-only
experiments use `--data` to inspect a saved run. Full GD and Monte Carlo sweeps
can be expensive; consult the experiment README and `--help` before running.
Published data and paper figures are never overwritten.

Released under the [MIT License](LICENSE).
