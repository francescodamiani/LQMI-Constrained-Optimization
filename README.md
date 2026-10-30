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
and data/output locations.

## Reproduce results

From the repository root, plot both main figures from saved data:

```bash
julia scripts/Main/Figure1/fig1.jl
julia scripts/Main/Figure2/fig2.jl
```

These commands run no simulations and write PDFs to `figures/generated/Main/`.
See the [Figure 1](scripts/Main/Figure1/README.md) and
[Figure 2](scripts/Main/Figure2/README.md) guides for panel selections.
Figure 2 studies a calibrated elbow-reaching task; the appendix stabilization
task uses the same redundant three-joint arm to maintain a reference posture
under stochastic disturbances.

Appendix commands are listed in their experiment READMEs. Most load saved data;
the [appendix stabilization figure workflow](scripts/Appendix/Arm_redundant_task/README.md)
runs full simulations because its figure datasets are not included.

## License

The code is licensed under the [MIT License](LICENSE).
