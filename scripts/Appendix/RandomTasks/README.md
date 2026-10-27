# LCO local-perturbation experiments

Related numerical experiments are organized separately:
[under-observed benchmarks](Under_observed/README.md) and
[initialization sensitivity](Initialization_sensitivity/README.md).
Their defaults print tables from archived results and create no figures.

This experiment tests changes in expected cost
around one LCO solution under relative perturbations of control gains L,
estimation gains K, and the joint vector `(vec(L), vec(K))`.

## Reproduce the plots from saved data

With Julia **1.8.5**, run from the repository root:

```bash
# Install pinned dependencies once:
julia --project=. -e 'using Pkg; Pkg.instantiate()'
# Plot all three panels without simulations or optimization:
julia scripts/Appendix/RandomTasks/perturbations.jl
```

The complete numerical archive is included under
`data/published/Appendix/RandomTasks/`. It is gzip-compressed to fit ordinary
GitHub file limits; the loader decompresses it temporarily and removes the
temporary file automatically. No manual extraction is needed.

Only three PDFs are generated, under
`figures/generated/Appendix/RandomTasks/`:

| Perturbation | PDF | Content |
|---|---|---|
| L | `Fig_LCO_local_perturbation_violin_deltaC_L_only.pdf` | Control-gain perturbations |
| K | `Fig_LCO_local_perturbation_violin_deltaC_K_only.pdf` | Estimation-gain perturbations |
| Joint | `Fig_LCO_local_perturbation_violin_deltaC_joint.pdf` | Joint perturbations |

Each plot uses red violins, gray median markers/lines, logarithmic
ΔC axis, epsilon labels and legends. Repeated plotting replaces these PDFs.
Scripts never write to `figures/published/`.

## What each script does

| Script | Purpose |
|---|---|
| [perturbations.jl](perturbations.jl) | Command-line entry point; activates the project, loads stored results by default, or explicitly runs the full experiment before plotting. |
| [perturbation_helpers.jl](perturbation_helpers.jl) | Stable random-task generation, L/K packing, normalized Gaussian directions, finite-value summary statistics, epsilon labels and analytical cost wrapper. |
| [Test_LCO_perurbations_for_local_minimum.jl](Test_LCO_perurbations_for_local_minimum.jl) | Computes the TOD initializer and LCO solution, scans the three perturbation families and saves every cost evaluation and summary. |
| [perturbation_plots.jl](perturbation_plots.jl) | Loads native or gzip-wrapped JLD2 data and renders the three panels; performs no optimization or cost recomputation. |

The latter three files are helpers loaded by the entry point, rather than
independent command-line programs. All inputs are included in the repository.

## Optional full simulation

```bash
# Recompute at full publication settings, save results, then plot:
julia scripts/Appendix/RandomTasks/perturbations.jl --simulate --tag=new_scan
# Later, plot the generated run without recomputing it:
julia scripts/Appendix/RandomTasks/perturbations.jl --data=data/generated/Appendix/RandomTasks/new_scan
```

The full scan evaluates **2.4 million analytical costs** and can be expensive;
it can take a substantial time. New results go to
`data/generated/Appendix/RandomTasks/TAG/` as a JLD2 file. These local results
are Git-ignored and never overwrite published data. Use a new tag or explicit
`--force` to replace an existing generated run. Relative `--data` paths resolve
from the repository root. `--data` cannot be combined with `--simulate`.
Use `--help` to see the available options.

The default numerical settings remain n=5, T=10, internal noise 0.2, 500 TOD
iterations and 500 LCO cycles. Each family uses 100,000 independent Gaussian
directions, normalized in its full parameter space. Both signs are evaluated
at relative amplitudes `[1e-4, 1e-3, 1e-2, 1e-1]`, reusing each direction across
amplitudes. L-only and K-only displacements scale with their own baseline
Frobenius norms; joint displacements scale with the packed vector norm.
The scan uses analytical expected cost, not Monte Carlo. The
regularization lists and disabled GD routes are preserved.

Every archived cost increment is positive. This is empirical evidence from
sampled directions, not a proof of local optimality or a Hessian test.

See the [data guide](../../../data/README.md) for archive provenance and availability.
