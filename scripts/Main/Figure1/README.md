# Figure 1: random LQMI tasks

This folder reproduces numerical Figure 1 panels **b–g**. All required numerical
data, including the complete convergence histories for panel d, are included in
`data/published/Main/Figure1/`. Panel a is manually prepared artwork.

## Reproduce the plots

Run these commands from the repository root with **Julia 1.8.5**:

```bash
# Generate every numerical panel from the included data:
julia scripts/Main/Figure1/fig1.jl
```

The command runs **no simulations or optimizations**. It writes
`fig1b.pdf` through `fig1g.pdf` and the combined `fig1_b-g.pdf` (b,c,d / e,f,g)
under `figures/generated/Main/Figure1/`. Only PDFs are generated; repeated runs
replace these derived files. Nothing is written to `figures/published/`.

## What each script does

Use `fig1.jl` as the command-line entry point. The other files are helpers loaded
by the entry point, rather than separate commands.

| Script | Purpose |
|---|---|
| [fig1.jl](fig1.jl) | Activates the project, parses options, loads saved data by default, and dispatches plotting or explicitly requested simulations. |
| [fig1_plots.jl](fig1_plots.jl) | Loads the datasets, calculates plot statistics and gain distances, and renders panels b–g and their composite. Runs no optimization. |
| [fig1_tasks.jl](fig1_tasks.jl) | Defines the random square task construction and gain packing/matrix helpers used by the optional numerical experiments. |
| [fig1_dimensions.jl](fig1_dimensions.jl) | Recomputes the dimension benchmark for b–c: TOD, LCO and numerical GD costs and solver runtimes over 10 systems per dimension. Saves the benchmark dataset. |
| [fig1_convergence.jl](fig1_convergence.jl) | Recomputes d: the fixed five-dimensional task at seven internal-noise levels, with TOD initialization and 500 LCO cycles. Saves all 1,001 cost values per history, plus settings and diagnostics. |
| [fig1_noise.jl](fig1_noise.jl) | Recomputes e–g: the fixed task at seven internal-noise levels, TOD/LCO/GD gains, analytical costs and 1,000-trial Monte Carlo summaries. Saves the noise-sweep dataset. |

## Panel and data mapping

Dataset names below are relative to `data/published/Main/Figure1/`.

| Panels | Plotted quantity | Dataset |
|---|---|---|
| b | Expected cost versus dimension; mean ± SEM over systems | `Fig1_ab_data_10000_GD.jld2` |
| c | Mean solver runtime versus dimension | `Fig1_ab_data_10000_GD.jld2` |
| d | LCO cost after each alternating control/filter update | `Fig1d_LCO_cost_trajectories.jld2` |
| e | Monte Carlo mean cost ± SEM versus internal noise | `Fig_internal_noise_sweep_full_solutions_IC_at_LCO_solutions_10000_GD.jld2` |
| f | Time-averaged Frobenius distance of filter gains K to GD | Same noise-sweep file |
| g | Time-averaged Frobenius distance of control gains L to GD | Same noise-sweep file |

Panel c excludes the dimension-1, seed-1 timing from all methods' plotted
means; the saved data are complete. See the
[data README](../../../data/README.md) for schemas, numerical settings and seeds.

## Optional experiment reruns

```bash
# Full publication settings, including expensive GD optimizations:
julia scripts/Main/Figure1/fig1.jl --simulate --tag=new_run
# Only panel d's full convergence experiment; no numerical GD:
julia scripts/Main/Figure1/fig1.jl --simulate --panels=d --tag=new_convergence
# Plot a previously generated full-settings run:
julia scripts/Main/Figure1/fig1.jl --data=data/generated/Main/Figure1/new_run
# Plot only the saved convergence panel:
julia scripts/Main/Figure1/fig1.jl --panels=d
# Reduced numerical validation; saves data but generates no figures:
julia scripts/Main/Figure1/fig1.jl --simulate --smoke --tag=test_run
```

Full GD runs can take days. New numerical data are saved under
`data/generated/Main/Figure1/TAG/`, which is Git-ignored; published datasets
are protected. Use a new tag, or `--force` to replace an existing generated run.
The tag names numerical data, while PDF filenames remain fixed. For a partial
simulation cache, select its corresponding panels when plotting with `--data`.
Reduced smoke datasets cannot be used to produce publication plots.

`--panels=b,c,d,e,f,g` selects panels and the simulations needed for them.
Relative `--data` paths resolve from the repository root. Run the entry point
with `--help` for its options. Scripts also work from another working directory
when invoked by absolute path.

See the [data guide](../../../data/README.md) for archive provenance and availability.
