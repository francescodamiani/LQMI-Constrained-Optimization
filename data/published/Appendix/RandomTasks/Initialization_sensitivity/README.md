# Initialization-sensitivity data

`lco_initialization_study_results.jld2` contains all 202 runs: 501-value cost
histories, final L/K gains, TOD gains, every random direction, complete task
matrices, initial conditions, settings and the result table. Task/direction
seeds are 4321/4322; runs use 500 complete L/K passes.

Tables are calculated from stored rows. All final costs are finite and below
TOD; 183/200 random starts finish within 0.1% of TOD-initialized LCO, and
199/200 within 0.2%. No saved run reaches the `1e-10` relative-change tolerance.

Run `julia scripts/Appendix/RandomTasks/Initialization_sensitivity/initialization_sensitivity.jl`
from the repository root. The
[experiment guide](../../../../../scripts/Appendix/RandomTasks/Initialization_sensitivity/README.md)
gives rerun commands and reproducibility limits. Generated results go to
`data/generated/Appendix/RandomTasks/Initialization_sensitivity/TAG/` and can
be inspected with `--data`; no figures are produced.
