# LCO initialization-sensitivity results

`lco_initialization_study_results.jld2` contains the complete initialization
study: task seed 4321, direction seed 4322, 50 directions, scales 1/10,
tolerance 1e-10 and 500 passes. The simulation script is
`scripts/Appendix/RandomTasks/Initialization_sensitivity/run_lco_initialization_study.jl`.
The file is about 1.7 MiB; `SHA256SUMS` records its checksum.

All fields are retained: 202 full 501-value cost histories, 202 final L/K
gain pairs, TOD gains, joint TOD vector, every direction, complete task
matrices/initial conditions, metadata and the 202-row result table. Tables are calculated directly from these complete numerical data.

```bash
julia scripts/Appendix/RandomTasks/Initialization_sensitivity/initialization_sensitivity.jl
```

Tables are recalculated from stored rows, with no solver or plotting code
loaded. There are **183/200** random finals within 0.1% of TOD-initialized LCO,
and **199/200** within 0.2%. The largest absolute gap is 0.2193836338%. All 202 finals are
finite and below TOD; no archived run reached its 1e-10 stopping tolerance.
