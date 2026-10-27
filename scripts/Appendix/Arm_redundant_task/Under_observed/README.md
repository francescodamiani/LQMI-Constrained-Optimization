# Under-observed redundant-arm stabilization

This is the three-joint stabilization model with only elbow angle
observed. It is separate from the calibrated reaching task in main Figure 2.
Using Julia **1.8.5**, run from the repository root:

```bash
# Print the saved six-condition numerical table; no simulation:
julia scripts/Appendix/Arm_redundant_task/Under_observed/under_observed.jl
# Explicitly recompute TOD/LCO at publication settings:
julia scripts/Appendix/Arm_redundant_task/Under_observed/under_observed.jl --simulate --tag=new_arm
# Reduced numerical validation, without figures:
julia scripts/Appendix/Arm_redundant_task/Under_observed/under_observed.jl --simulate --validate --tag=check_arm
# Print an existing generated table:
julia scripts/Appendix/Arm_redundant_task/Under_observed/under_observed.jl --data=data/generated/Appendix/Arm_redundant_task/Under_observed/new_arm
```

`under_observed.jl` is the entry point; default output is Markdown on stdout.
`run_arm_redundant_task_severely_underobserved_benchmark.jl` preserves the
solve, cost, storage and diagnostics body. It reuses
[arm_task.jl](../arm_task.jl)'s plant and changes only H, D and the scalar
sensory covariance. Shared path/table helpers are in
[benchmark_tables.jl](../../benchmark_tables.jl).

Settings: m=6, p=9, k=1, T=300, Δt=0.01,
`H=[0,1,0,0,0,0]`, inertia `diag(1.2,0.8,0.5)`, damping 2I,
nine-actuator map, Euler dynamics, zero initial means/covariances,
angle/velocity weights 1/0.001, and control weight 0.01 with zero terminal
control weight. All four external noise amplitudes are 0.1 and
`σ_η=[0.05,0.10,0.15,0.20,0.30,0.40]`.

TOD uses 500 iterations and initializes LCO independently at each noise level.
LCO uses 500 passes, one L/K update per pass, `rho=[0]`, fallbacks
`[1e-10,1e-8,1e-6]`, and no internal GD routes. Costs use the validated
general analytical routines. **No GD or Monte Carlo experiment is added.**
Reduced validation uses T=12, σ_η=0.1, and three TOD/LCO passes.

The plant has observability rank **2/6**, and is retained without modification.
Archived LCO improvements are **15.0056%, 41.3492%, 61.2825%, 73.7357%,
86.2803%, 91.7578%**. All policies and full histories are included,
so all twelve analytical objectives can be independently checked.

New results go only to `data/generated/Appendix/Arm_redundant_task/Under_observed/TAG/`.
Existing tags are protected; `--force` allows resuming saved rows. Reduced and
full runs cannot be mixed. No figures or duplicate table files are produced.
Saved timings are machine-dependent and can include compilation overhead.
