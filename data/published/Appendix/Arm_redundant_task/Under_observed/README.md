# Under-observed arm results

`arm_redundant_task_severely_underobserved_results.jld2` contains the complete
under-observed arm benchmark. The simulation script is
`scripts/Appendix/Arm_redundant_task/Under_observed/run_arm_redundant_task_severely_underobserved_benchmark.jl`.
The complete file is about 2 MiB; verify it with `SHA256SUMS`.

All top-level fields are retained, including the six-row numerical table,
all twelve TOD/LCO policies, six full 1,001-value LCO histories, plant/cost/noise
matrices, initial conditions, solver settings and observability diagnostics.
This is the T=300 stabilization model with k=1 and elbow-angle
observation; it is separate from both the k=3 arm experiments and main Figure 2.

```bash
julia scripts/Appendix/Arm_redundant_task/Under_observed/under_observed.jl
```
