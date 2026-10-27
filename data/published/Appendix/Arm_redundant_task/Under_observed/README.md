# Under-observed arm data

`arm_redundant_task_severely_underobserved_results.jld2` contains six noise
conditions, all twelve TOD/LCO policies, six full 1,001-value LCO histories,
model/cost/noise matrices, initial conditions and solver diagnostics.
The task has T=300 and one elbow-angle observation; all analytical objectives
can be reevaluated from the saved gains.

Run `julia scripts/Appendix/Arm_redundant_task/Under_observed/under_observed.jl`
from the repository root to print tables. The
[experiment guide](../../../../../scripts/Appendix/Arm_redundant_task/Under_observed/README.md)
gives settings and rerun commands. Generated runs are saved under
`data/generated/Appendix/Arm_redundant_task/Under_observed/TAG/` and read with
`--data`; no figures are produced.
