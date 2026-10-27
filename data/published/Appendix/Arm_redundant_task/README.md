# Redundant-arm data

Full trajectory datasets for the k=3 arm figures are not bundled.
The following analyses have complete numerical archives:

- [Under_observed/](Under_observed/README.md): k=1 TOD/LCO benchmark.
- [Noise_ablation/](Noise_ablation/README.md): k=3 full-model gains,
  behavioral variance summaries, fixed-policy ablations and alternative noise regime.

Their entry points print numerical tables from saved data without figures.
The existing k=3 figure workflow retains its simulation/plot entry point:

```bash
julia scripts/Appendix/Arm_redundant_task/arm_redundant_task.jl
```

See the [arm experiment guide](../../../../scripts/Appendix/Arm_redundant_task/README.md) for parameters, cache schemas and commands.
