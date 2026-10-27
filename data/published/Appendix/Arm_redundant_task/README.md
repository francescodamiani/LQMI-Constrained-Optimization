# Historical redundant-arm data

The complete original k=3 arm-figure trajectory caches are unavailable; no
replacement figure datasets are generated or bundled. Separate rebuttal
experiments have complete numerical archives:

- [Under_observed/](Under_observed/README.md): k=1 TOD/LCO benchmark.
- [Noise_ablation/](Noise_ablation/README.md): original k=3 full-model gains,
  behavioral variance summaries, fixed-policy ablations and alternative noise regime.

Their entry points print numerical tables from saved data without figures.
The existing k=3 figure workflow retains its simulation/plot entry point:

```bash
julia scripts/Appendix/Arm_redundant_task/arm_redundant_task.jl
```

See the [arm experiment guide](../../../../scripts/Appendix/Arm_redundant_task/README.md) for parameters, cache schemas and commands.
