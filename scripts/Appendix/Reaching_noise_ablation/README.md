# Calibrated reaching: fixed-gain noise ablation

Remove noise sources while holding the Figure 2 LCO gains fixed to quantify
their effects on endpoint and joint variability.

## Run

With Julia **1.8.5**, run from the repository root:

```bash
julia scripts/Appendix/Reaching_noise_ablation/reaching_noise_ablation.jl
julia scripts/Appendix/Reaching_noise_ablation/reaching_noise_ablation.jl --simulate --tag=my_run
julia scripts/Appendix/Reaching_noise_ablation/reaching_noise_ablation.jl --data=data/generated/Appendix/Reaching_noise_ablation/my_run
```

The default prints tables from `data/published/Appendix/Reaching_noise_ablation/`.
`--simulate` recomputes all 35 cases by analytical moment propagation only:
no optimization, Monte Carlo or figures. It reuses `data/published/Main/Figure2/sweep.jld2`
and saves one result file to `data/generated/Appendix/Reaching_noise_ablation/TAG/`.

## Settings and interpretation

The [Figure 2 model](../../Main/Figure2/README.md) uses external amplitudes
`σ_ξ=0.01`, `σ_ω=deg2rad(0.9)`, `σ_ϵ=0.10`, `σ_ρ=0.01`.
Five selected internal amplitudes are `[0,.011,.019,.029,.040]`, nearest the
nominal targets `[0,.01,.02,.03,.04]`; gains are not interpolated.

Seven conditions comprise the full model, removal of process, additive sensory,
motor multiplicative, sensory multiplicative or internal noise individually,
and combined removal of both sensory channels. K/L stay fixed throughout.

Terminal elbow SD is `100*sqrt(Var(x[2,T]))/(π/9)`. Mean angle/velocity SD
averages instantaneous SD across three joints and times `2:41`, in rad and
rad/s. Each reduction is `100*(1-SD_ablated/SD_full)`; combined sensory effects
are evaluated together, and reductions are not additive variance shares.

For nonzero internal noise, terminal-SD reductions are 65.3551–78.2432% for
process noise and 1.1623–5.1861% for internal noise. These ranges exclude η=0.
Analytical full-model endpoint SD is 8.63464–10.54364%; Figure 2 uses Monte Carlo
sample SD instead, explaining its small numerical difference.

`reaching_noise_ablation.jl` is the entry point; `ablation_tables.jl` reads and
prints results, and `fixed_gain_ablation.jl` recomputes them. See the
[dataset guide](../../../data/published/Appendix/Reaching_noise_ablation/README.md)
and `--help` for options.
