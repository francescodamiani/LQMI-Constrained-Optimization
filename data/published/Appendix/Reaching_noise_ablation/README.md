# Calibrated reaching noise-ablation data

`calibrated_fixed_gain_noise_ablation_subset.csv` contains all 35 rows:
seven fixed-gain conditions at internal amplitudes `[0,.011,.019,.029,.040]`.
It reuses the calibrated Figure 2 gains in `data/published/Main/Figure2/sweep.jld2`.

Columns contain `sigma_eta`, `condition`, terminal elbow SD (% of the 20°
reach), mean angle SD (rad), mean angular-velocity SD (rad/s), and the percentage
reduction of each SD relative to the full model. Angle/velocity SDs average
instantaneous SD across three joints and times `2:41`.

All LCO gains stay fixed. Combined sensory removal zeros both sensory
channels together; reductions are not additive variance shares. Analytical
endpoint SDs differ slightly from the Monte Carlo SDs in Figure 2.

Run `julia scripts/Appendix/Reaching_noise_ablation/reaching_noise_ablation.jl`
from the repository root to print tables. The
[experiment guide](../../../../scripts/Appendix/Reaching_noise_ablation/README.md)
gives settings and analytical rerun commands. Generated results go to
`data/generated/Appendix/Reaching_noise_ablation/TAG/calibrated_reaching_ablation.jld2`
and can be read with `--data` pointing to the run directory. No figures,
optimization or Monte Carlo runs are produced.
