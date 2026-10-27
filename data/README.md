# Numerical datasets

`published/` contains the complete available numerical records used by the
paper. These files are immutable inputs to reproduction commands. `generated/`
is initially empty and Git-ignored; explicit reruns write only there.
No public script requires historical/private source material.

| Published folder | Contents and provenance |
|---|---|
| `Main/Figure1/` | Three complete JLD2 archives: dimension benchmark (240 rows), seven full 1,001-value convergence histories, and the fixed-task noise sweep. The convergence archive is the validated full-settings reconstruction; the other two are original submission files. See the [dataset schema](published/Main/Figure1/README.md) and [Figure 1 guide](../scripts/Main/Figure1/README.md). |
| [Main/Figure2](published/Main/Figure2/README.md) | All 75 fields of the final calibrated reaching archive, split losslessly into metadata/gains and 20 trajectory files. |
| `Appendix/RandomTasks/` | Original Appendix Figure 3 perturbation scan: 2.4 million costs, compressed byte-preservingly as `.jld2.gz`. The loader decompresses temporarily; allow about 307 MB temporary space. |
| [RandomTasks/Under_observed](published/Appendix/RandomTasks/Under_observed/README.md) | Complete random-task costs, timings and LCO histories; historical final gain arrays were not saved. |
| [RandomTasks/Initialization_sensitivity](published/Appendix/RandomTasks/Initialization_sensitivity/README.md) | Seed-4321 study: all 202 starts/final policies, directions, histories and diagnostics. |
| [Arm_redundant_task](published/Appendix/Arm_redundant_task/README.md) | Availability notes for historical arm figures; their original full trajectory caches were not supplied. |
| [Arm_redundant_task/Under_observed](published/Appendix/Arm_redundant_task/Under_observed/README.md) | Complete elbow-observed stabilization policies, objectives and histories. |
| [Arm_redundant_task/Noise_ablation](published/Appendix/Arm_redundant_task/Noise_ablation/README.md) | Four losslessly converted original archives, all 82 fields: policies, component variability, fixed-policy ablations and alternative external-noise regimes. |
| [Reaching_identifiability](published/Appendix/Reaching_identifiability/README.md) | Five complete legacy archives (229 logical fields), including shared raw reference trials split by noise level. This uses the earlier .01-amplitude reaching configuration. |
| [Reaching_noise_ablation](published/Appendix/Reaching_noise_ablation/README.md) | Original 35-row calibrated reaching ablation CSV; reuses the existing main Figure 2 gains without duplicating its large archive. |

Dataset READMEs retain hashes, schemas, units and source identification.
Lossless conversions were checked field by field; compression/splitting does
not remove numerical precision or trials. Individual distributed numerical
files remain below 100 MiB. Figure 2 and identifiability together account for
most of the roughly 1 GB download.

## Plot or inspect a generated run

Every plotting entry point accepts `--data=FOLDER`, using its existing archive
schema. For example:

```bash
julia scripts/Main/Figure2/fig2.jl --data=data/generated/Main/Figure2/my_run
julia scripts/Appendix/RandomTasks/perturbations.jl --data=data/generated/Appendix/RandomTasks/my_run
```

Plots always go to `figures/generated/`. Table-only commands print generated
results with `--data` and do not create plots. See each experiment README for
its explicit simulation command, cost and scientific qualifications.

The arm figure entry point can save new caches with `--save-data`; its default
runs simulations because original figure caches are absent. A full arm cache
is about 1 GB and should remain local. No missing historical trial, gain or
history is fabricated from rounded paper text. Saved Monte Carlo samples and
new seeded simulations need not be bit-identical because sampler draw order
has changed; deterministic archive evaluation remains available.
