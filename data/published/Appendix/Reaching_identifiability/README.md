# Reaching-identifiability numerical data

These five **JLD** archives are stored as JLD2 for the pinned Julia 1.8.5
environment. The complete dictionaries
passed fieldwise `isequal` after conversion and reload: **229 logical fields**
in total. No numerical values, precision, gains, trials or histories were
removed. Serialization bytes differ because the source format requires
JLD/HDF5, which are not dependencies of the public project.

The table identifies the source containers for each dataset. Public scripts
load only the JLD2 files included in this directory.

| Public file/directory | Source container | Fields |
|---|---|---:|
| `sensory_overestimation.jld2` | `elbow_reaching_wrong_sensory_model_sweep/wrong_sensory_model_sweep_results.jld` | 47 |
| `matched_noise_sources.jld2` | `elbow_reaching_matched_all_noise_sources_31_point_grid/matched_all_noise_sources_extended_results.jld` | 44 |
| `candidate_noise_grids.jld2` | `elbow_reaching_matched_all_noise_sources_31_point_grid/dense_grid_31_LCO_solutions_all_noise_sources.jld` | 24 |
| `sensory_underestimation.jld2` | `elbow_reaching_wrong_sensory_underestimation_dense_raw_sweeps/wrong_sensory_underestimation_results.jld` | 53 |
| `reference/` | `elbow_reaching_internal_noise_sweep_baseline_configuration/elbow_reaching_internal_noise_sweep_baseline_configuration.jld` | 61 |

Source-container SHA-256 values, in the same order:

```text
ce8d5f90a45bdebeee33bacc3b768fed954f175e4afc7bc4796efbc8e78eb6da
818036bd37a859fa48c40af8a294948a27db91496319867a20dc42bbcfc19bdb
348fbb042199529fe4c93e8cee33abada02affad5773d42d9f53b0dde1acf616
7cd5976d93a2c270424d35289a6b49451684374232ecfc81059ff15e9e5c99f8
7e9320d99dd8f980fa5ebcd8f18a3bca83ed451c8d3fcc5a548ada2be5661138
```

[`SHA256SUMS`](SHA256SUMS) records all **25 public numerical files**. Their
combined size is **508,512,508 bytes** (about 485 MiB); the largest file is
**25,099,879 bytes**, below 100 MiB. No LFS or external download is required.

The shared reference is stored once. `reference/sweep.jld2` retains all 55
non-trajectory fields, including the model, seed, trial costs,
TOD/LCO gains, solver histories and diagnostics. Twenty compressed
`reference/trajectories_sigma_eta_*.jld2` files each retain all six x/z/u
TOD/LCO arrays for that noise level. Their axes remain the
**(dimension, trial, time)**. Loading all parts reconstructs the complete
61-field dictionary; no trajectory downsampling occurs. The unchanged generic
Figure 2 storage helper is reused solely for this split format.

The shared reference uses **0.01 for all external-noise amplitudes**,
20 internal-noise levels and 2,000 trials per level. A/B reuse this reference.
The experiment-specific noise settings are documented in the
[experiment guide](../../../../scripts/Appendix/Reaching_identifiability/README.md).

A retains all 21 assumed-noise conditions, including the lower-noise branch,
true/assumed costs, wrong-model policies, analytical/MC endpoint SDs,
per-trial-peak means/SEMs, timing, biases and Frobenius gain metrics. Its
matching statistics are derived by overlap-only interpolation;
no invented matched policies are stored.

B retains all 125 dense candidate policies, costs, SDs and completion flags,
plus all 100 internal/external comparison rows, 80 reoptimized external
policies, matched amplitudes, uncertainty estimates, profiles, PCA directions,
readout/component measures, control effort and completion flags. The sensory
grid has 32 points after inserting .01. The multiplicative grids do not
bracket 13/20 endpoint targets each; 12/20 achieved comparisons in each exceed
1% SD error. Saturated nearest-grid comparisons must not be described as
successfully matched. Gain log volumes are recomputed from saved policies with the 1e-6 threshold by `--full`.

C is the **dense** 10-assumed-level / 13-internal-level archive with
true sensory amplitude .05 and **1,000** MC trials. It retains all policies,
true and assumed costs, endpoint and peak statistics, mean profiles, norm and
log-volume measures, threshold and completion flags.

The A/B/C summary archives did not save wrong/matched/C per-trial trajectories
or condition-level optimization histories; these cannot be recovered from
rounded text or CSV exports. All trial-level trajectories and histories that
were available in the shared internal-reference archive are retained. Missing
raw data are not fabricated. CSV/text exports duplicate stored or
derivable summaries and are not copied; no plots are included.

```bash
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl --full
julia scripts/Appendix/Reaching_identifiability/reaching_identifiability.jl --simulate --tag=new_run
```

Default tables read only the small summary archives and reference metadata;
large trajectories are needed only for explicit A/B reruns or trial-level
validation. Generated runs go exclusively under
`data/generated/Appendix/Reaching_identifiability/TAG/` and never overwrite
these archives. Full runs record reference metadata without duplicating its
trajectories. Reduced runs generate their own marked small reference.
