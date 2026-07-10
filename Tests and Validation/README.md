# Tests and Validation

All Viking HPC validation work for the custom laser-injection feature,
consolidated into one place (previously split across `Viking_results/` and
this folder — `Viking_results/` no longer exists).

- **`campaign_A_injector_2x2_validation/`** and **`campaign_B_tight_focus_f1/`**
  — the two newest (July 2026) campaigns testing the 2D/3D file-injected
  laser profile pipeline (`use_custom_profile`/`use_spatiotemporal_profile`
  on `epoch_dev`). Full narrative, debugging log, and job IDs are in
  `VIKING_PROMPT_urgent_2x2_and_f1_demo.md` in this same folder. See
  `results_summary.ipynb` for a walkthrough with commentary. Each dimension's
  `decks/`+`scripts/`+top-level files are the distilled results (as pulled
  from Viking); the per-cell run directories alongside them (e.g.
  `amp_file_phase_file/`, `lasy/`, `paraxial/`) plus `generate_profile.py` /
  `generate_lasy.py` / `generate_paraxial.py` and the campaign-level
  `physics_params.py` are the actual profile-generator scripts and injected
  `.dat` files that produced those decks, for anyone who wants to regenerate
  or adapt them (kept at their original relative depth from `physics_params.py`
  — the scripts locate it via `Path(__file__).parent.parent`, so don't move
  them further apart). Note: `scripts/` contains the *current* `analyse.py`
  (with the snapshot-timing-interpolation fix below); an older copy that used
  to live alongside the generators has been dropped as stale.
- **`2D/`** and **`3D/`** — older (June 2026) validation milestones, each
  test folder self-contained (generator + decks + results + figures
  together, unlike the `decks`/`scripts`/generator split above). See
  `2D/README.md` for that batch's own workflow notes, and `3D/HANDOVER.md`
  for the two 3D tests that were designed but postponed (not yet run).
- **`shared_libraries/`** — common Python helpers (`viking_analysis_lib*.py`,
  `analytical_field.py`, `physics_params_*.py`) used by the two newest
  campaigns' analysis scripts.
- **`lasy_epoch_export/`** — pytest suite (July 2026) for the
  `file_format="epoch"` (3D) and `"epoch2d"` (lasy-y slice) exporters
  added to the lasy fork, which write the headerless amplitude/phase
  `.dat` pair consumed by `use_spatiotemporal_profile`. Pure Python (no
  EPOCH build or Viking runs needed): verifies the on-disk Fortran
  layouts, [0, 1] normalisation, the
  `phase = −(φ − φ_ref) + π/2 + carrier_phase_ref` cos→sin convention via
  pointwise field reconstruction, seam-free phase unwrapping (the file is
  bilinearly interpolated by EPOCH), 2D/3D cross-consistency, and the
  deck-parameter metadata sidecar. See its own `README.md` for the full
  test-by-test breakdown.
- **`campaign_C_lasy_exporter_vs_native/`** — the simulation-level
  complement (July 2026, designed; Viking runs pending): the in-lasy
  exporter end-to-end in real EPOCH runs vs EPOCH's native deck Gaussian,
  recycling campaign A's beam/grids/native decks and campaign B's lasy
  pipeline. Includes a local no-EPOCH gate
  (`2d/crosscheck_legacy_pipeline.py`: new exporter vs campaign B's
  hand-rolled conversion — measured agreement ~6e-4 of peak) and
  `VIKING_PROMPT_campaign_C.md` with resource estimates. Expected
  cell-vs-cell residual is the ~1–2 % lasy-vs-paraxial physics gap, not
  campaign A's ~1e-3 % injector floor — see its `README.md` for the
  calibrated failure signatures.

## Campaign A / B methodology update (3 July 2026): snapshot-timing interpolation

## Methodology update (3 July 2026): snapshot-timing interpolation

Every `x_focus`/`w0` number below now uses a corrected method. Output
snapshots only exist every `dt_snapshot` (10fs for Campaign B, 25fs for
Campaign A), so the ideal analysis time (`t_peak@focus = t_centre +
x_spot/c`, when the pulse envelope's peak is predicted to be exactly at
the geometric focus) essentially never lands exactly on an available
snapshot. The OLD method snapped to whichever snapshot was nearest; the
NEW method linearly interpolates the FITTED w(x) curve between the two
bracketing snapshots to the exact ideal time (interpolating the raw
oscillating field itself would be unsound — the carrier completes several
cycles per dt_snapshot, so linear interpolation between distant-in-phase
snapshots wouldn't reconstruct the intermediate field correctly; the
fitted width varies slowly over the pulse duration and interpolates
cleanly). Each `summary.txt` reports both the new interpolated number and
the old nearest-snapshot one side by side. Effect size: negligible for
Campaign A (sub-nm shifts, RMS numbers unchanged at the ~1e-4% level);
modest for Campaign B 2D (paraxial-vs-LASY RMS: 10.9% → 10.0%) — real,
but not the dominant effect on its own.

## Campaign A — `campaign_A_injector_2x2_validation/`

**Question**: does file-based injection (binary amplitude/phase files)
reproduce EPOCH's own native deck-expression evaluation? Moderate NA
(~0.265), so the native-deck result is close to ground truth — isolates
injector-pipeline error, not paraxial-approximation error.

**Result**: clean validation in both 2D and 3D. In 3D, all three cells
(native, file-amplitude, file-amplitude+phase) agree to **~4.5e-4% RMS**
in beam waist w(x). See `2d/summary.txt` and `3d/summary.txt`.

2D also has a second, independent confirmation: `2d/field_comparison_summary.txt`
compares the full raw oscillating E-field (not just the beam-waist
envelope) between the native-deck and file-injected runs directly —
**0.016% max relative difference**, animated + a couple of screenshots
included (`field_comparison_t00.png`/`_t14.png`). 3D version of this same
check hasn't been run yet (not blocking — the envelope-based 3D result
above is already conclusive).

## Campaign B — `campaign_B_tight_focus_f1/`

**Question**: at f/1 (NA=0.5), does the closed-form paraxial formula still
describe the beam, or does non-paraxial (LASY angular-spectrum)
propagation differ measurably? Injector-only — no native-deck baseline.

- **`2d/`** — full-resolution, complete and trustworthy. **Key number:
  10.0% RMS difference** between paraxial and LASY beam-waist curves
  (was 10.9% before the timing-interpolation fix above), confirming the
  paraxial approximation breaks down at f/1 as expected.
- **`3d_smoke_test/`** — coarse-resolution (~2M cells vs. 147.2M full-res)
  sanity check, run before committing the full 3D allocation. Both
  generators show a genuine converging waist minimum (not the
  "diverges instead of converges" signature of an earlier sign bug) —
  pass, cleared the full-resolution 3D runs to submit. (Historical —
  predates the fixes below; not being rerun.)
- **`3d_full_resolution/`** — **complete and trustworthy** (rerun 3-4 July
  2026 with the fixes below). **Key numbers**: paraxial w0=0.679µm
  (+6.6% vs theory), LASY w0=0.643µm (+1.0%, essentially exact),
  LASY x_focus=5.095µm (+0.03%, near-perfect) — **paraxial-vs-LASY RMS =
  14.7%**, confirming the paraxial approximation breaks down at f/1 in
  3D too, consistent with the 2D result (10.0%). Some quantitative
  difference between 2D/3D is expected (different Gouy-phase/mode
  structure), same qualitative conclusion.

## Resolved: Campaign B 3D data-quality issues (3-4 July 2026)

The original 3D full-resolution numbers (w0 ~36-41% off theory) had TWO
separate bugs, now fixed and confirmed by the rerun above:

1. **Under-resolved injected file.** The `.dat` profile files (`N_TR=400`,
   52.5nm grid spacing) were generated at the SAME resolution regardless
   of simulation grid — fine at coarse smoke-test resolution (sim
   dy=131nm, file still 2.5× finer) but at full resolution the sim grid
   (dy=31.25nm) became FINER than the file (1.68× coarser), forcing
   EPOCH's boundary sampler to interpolate this f/1 beam's
   rapidly-curving wavefront from an under-sampled file. Diagnosed by
   comparing w0/theory ratios across every cell in both campaigns: only
   Campaign B's 3D cells were anomalous (ratio ~1.4, suspiciously close
   to √2), and the smoke test (same code, correctly-resolved file at
   that resolution) showed no such anomaly — ruling out a fitting-code
   bug. Fixed: regenerated both `.dat` files at `N_TR=1200` (17.5nm
   spacing, ~1.8× finer than the full-res sim grid), matching the margin
   convention already used in Campaign A's generator. **Confirmed
   fixed** — w0 errors dropped from 36-41% to 1-7% after rerunning.
2. **`analytical_field.py` bug** (used by `field_comparison.py`, not
   `analyse.py`): a helper function's return values were mislabelled,
   causing the closed-form beam width to SHRINK away from focus instead
   of growing — a 17× error at the injection boundary. Fixed (not yet
   re-run against Campaign B's data, see caveat below).

Regenerating the files bigger also surfaced a separate memory bug when
rerunning at full resolution — see the `use_pre_balance` finding further
down, which was blocking and is now also resolved.

**Side finding, relevant beyond this campaign — `use_pre_balance` vs.
large injected files.** Regenerating the `.dat` files bigger (N_TR=1200)
triggered an OOM even at `--mem 480G`: EPOCH's per-rank-slab optimisation
for large spatiotemporal injector files (`custom_laser.f90`,
`local_slab_window`) has a documented fallback to storing the FULL file
on EVERY rank whenever `use_pre_balance` (defaults to `.TRUE.`, startup
load balancing) is active — `64 ranks x 2 files x 5.76GB ~= 737GB`, more
than a node's 515GB. Fixed here by setting `use_pre_balance = F` in both
decks, safe only because these are particle-free (vacuum) runs with
nothing for pre-balance to do. **This will NOT generalise** to a future
particle simulation that needs both real load balancing (non-uniform
particle density) AND a large injected file — that combination has no
current fix short of a coarser file, more memory/nodes, or a genuine
source change in `custom_laser.f90` to compute the per-rank window after
pre-balance settles instead of needing a defensive fallback. Full
technical detail in `VIKING_PROMPT_urgent_2x2_and_f1_demo.md`.

## Resource usage (actual, not requested)

Read directly from `sacct` for each job's final/successful run. "Time %"
and "Mem %" are usage against the `--time`/default-`--mem` request
(`5.2G/core` unless noted) — i.e. how much of the requested allocation
was actually needed, not an estimate.

### EPOCH simulation jobs ("Test A"/"Test B" cells)

| Test | Job ID | ntasks | Time used | Time % | Mem used | Mem % |
|---|---|---|---|---|---|---|
| Campaign A 2D — amp_deck_phase_deck | 35476665 | 8 | 00:01:08 | 0.9% | 0.32 GB | 0.8% |
| Campaign A 2D — amp_file_phase_deck | 35476666 | 8 | 00:02:48 | 2.3% | 0.32 GB | 0.8% |
| Campaign A 2D — amp_file_phase_file | 35476667 | 8 | 00:02:08 | 1.8% | 0.32 GB | 0.8% |
| Campaign A 3D — amp_deck_phase_deck | 35476658 | 48 | 00:48:20 | 20.1% | 11.16 GB | 4.4% |
| Campaign A 3D — amp_file_phase_deck | 35477451 | 48 | 00:42:17 | 17.6% | 11.76 GB | 4.7% |
| Campaign A 3D — amp_file_phase_file | 35477452 | 48 | 00:47:42 | 19.9% | 11.68 GB | 4.6% |
| Campaign B 2D — paraxial | 35477450 | 8 | 00:02:18 | 1.9% | 0.41 GB | 1.0% |
| Campaign B 2D — lasy | 35476969 | 8 | 00:08:59 | 7.5% | 0.41 GB | 1.0% |
| Campaign B 3D smoke test — paraxial | 35476655 | 16 | 00:07:10 | 6.0% | 20.15 GB | 24.0% |
| Campaign B 3D smoke test — lasy | 35476656 | 16 | 00:09:02 | 7.5% | 20.14 GB | 24.0% |
| Campaign B 3D full-res — paraxial | 35504423 | 64 | 00:43:15 | 12.0% | 30.00 GB | 8.9% |
| Campaign B 3D full-res — lasy | 35504424 | 64 | 00:44:55 | 12.5% | 30.00 GB | 8.9% |

Every simulation finished well under its requested walltime and memory —
consistent with this project's own guidance to pad generously rather
than trim close to the estimate. The two Campaign B 3D full-res jobs are
the ones that previously OOM'd at `--mem=480G` (69GB+ used, climbing)
before the `use_pre_balance` fix; post-fix they need only 30GB, a ~16×
reduction, confirming that diagnosis directly.

### Analysis/post-processing jobscripts

| Job | Job ID | ntasks | Time used | Time % | Mem used | Mem % |
|---|---|---|---|---|---|---|
| Smoke-test waist-scan analysis (both cells) | 35477079 | 4 | 00:01:44 | 5.8% | 1.42 GB | 6.8% |
| Campaign A 2D `field_comparison.py` | 35492363 | 4 | 00:02:33 | 5.7% | 0.38 GB | 1.3% |
| Campaign B 2D `field_comparison.py` | 35492439 | 4 | 00:02:07 | 4.7% | 0.63 GB | 2.1% |
| Campaign A 2D `analyse.py` (final) | 35495345 | 2 | 00:03:33 | 17.8% | 0.29 GB | 2.7% |
| Campaign A 3D `analyse.py` (final) | 35495346 | 4 | 00:21:03 | 46.8% | 51.81 GB | 43.2% |
| Campaign B 2D `analyse.py` (final) | 35495347 | 2 | 00:03:32 | 17.7% | 0.19 GB | 1.8% |
| Campaign B 3D `analyse.py` (final, both fixes) | 35513347 | 4 | 00:20:55 | 17.4% | 92.63 GB | 51.5% |

The 3D `analyse.py` jobs are the memory-hungry ones — SDF's per-file
block storage means even loading a single field snapshot pulls the
whole multi-GB block into memory before any slicing happens (see the
execution log for the two OOM-and-refix cycles this caused during the
session: Campaign A 3D needed bumping from the login-node default up to
120G, Campaign B 3D full-res up to 180G).

## File key (each subfolder)

- `summary.txt` — plain-text write-up of the result
- `metrics.csv` — x_focus/w0 numbers per cell, machine-readable
- `waist_vs_x.png` (or `smoke_test_waist_vs_x.png`) — the beam-waist-vs-x plot
- `scripts/` — the analysis script(s) that produced everything in this
  folder (`analyse.py` = the waist-scan method above; `field_comparison.py`
  = a newer analytical-vs-numerical raw E-field comparison, see caveat
  below). Needs `sdf_xarray`/`xarray`/`scipy`/`matplotlib` and the SDF
  data still sitting on Viking to re-run — not standalone on a laptop.
- `decks/` — each cell's `input.deck`, renamed `<cell>_input.deck`

`shared_libraries/` (top level, not per-campaign) holds the common code
the scripts above import: `viking_analysis_lib.py` (2D plotting/analysis
helpers), `viking_analysis_lib_3d.py` (3D equivalent), `analytical_field.py`
(closed-form paraxial field generator, see caveat below), and each
campaign's `physics_params_*.py` (beam parameters).

## Caveat: `field_comparison.py` / `analytical_field.py`

A newer analysis (animated analytical-vs-numerical-vs-difference E-field
triptychs, normalised to a0) is included as scripts. Three bugs in
`analytical_field.py` were found and fixed: (1) boundary-vs-focus
amplitude normalisation, (2) a causality violation in the temporal
envelope, (3) a width/amplitude-prefactor conflation causing the
closed-form beam width to shrink away from focus instead of growing (17×
error at the boundary). Campaign A's version of this analysis (comparing
two EPOCH runs against each other, not a closed-form formula, so
unaffected by any of the three bugs above) works cleanly (0.016% max
relative difference) and is trustworthy.

**Campaign B has now been re-run against all three fixes (jobs 35516982
2D, 35516983 3D) and the raw-field gap did *not* close** — max
|difference| is still ~148% of peak amplitude at the pulse peak, for
*both* the paraxial cell (which should closely track the formula) and
the lasy cell. Visual inspection of the comparison frames
(`results/figures/field_comparison_paraxial_t04.png` etc.) shows why:
the analytical and numerical envelopes and spatial extent genuinely
match, but the two fields drift in and out of carrier phase along the
propagation direction (e.g. compare x≈0–2µm vs x≈4–6µm in the t=40fs
frame — the sign pattern flips). This is **not** consistent with
ordinary PIC numerical dispersion: at this deck's resolution (32
points/wavelength, 2nd-order FDTD) the expected phase-velocity error is
~0.16%, which over ~5 wavelengths of propagation only accumulates to
~0.05 rad — far too small to produce the near-half-cycle phase flip
seen here.

This points to a **4th, previously undiagnosed bug**, most likely in
`analytical_field.py`'s carrier propagation phase (`k0·ξ`) or
wavefront-curvature phase (from `R(ξ)`), since those are the terms that
would make the analytical/numerical phase relationship drift *with x*
rather than stay fixed — unlike the three already-fixed bugs, which
affected the envelope/amplitude, not the carrier phase. This is why
`analyse.py`'s waist-fitting numbers (Hilbert-envelope based, insensitive
to carrier phase) remain correct and are unaffected by this issue — the
headline Campaign B validation numbers in this README stand. Only the
raw-field `field_comparison.py` triptychs/animations for Campaign B are
affected, and should currently be treated as illustrative (envelope
shape/extent) rather than quantitative confirmation. Not yet
investigated further — flagged here per project convention rather than
chased down, since it needs a careful read of the phase-term derivation
in `analytical_field.py` rather than another rerun.
