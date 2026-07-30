# Handover — postponed 3D validation tests

> **Status:** historical design hand-over. Campaigns A/B/C and the current
> validation conclusions are summarised in
> [`../results_summary.ipynb`](../results_summary.ipynb). The two tests in
> this folder remain postponed.

Built on 2 July 2026, but **postponed** in favour of the urgent
[`Campaign A`](../campaign_A_injector_2x2_validation/) and
[`Campaign B`](../campaign_B_tight_focus_f1/) runs. Their execution history
is retained in [`../README.md`](../README.md). The tests below have not been
submitted to Viking; everything is designed and, where noted, was
smoke-tested locally at low resolution.

## test1_converging_gaussian/ — moderate-NA 3D field-level benchmark

**Status: fully built, smoke-tested locally (low-res, converges correctly
to the right focus/waist), NOT yet run at full Viking resolution.**

What it tests: quantitative validation of the file-injection pipeline
against closed-form Gaussian-beam theory (waist position, waist size,
transverse-power energy conservation), at moderate NA (~0.265, the same beam
as `../campaign_A_injector_2x2_validation/`) — the "does the whole pipeline
reproduce known physics" check, one level more demanding than the
file-vs-native comparison in Campaign A (which only checks internal
consistency, not absolute correctness against theory).

Files: `physics_params.py` (single source of truth), `generate_profile.py`
(writes the static (y,z) amplitude+phase planes), `numerical/input.deck`
(file-injected) + `analytical/input.deck` (native deck expression, for the
3-way numerical/analytical/theory decomposition), `analyse.py` (waist
scan, energy-conservation ratio corrected for the pulse's own temporal
envelope shape — see `common/viking_analysis_lib_3d.expected_envelope_
shape`'s docstring for why raw P(x) is not flat for a pulsed beam).

Resource estimate: 434x451x451 = 88.3M cells, vacuum, `t_end`=350fs ->
**`--ntasks 48 --time 4:00:00`** (same sizing logic as Campaign A's 3D
cells, since it's the identical box/beam).

To run: `python generate_profile.py` (writes into `numerical/`), then run
both decks, then `python analyse.py`.

## test2_two_channel_polarisation/ — 3D two-channel polarisation independence

**Status: built (generator + deck), NOT smoke-tested successfully yet** —
the first attempt hit the per-rank-slab MPI abort (see below), which is
now fixed. Re-attempt the smoke test before trusting this at Viking scale.

What it tests: repeats, at 3D/Viking scale, the exact recipe that caught a
real per-laser-storage bug in epoch2d (two `x_min` lasers, `pol_angle` 0
and pi/2, same-grid-size-but-different-content spatiotemporal files,
deliberately avoiding the grid-mismatch abort to isolate the more
dangerous silent-wrong-answer failure mode). The epoch3d port was built
per-laser from the start (code-reading inference), but this test measures
it rather than just trusting that inference — matching this project's
general practice of independent verification.

Files: `generate_profiles.py`, `run/input.deck` (two laser blocks, ids 1
and 2), `analyse.py` (fits Ey against laser A's own input AND laser B's
input, and Ez against both, to make a storage-independence regression
structurally loud — if Ez matches laser A instead of laser B, that's the
exact old bug signature).

**Directly relevant history from this session**: the first local smoke
test of this exact deck (4 MPI ranks, spatiotemporal path, file grid
spanning the full box) hit `ERROR: custom laser profile sampled outside
the stored per-rank slab` — a real bug in epoch3d commit `9d50a6a6`'s
per-rank memory optimisation, root-caused and fixed same session in
commit `aa605cca` (see the epoch3d-per-rank-slab-mpi-bug memory / the
"Known MPI limitation" section that was added then updated in
`epoch_dev/TIGHT_FOCUSING_INJECTION_PROMPT.md`). This test's own deck is
exactly the kind of two-laser, full-box-spanning spatiotemporal case that
triggered it, so it's worth re-running the smoke test now that the fix is
in before trusting a full Viking submission — should just work, but
hasn't been re-verified.

Resource estimate: 480x480x480 = 110.6M cells, vacuum, `t_end`=30fs ->
per the sizing table (3D vacuum 50-200M) **`--ntasks 48 --time 2:00:00`**
(short run, small `t_end`, generous padding).

## Open items

- Campaign B's analysis is now complete and retained under
  `../campaign_B_tight_focus_f1/`; no analysis implementation is outstanding
  there.
- `dev_test/laser_profile_injection/tight_focusing_f1/` (30 June, 2D,
  array-only, text format) is now stale against current epoch2d (binary-
  only) — not updated this session; superseded in spirit by
  Campaign B but the old scripts/results remain useful historical reference
  material.
- The broader local epoch3d custom-laser regression suite has since advanced
  beyond the partial state described in the original hand-over. Consult the
  current `epoch_dev` and `Project_EPOCH` validation reports before using
  this dated checklist as an execution authority.
