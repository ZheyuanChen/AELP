# Campaign C — lasy-fork exporter vs EPOCH's native laser injector

End-to-end validation of the **in-lasy EPOCH exporter**
(`Laser.write_to_file(file_format="epoch" / "epoch2d")` on the
[ZheyuanChen/lasy fork](https://github.com/ZheyuanChen/lasy)) in real EPOCH
runs, against EPOCH's own native deck-expression Gaussian. This is the
simulation-level complement to the pure-Python pytest suite in
`../lasy_epoch_export/`.

## Relation to campaigns A and B (what is recycled, what is new)

- **Campaign A** proved the *file-injection path* is faithful: file vs deck
  with the **identical closed-form source** on both sides → residual
  ~1e-3 % (the injector floor). Campaign C reuses campaign A's beam, grids
  and native-cell decks verbatim (moderate NA ≈ 0.265, so the paraxial
  native deck is close to ground truth).
- **Campaign B** validated the *lasy-generated* profiles in EPOCH, but via
  a hand-rolled conversion script (`generate_lasy.py`: rt propagate +
  interpolation + skimage unwrap + manual `.tofile()`).
- **Campaign C** replaces that hand-rolled step with the exporter now built
  into the lasy fork, and compares against the native injector. The new
  things under test: the exporter's normalisation, phase referencing
  (`phi_ref` piston removal), seam-free unwrapping, `carrier_phase_ref`
  pin, Fortran file layout, and the 2D slice path.

## Design (per dimensionality: one pair of runs)

| Cell | Boundary source |
|---|---|
| `native/` | EPOCH deck expressions: `gauss(...)` amplitude, quadratic phase, `phase_const = π/2 + ψ_bnd` (half Gouy in 2D, full in 3D) |
| `lasy_file/` | Binary pair from `write_to_file`, spatiotemporal path; generator passes `carrier_phase_ref = ψ_bnd` so both cells share the same carrier-envelope reference |

Both cells use `amp = 3.2e12` — the file amplitude is normalised to 1, so
the deck carries the physical scale identically in both.

## Expected outcome (pass/fail calibration)

The two sources **deliberately differ in physics** (angular-spectrum lasy
envelope vs paraxial closed form), so the expected cell-vs-cell residual is
**~1–2 % of peak, smooth and amplitude-dominated** — the June 2025
lasy-vs-analytical numbers (~1.9 % amplitude, ~0.5 % phase), *not*
campaign A's ~1e-3 % injector floor. Failure signatures that would indicate
an exporter bug instead:

- **~200 % pointwise field residual** → sign-flipped carrier (the −π/2
  carrier-offset bug the pytest suite guards against);
- **spatially uniform few-tens-of-% residual** → CEP/carrier-pin mismatch
  (`carrier_phase_ref` vs the native deck's `phase_const`);
- **localised streaks in the residual** → phase-unwrap seams hitting
  EPOCH's bilinear interpolation.

## Local gates (run before spending Viking hours)

1. `pytest "../lasy_epoch_export/test_lasy_epoch_export.py"` — 16 tests,
   format/convention level.
2. `python 2d/generate_lasy_exporter.py` then
   `python 2d/crosscheck_legacy_pipeline.py` — checks the new exporter
   against campaign B's validated hand-rolled pipeline on the same beam.
   **Measured 10 July 2026: amplitude max diff 5.8e-4 of peak,
   amplitude-weighted field error max 5.9e-4 — PASS.**

## Files

- `physics_params.py` — imports campaign A's constants (single source) and
  adds the file-grid parameters.
- `2d/`, `3d/` — per-dimensionality: `generate_lasy_exporter.py`,
  `native/input.deck`, `lasy_file/input.deck`, `scripts/analyse.py`
  (waist-vs-x, snapshot-timing-interpolated), and in 2D
  `scripts/field_comparison.py` (pointwise E_y) plus
  `crosscheck_legacy_pipeline.py` (local, no EPOCH).
- The 2D `.dat` pair (~2.4 MB each) is committed; the 3D pair (~91 MB
  each) is **not** — regenerate with `3d/generate_lasy_exporter.py`
  (locally and rsync, or on Viking via sbatch after installing the fork).
- `VIKING_PROMPT_campaign_C.md` — run instructions + resource estimates.
