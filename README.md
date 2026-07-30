# AELP (Arbitrary EPOCH Laser Profile)

Tutorials, validation tests, and utilities for a custom laser injection
feature added to the [EPOCH](https://github.com/epochpic/epoch) PIC code:
loading a laser's amplitude and (optionally) phase profile from an external
raw binary data file, instead of specifying it as an analytic deck
expression. This lets a laser profile come from an external beam
propagation tool (e.g. [LASY](https://github.com/LASY-org/lasy)) or from
any shape that isn't one of EPOCH's built-in analytic primitives.

The EPOCH-side modification lives on the
[`upstream-pr-custom-laser-injection` branch of `epoch_dev`](https://github.com/ZheyuanChen/epoch_dev/tree/upstream-pr-custom-laser-injection);
this repository is where it is tested, validated, and documented for end
users.

## Status

Implemented and validated in both **epoch2d** and **epoch3d** (epoch1d is
unmodified). The proposed upstream changes remain on the
`upstream-pr-custom-laser-injection` branch rather than `main`. For the full,
current, session-by-session history see the
[Issues](https://github.com/ZheyuanChen/AELP/issues) tab of this repo,
which is where day-to-day development notes and bug reports now live
rather than in this README.

## Getting the modified EPOCH

Clone [`epoch_dev`](https://github.com/ZheyuanChen/epoch_dev) and check out
the `upstream-pr-custom-laser-injection` branch (the modification is **not**
on `main`). Build it in the usual way
(`make COMPILER=gfortran` from `epoch2d/` or `epoch3d/`). Once the pending
changes are available upstream, this fork-specific step will no longer be
necessary.

## Installation

If you already have a working Python environment, you can go straight to
the tutorial notebooks. Otherwise, from the top of this repository:

```
uv pip install -e .
```

(or `pip install -e .` if you don't use `uv`).

## Documentation

[DOCUMENTATION_LASER_INJECTION.pdf](DOCUMENTATION_LASER_INJECTION.pdf) is
the full reference: deck elements for both epoch2d and epoch3d, the binary
file format and axis ordering, the phase sign convention (including the LASY
conversion), two-channel polarisation, and current limitations. Read that
before writing a profile generator — the array-ordering convention in
particular is easy to get backwards silently.

## Tests and Validation

Start with the
[Campaign A/B/C results notebook](Tests%20and%20Validation/results_summary.ipynb)
and the [validation index](Tests%20and%20Validation/README.md).

- **Campaign A — injector isolation:** file amplitude/phase injection agrees
  with native deck evaluation to about $4.5\times10^{-4}\%$ RMS in fitted
  beam waist in both 2D and 3D. The independent 2D raw-field comparison is
  bounded by 0.016% maximum relative difference (0.0011% of the reference
  peak).
- **Campaign B — $f/1$ LASY demonstration:** paraxial-vs-LASY waist curves
  differ by 10.0% RMS in 2D and 14.7% in 3D. In 3D, the LASY focus is within
  0.03% of target and its waist within 1%; the 2D LASY result is illustrative
  because an axisymmetric LASY field is propagated in EPOCH2D's slab
  geometry.
- **Campaign C — LASY exporter:** the reusable LASY-to-EPOCH exporter passes
  end-to-end validation. Against Campaign B's hand-written export, propagated
  fields agree within 0.47% of peak in 2D and 0.28% on the 3D mid-plane. The
  [Campaign C report](Tests%20and%20Validation/campaign_C_lasy_exporter_vs_native/README.md)
  separates this exporter-isolating result from the moderate-NA
  LASY-vs-paraxial physics comparison.

These are vacuum field-propagation tests. They validate the profile delivery
and propagation path, not plasma coupling, particle dynamics, or QED
observables.

## Repository layout

- `tutorial/` — worked examples and Jupyter notebooks comparing the
  custom-profile injection against EPOCH's analytic laser profiles.
- `Tests and Validation/` — HPC validation campaigns (LASY-vs-paraxial
  comparisons, memory/load-balancing tests, etc.) and their analysis
  notebooks.
- `src/` — shared Python utilities (currently an `sdf-xarray` helper for
  loading EPOCH output in Jupyter).

## LASY integration

Several of the tests under `tutorial/` and `Tests and Validation/` drive EPOCH
directly from a [LASY](https://github.com/LASY-org/lasy)-generated beam,
including the amplitude/phase extraction and the LASY-to-EPOCH phase
conversion documented in
[the laser-injection reference](DOCUMENTATION_LASER_INJECTION.pdf). The
reusable exporter lives in the LASY fork as
`Laser.write_to_file(file_format="epoch" / "epoch2d")`; AELP retains its
[format/convention tests](Tests%20and%20Validation/lasy_epoch_export/README.md)
and [simulation-level validation](Tests%20and%20Validation/campaign_C_lasy_exporter_vs_native/README.md).
Some older tutorial scripts still show the pre-exporter, hand-written
conversion and are retained as historical worked examples.

## Known issues

Tracked in this repo's [Issues](https://github.com/ZheyuanChen/AELP/issues),
not here.
