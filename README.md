# AELP (Arbitrary EPOCH Laser Profile)

Tutorials, validation tests, and utilities for a custom laser injection
feature added to the [EPOCH](https://github.com/epochpic/epoch) PIC code:
loading a laser's amplitude and (optionally) phase profile from an external
raw binary data file, instead of specifying it as an analytic deck
expression. This lets a laser profile come from an external beam
propagation tool (e.g. [LASY](https://github.com/LASY-org/lasy)) or from
any shape that isn't one of EPOCH's built-in analytic primitives.

The EPOCH-side modification lives in a separate repository,
[epoch_dev]([https://github.com/ZheyuanChen/epoch_dev](https://github.com/ZheyuanChen/epoch_dev/tree/upstream-pr-custom-laser-injection)); this repository is
where it gets tested, validated, and documented for end users.

## Status

Implemented and validated in both **epoch2d** and **epoch3d** (epoch1d is
unmodified). A pull request to merge this upstream into EPOCH is in
preparation. For the full, current, session-by-session history see the
[Issues](https://github.com/ZheyuanChen/AELP/issues) tab of this repo,
which is where day-to-day development notes and bug reports now live
rather than in this README.

## Getting the modified EPOCH

Clone [epoch_dev]([https://github.com/ZheyuanChen/epoch_dev](https://github.com/ZheyuanChen/epoch_dev/tree/upstream-pr-custom-laser-injection)) and check out
the `upstream-pr-custom-laser-injection` branch (this is the branch sent to official EPOCH repo for a PR — the
modification is **not** on `main`). Build it in the usual way
(`make COMPILER=gfortran` from `epoch2d/` or `epoch3d/`). Once the pending
pull request is merged upstream, this step will no longer be necessary.

## Installation

If you already have a working Python environment, you can go straight to
the tutorial notebooks. Otherwise, from the top of this repository:

```
uv pip install -e .
```

(or `pip install -e .` if you don't use `uv`).

## Documentation

[DOCUMENTATION_LASER_INJECTION.pdf](DOCUMENTATION_LASER_INJECTION.pdf) is the full reference: deck elements
for both epoch2d and epoch3d, the binary file format and axis ordering,
the phase sign convention (including the LASY conversion), two-channel
polarisation, and current limitations. Read that before writing a profile
generator — the array-ordering convention in particular is easy to get
backwards silently.

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
convention conversion documented in `DOCUMENTATION.md`. This is not yet
packaged as a standalone, reusable wrapper — each test script currently
does its own LASY setup — but the conversion logic itself is stable and
validated (see the phase convention section of the documentation).

## Known issues

Tracked in this repo's [Issues](https://github.com/ZheyuanChen/AELP/issues),
not here.
