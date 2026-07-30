# LASY → EPOCH raw-binary export baseline tests

Validation suite for the `file_format="epoch"` (3D, epoch3d) and
`file_format="epoch2d"` (lasy-y slice, epoch2d) exporters added to the
[lasy fork](https://github.com/ZheyuanChen/lasy) (`lasy/utils/epoch_helper.py`
plus the `Laser.write_to_file` dispatch). The exporters convert a lasy
complex envelope into the headerless amplitude/phase `.dat` pair consumed
by the `use_spatiotemporal_profile` injection path of the
[epoch_dev](https://github.com/ZheyuanChen/epoch_dev) fork.

> **Version boundary:** this 16-test AELP suite records the initial Cartesian
> (`xyt`) exporter at LASY commit `83d2b3e`. Its
> `test_rt_geometry_is_rejected` is intentionally historical. Commit
> `aaacc68` subsequently added supported mode-0 `rt` export, with seven
> current `rt`/`xyt` tests in
> [`lasy/tests/test_epoch_export.py`](https://github.com/ZheyuanChen/lasy/blob/dev-zheyuan/tests/test_epoch_export.py).
> Run the suite in this folder against `83d2b3e`, or run the LASY fork's own
> tests against `aaacc68` and later. Campaign C provides the current
> simulation-level `rt` validation.

## What is verified

3D path (`file_format="epoch"`):

| Test | Convention checked |
|---|---|
| `test_file_sizes_match_headerless_convention` | file size = n_x·n_y·n_t × 8 bytes (EPOCH's only load-time check) |
| `test_amplitude_normalised_to_unit_peak` | profile ∈ [0, 1], peak exactly 1 (deck `amp` carries the physical scale) |
| `test_fortran_layout_x_fastest` | C-order (t, y, x) on disk = Fortran (x, y, t), matching epoch3d's `ALLOCATE(matrix(n_tr1, n_tr2, n_t))` |
| `test_field_reconstruction_matches_lasy` | `amp·sin(ω₀t + phase)` rebuilt from the files reproduces lasy's `Re[E_env e^{−iω₀t}]` pointwise (catches phase-sign and ±π/2 carrier-offset errors) |
| `test_phase_at_peak_is_carrier_offset` | stored phase at the amplitude peak is exactly +π/2 (cos → sin shift, φ_ref referencing) |
| `test_phase_is_unwrapped_where_it_matters` | no 2π wrap seams between neighbouring samples where amplitude > 10⁻³ of peak (EPOCH bilinearly interpolates the phase file) |
| `test_metadata_sidecar_reports_deck_parameters` | sidecar reports peak field (V/m), epoch3d deck element names (n_tr1/n_tr2/n_t_points, tr extents, time window) — everything the deck must declare |
| `test_rt_geometry_is_rejected` | cylindrical (`rt`) grids refuse loudly rather than writing a misshapen file (both formats) |

2D path (`file_format="epoch2d"`, for epoch2d's `(n_transverse, n_t)` reader —
the retained lasy x axis maps to EPOCH2d's transverse y):

| Test | Convention checked |
|---|---|
| `test_2d_file_sizes_match_headerless_convention` | file size = n_transverse·n_t × 8 bytes |
| `test_2d_layout_matches_y0_slice` | C-order (t, x) on disk = the y=0 plane of the envelope, Fortran (n_transverse, n_t), transverse fastest |
| `test_2d_field_reconstruction_matches_lasy` | pointwise physics check on the slice |
| `test_2d_phase_at_peak_is_carrier_offset` | +π/2 at the slice's peak |
| `test_2d_phase_is_unwrapped_where_it_matters` | seam-free 2D phase |
| `test_2d_consistent_with_3d_export` | 2D export ≡ y=0 plane of the 3D export, sample by sample (mod 2π in phase) |
| `test_2d_metadata_reports_slice_and_deck_parameters` | epoch2d deck element names (n_y, y_min/max, n_t_points, t window) + the actual slice coordinate used |
| `test_carrier_phase_ref_shifts_stored_phase` | `carrier_phase_ref` adds exactly that constant (the only way to set the CEP/Gouy pin, since EPOCH ignores the deck `phase` expression when the phase comes from file) |

## Conventions under test

```
lasy : E_phys = Re[E_env · exp(−iω₀t)] = |E_env| · cos(ω₀t − φ)
EPOCH: E_phys = amp · profile · sin(∫ω dt + phase)   (laser.f90)

⇒ profile = |E_env| / max|E_env|
  phase   = −(φ − φ_ref) + π/2 + carrier_phase_ref
```

φ_ref (envelope phase at the amplitude peak) removes the arbitrary global
phase that LASY's angular-spectrum propagator adds; the absolute
carrier-envelope phase is set by `carrier_phase_ref` in the file because
EPOCH ignores the deck `phase` expression when `use_phase_from_file = T`.

## Running

Requires the LASY fork at `83d2b3e` importable (editable install) plus
`numpy`, `scipy`, and `pytest`. No EPOCH build or SDF files are needed — the
suite is pure Python and runs in about a second:

```bash
pytest test_lasy_epoch_export.py -v
```

For the current `dev-zheyuan` branch (`aaacc68` or later), run the fork's
version-matched suite instead:

```bash
pytest tests/test_epoch_export.py -v
```

The suite was negative-control checked: flipping the carrier offset back
to −π/2 (the bug it was written to guard against — a globally sign-flipped
field, invisible in intensity but a ~200% pointwise field error) fails
`test_field_reconstruction_matches_lasy` and
`test_phase_at_peak_is_carrier_offset`.
