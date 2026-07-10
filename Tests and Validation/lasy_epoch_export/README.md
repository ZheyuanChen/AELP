# lasy → EPOCH raw-binary export tests

Validation suite for the `file_format="epoch"` exporter added to the
[lasy fork](https://github.com/ZheyuanChen/lasy) (`lasy/utils/epoch_helper.py`
plus the `Laser.write_to_file` dispatch). The exporter converts a lasy
complex envelope into the headerless amplitude/phase `.dat` pair consumed
by the `use_spatiotemporal_profile` injection path of the
[epoch_dev](https://github.com/ZheyuanChen/epoch_dev) fork.

## What is verified

| Test | Convention checked |
|---|---|
| `test_file_sizes_match_headerless_convention` | file size = n_x·n_y·n_t × 8 bytes (EPOCH's only load-time check) |
| `test_amplitude_normalised_to_unit_peak` | profile ∈ [0, 1], peak exactly 1 (deck `amp` carries the physical scale) |
| `test_fortran_layout_x_fastest` | C-order (t, y, x) on disk = Fortran (x, y, t), matching epoch3d's `ALLOCATE(matrix(n_tr1, n_tr2, n_t))` |
| `test_field_reconstruction_matches_lasy` | `amp·sin(ω₀t + phase)` rebuilt from the files reproduces lasy's `Re[E_env e^{−iω₀t}]` pointwise (catches phase-sign and ±π/2 carrier-offset errors) |
| `test_phase_at_peak_is_carrier_offset` | stored phase at the amplitude peak is exactly +π/2 (cos → sin shift, φ_ref referencing) |
| `test_phase_is_unwrapped_where_it_matters` | no 2π wrap seams between neighbouring samples where amplitude > 10⁻³ of peak (EPOCH bilinearly interpolates the phase file) |
| `test_metadata_sidecar_reports_deck_parameters` | sidecar reports peak field (V/m), point counts, transverse extents, time window — everything the deck must declare |
| `test_rt_geometry_is_rejected` | cylindrical (`rt`) grids refuse loudly rather than writing a misshapen file |

## Conventions under test

```
lasy : E_phys = Re[E_env · exp(−iω₀t)] = |E_env| · cos(ω₀t − φ)
EPOCH: E_phys = amp · profile · sin(∫ω dt + phase)   (laser.f90)

⇒ profile = |E_env| / max|E_env|
  phase   = −(φ − φ_ref) + π/2
```

φ_ref (envelope phase at the amplitude peak) removes the arbitrary global
phase that lasy's angular-spectrum propagator adds; the absolute
carrier-envelope phase is then set by the deck's own `phase` parameter.

## Running

Requires the lasy fork importable (editable install) plus `numpy`, `scipy`
and `pytest`. No EPOCH build or SDF files needed — the suite is pure
Python and runs in about a second:

```bash
pytest test_lasy_epoch_export.py -v
```

The suite was negative-control checked: flipping the carrier offset back
to −π/2 (the bug it was written to guard against — a globally sign-flipped
field, invisible in intensity but a ~200% pointwise field error) fails
`test_field_reconstruction_matches_lasy` and
`test_phase_at_peak_is_carrier_offset`.
