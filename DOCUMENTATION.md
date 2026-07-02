# Custom Laser Profile Block

This document describes the modifications made to EPOCH's `laser` block to
support loading arbitrary laser amplitude and phase profiles from [raw binary
data files](https://epochpic.github.io/documentation/input_deck/binary_files.html). 
See the [EPOCH documentation](https://epochpic.github.io/documentation/input_deck/input_deck_laser.html) 
for the standard `laser` block parameters that are not covered here.

The custom profile injection feature allows a laser field of the form
`E = A · sin(ωt + φ)` — with `A` and `φ` arbitrary functions of the
transverse boundary coordinates and (optionally) time — to be driven
entirely from externally generated data files, bypassing EPOCH's built-in
analytic profile functions. This is particularly useful when interfacing
with external beam propagation codes such as
[LASY](https://github.com/LASY-org/lasy), or when the target profile cannot
be expressed as one of EPOCH's analytic primitives. A Python wrapper for LASY 
written by the same author can be found [here](https://github.com/ZheyuanChen/AELP/tree/main)

The feature is implemented in **epoch2d** (released as epoch_dev `v2.0.0`)
and **epoch3d** (released as epoch_dev `v2.1.0`). epoch1d is unmodified.

Two modes exist, selected per laser block:

- **Spatiotemporal** (`use_spatiotemporal_profile = T`): the full
  amplitude `A(transverse, t)` is interpolated from file at every time
  step. In 2D the data is `A(y, t)`; in 3D it is `A(tr1, tr2, t)`.
- **Static spatial** (`use_spatiotemporal_profile = F`): a
  time-independent spatial profile is interpolated onto the boundary once
  at start-up; the temporal envelope is set analytically in the deck via
  `t_profile`.

In both modes the phase can also be read from file with
`use_phase_from_file = T`, on the same grid as the amplitude.

**Note on defaults:** when `use_custom_profile = T` and
`use_spatiotemporal_profile` is not given, epoch2d defaults to the
spatiotemporal mode (`T`) while epoch3d defaults to the static spatial mode
(`F`). Setting the flag explicitly is recommended.

## EPOCH2D

An example 2D `laser` block using all custom profile features:

```
begin:laser
  boundary = x_min
  intensity_w_cm2 = 1e22
  lambda = 1.0 * micron
  pol = 0

  use_custom_profile = T
  use_spatiotemporal_profile = T
  profile_data_file = amplitude.dat

  n_t = 50
  n_y = 30
  y_min = -10.0 * micron
  y_max =  10.0 * micron
  t_start = 0.0
  t_end = 30e-15

  use_phase_from_file = T
  phase_data_file = phase.dat
end:laser
```

A 2D custom-profile laser block accepts the following parameters:

`use_custom_profile` — Enable custom profile injection. Must be set to `T`
to activate any of the parameters below. The default is `F`, which leaves
all standard EPOCH profile behaviour unchanged.

`use_spatiotemporal_profile` — When `T`, the full two-dimensional
`A(y, t)` amplitude is read from a raw binary file. When `F`, a
one-dimensional static spatial profile `A(y)` is used instead (in epoch2d
this legacy path still reads the old **text** format of coordinate–value
pairs, not binary). Defaults to `T` in epoch2d.

`profile_data_file` — Path to the amplitude file. Relative paths resolve
from the deck's data directory; absolute paths are used as-is. Defaults to
`temporal_spatial_profile.dat` (spatiotemporal) or `spatial_profile.dat`
(static) if omitted.

`n_t_points` / `n_t` — Number of time points in the profile array.
Mandatory when `use_spatiotemporal_profile = T`.

`n_transverse_points` / `n_y` — Number of transverse spatial points in the
profile array. Mandatory when `use_spatiotemporal_profile = T`.

`profile_transverse_min` / `y_min` — Minimum transverse coordinate of the
profile grid (metres). Mandatory when `use_spatiotemporal_profile = T`.
(For a laser on `y_min`/`y_max` the transverse axis is physically `x`; the
`y_min` alias name is kept for backward compatibility, but the
boundary-agnostic `profile_transverse_min` reads better there.)

`profile_transverse_max` / `y_max` — Maximum transverse coordinate of the
profile grid (metres). Mandatory when `use_spatiotemporal_profile = T`.

`t_start`, `t_end` — Start and end times of the profile grid (seconds).
These reuse EPOCH's standard laser elements; the profile time axis runs
from `t_start` to `t_end` with `n_t` uniformly spaced points. Outside this
window the laser is inactive.

`use_phase_from_file` — When `T`, the phase `φ(y, t)` is also read from
file each time step. Any `phase = ...` expression in the deck is ignored.
Defaults to `F`, in which case the phase is set as normal via the `phase`
deck element.

`phase_data_file` — Path to the raw binary phase file (radians). Follows
the same path resolution as `profile_data_file`. The phase file must lie on
the same grid as the amplitude file (same `n_t`, `n_y`, bounds, `t_start`,
`t_end`). Defaults to `phase_profile.dat` if omitted.

## EPOCH3D

In 3D the boundary plane has two transverse axes, so the grid declaration
uses two transverse point counts and two pairs of bounds. The
boundary-agnostic names are `tr1` and `tr2`, whose physical meaning depends
on the boundary the laser is attached to:

| Boundary        | tr1 axis | tr2 axis |
|-----------------|----------|----------|
| `x_min`/`x_max` | y        | z        |
| `y_min`/`y_max` | x        | z        |
| `z_min`/`z_max` | x        | y        |

An example 3D `laser` block (spatiotemporal, amplitude + phase, on
`x_min`, written with the axis-named aliases):

```
begin:laser
  boundary = x_min
  intensity_w_cm2 = 1e22
  lambda = 1.0 * micron
  pol = 0

  use_custom_profile = T
  use_spatiotemporal_profile = T
  profile_data_file = amplitude.dat

  n_y_points = 48
  n_z_points = 40
  n_t_points = 16
  y_min = -8.0 * micron
  y_max =  8.0 * micron
  z_min = -8.0 * micron
  z_max =  8.0 * micron
  t_start = 0.0
  t_end = 60e-15

  use_phase_from_file = T
  phase_data_file = phase.dat
end:laser
```

The same block written with the boundary-agnostic names (useful for decks
generated by boundary-independent tooling):

```
  n_tr1 = 48
  n_tr2 = 40
  n_t = 16
  tr1_min = -8.0 * micron
  tr1_max =  8.0 * micron
  tr2_min = -8.0 * micron
  tr2_max =  8.0 * micron
```

A 3D custom-profile laser block accepts the following parameters (in
addition to `use_custom_profile`, `use_spatiotemporal_profile`,
`profile_data_file`, `use_phase_from_file`, `phase_data_file`, `t_start`
and `t_end`, which behave exactly as in 2D):

`n_transverse1_points` / `n_tr1` — Number of points along the first
transverse axis (see the table above). Mandatory (>= 2) whenever
`use_custom_profile = T`.

`n_transverse2_points` / `n_tr2` — Number of points along the second
transverse axis. Mandatory (>= 2) whenever `use_custom_profile = T`.

`n_t_points` / `n_t` — Number of time points. Mandatory (>= 2) when
`use_spatiotemporal_profile = T`; ignored on the static spatial path.

`profile_transverse1_min` / `tr1_min`, `profile_transverse1_max` /
`tr1_max` — Bounds of the profile grid along the first transverse axis
(metres). Mandatory whenever `use_custom_profile = T`.

`profile_transverse2_min` / `tr2_min`, `profile_transverse2_max` /
`tr2_max` — Bounds along the second transverse axis (metres). Mandatory
whenever `use_custom_profile = T`.

`n_{x,y,z}_points` / `n_{x,y,z}`, `{x,y,z}_min`, `{x,y,z}_max` —
Axis-named aliases for the point counts and bounds above. These are
resolved to tr1 or tr2 according to the boundary table, so a deck for an
`x_min` laser can say `n_y_points`/`n_z_points` and `y_min`...`z_max`
naturally. Naming the propagation axis (e.g. `n_x_points` on an `x_min`
laser) is an error — EPOCH reports the offending element and line and
suggests the transverse or boundary-agnostic names.

Unlike epoch2d, the epoch3d **static spatial path also uses the raw binary
format** (the pre-`v2.1.0` text format with embedded coordinates is no
longer read), and supports `use_phase_from_file`: the phase plane is
interpolated onto the boundary once at start-up alongside the amplitude.

## File Format

All data files are raw binary with no embedded header. Shape and grid
bounds are declared entirely in the deck (see above). This matches EPOCH's
documented binary-file convention used elsewhere in the code (e.g.
`particles_from_file`). EPOCH checks the file size against the declared
shape at load time and aborts with the expected and actual byte counts if
they disagree. Note that the size check cannot detect a transposed array of
the same total size — the axis ordering below must be respected.

Each file contains an array of 64-bit IEEE 754 floating point values (the
same precision as EPOCH's `REAL(num)` type), written in Fortran-compatible
column-major order:

- **2D spatiotemporal** `A(y, t)`: transverse `y` fastest-varying, time
  slowest. `n_y * n_t` values.
- **3D spatiotemporal** `A(tr1, tr2, t)`: tr1 fastest-varying, then tr2,
  time slowest. `n_tr1 * n_tr2 * n_t` values.
- **3D static spatial** `A(tr1, tr2)`: tr1 fastest-varying. `n_tr1 * n_tr2`
  values.

All grids must be uniform; EPOCH reconstructs the axes from the declared
bounds and point counts (`np.linspace` on the Python side). Amplitude
values should be normalised to the range `[0, 1]`; EPOCH multiplies by the
peak field derived from `intensity_w_cm2` (or `amp`). Phase values are in
radians and are used directly — EPOCH adds no constant offset.

Interpolation is bilinear (2D) or trilinear (3D spatiotemporal) between
grid points. On the spatiotemporal path, positions outside the declared
transverse bounds return **zero** amplitude and phase; on the 3D static
spatial path, positions outside the plane are **clamped** to the edge
values. In either case the profile grid should comfortably cover the
boundary region the laser is meant to illuminate.

### Python generators

The Fortran column-major layouts above correspond to NumPy arrays in
default C (row-major) order with the axes reversed, written with
`.tofile()` — **no transpose is needed** in either dimension:

| File | NumPy shape (C order) |
|------|----------------------|
| 2D spatiotemporal | `(n_t, n_y)` |
| 3D spatiotemporal | `(n_t, n_tr2, n_tr1)` |
| 3D static spatial | `(n_tr2, n_tr1)` |

2D example:

```python
import numpy as np

N_T   = 50
N_Y   = 30
Y_MIN, Y_MAX   = -10e-6, 10e-6
T_START, T_END = 0.0, 30e-15

y = np.linspace(Y_MIN, Y_MAX, N_Y)
t = np.linspace(T_START, T_END, N_T)

# meshgrid with indexing="ij" gives arr[t_index, y_index]
Tg, Yg = np.meshgrid(t, y, indexing="ij")

# --- fill amplitude (values in [0, 1]) ---
t0, tau = 15e-15, 5e-15
w0 = 3e-6
amplitude = np.exp(-((Tg - t0) / tau)**2) * np.exp(-(Yg / w0)**2)

assert amplitude.shape == (N_T, N_Y)
amplitude.astype(np.float64).tofile("amplitude.dat")

# --- fill phase (values in radians) ---
phase = np.zeros((N_T, N_Y))  # flat phase for a plane wave
phase.astype(np.float64).tofile("phase.dat")
```

3D example (for an `x_min` laser, so tr1 = y and tr2 = z):

```python
import numpy as np

N_Y, N_Z, N_T  = 48, 40, 16          # n_tr1, n_tr2, n_t in the deck
Y_MIN, Y_MAX   = -8e-6, 8e-6
Z_MIN, Z_MAX   = -8e-6, 8e-6
T_START, T_END = 0.0, 60e-15

y = np.linspace(Y_MIN, Y_MAX, N_Y)
z = np.linspace(Z_MIN, Z_MAX, N_Z)
t = np.linspace(T_START, T_END, N_T)

# indexing="ij" gives arr[t_index, z_index, y_index] = (n_t, n_tr2, n_tr1)
Tg, Zg, Yg = np.meshgrid(t, z, y, indexing="ij")

t0, tau = 30e-15, 10e-15
w0 = 3e-6
amplitude = (np.exp(-((Tg - t0) / tau)**2)
             * np.exp(-(Yg**2 + Zg**2) / w0**2))

assert amplitude.shape == (N_T, N_Z, N_Y)
amplitude.astype(np.float64).tofile("amplitude.dat")

phase = np.zeros((N_T, N_Z, N_Y))
phase.astype(np.float64).tofile("phase.dat")
```

## Phase Convention

EPOCH injects the field as:

```
E ∝ profile · sin( ω·t + phase )
```

where `ω·t` (the carrier phase) is accumulated internally by EPOCH. The
phase file must therefore contain only the slowly varying envelope phase,
not the carrier. This applies identically in 2D and 3D.

When using a complex envelope from LASY (where the physical field is
`E = Re[E · exp(−iω₀t)] = |E| · cos(ω₀t − φ)`), the required conversion is:

```python
phase_epoch = -(phase_lasy - phi_ref) - gouy_phase
```

where `phi_ref = angle(E[on axis, t == t_peak])` is the on-axis phase at
peak time (which removes the arbitrary global piston phase introduced by
the propagator), and `gouy_phase` pins the result to EPOCH's carrier
convention. A `π/2` offset converts `cos → sin` but as a constant it does
not affect field shape — the sign and reference subtraction are what matter
physically.

See `PHASE_INJECTION.md` (in the `epoch_dev` repository) for a detailed
derivation and validation, including a description of the 137% error that
results from omitting the reference subtraction.

## Two-Channel Polarisation

Two custom-file laser blocks on the same boundary with `pol = 0` and
`pol = 90` inject into the two transverse polarisation channels
independently (e.g. `Ey` and `Ez` for an `x_min` laser):

```
begin:laser
  boundary = x_min
  ...
  pol = 0
  profile_data_file = amplitude_H.dat   # drives Ey
  ...
end:laser

begin:laser
  boundary = x_min
  ...
  pol = 90
  profile_data_file = amplitude_V.dat   # drives Ez
  ...
end:laser
```

Each block declares its own grid and carries independent amplitude and
phase data — the file data is stored per laser block, so there is no
aliasing between the two channels even if the grids differ. The
`cos`/`sin(pol_angle)` channel split is identical in every boundary routine
of epoch2d and epoch3d, so the same pattern works in 3D. (The two-channel
configuration has been validated end-to-end in 2D; the 3D per-laser storage
is the same design, built in from the start.)

## Current Limitations

- **epoch1d** has no custom-profile support.
- The **epoch2d static spatial path** (`use_spatiotemporal_profile = F`)
  still reads the legacy 1D text format, not binary. In epoch3d both paths
  are binary.
- All file grids must be **uniform**; non-uniform coordinate axes are not
  supported.
- Memory: in **3D**, each MPI rank stores only the slab of the
  spatiotemporal array covering its own patch of the laser boundary
  (plus a one-cell margin); ranks that do not own the laser's boundary
  face store nothing, and loading streams the file one time-slice at a
  time, so no rank ever holds the full 3D array. Per-rank cost is
  roughly `(local patch fraction) * n_tr1 * n_tr2 * n_t * 8` bytes. Two
  fallbacks store the full plane on boundary ranks: dynamic load
  balancing, and a moving window when the laser is on a `y`/`z` boundary
  (both make the local patch time-dependent). In **2D** the (small) full
  array is still broadcast to every rank.
- The 3D implementation (epoch_dev `v2.1.0`) has been verified
  functionally (loading, axis-ordering, error paths); a quantitative
  field-level benchmark in 3D (matching the 2D LASY validation) has not
  yet been run.
