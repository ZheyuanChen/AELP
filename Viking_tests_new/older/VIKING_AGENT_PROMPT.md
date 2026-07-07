# Viking Agent Prompt — Laser-Injection Validation Tests

## Context

This directory (`lasers_in_epoch/Viking_tests/`) contains three high-resolution
laser-injection validation tests, packaged for head-less execution on Viking.
They were developed and tested locally in `Project_EPOCH`; you are now running
them on Viking HPC.

These tests validate the **custom numerical laser profile injector** that was
added to EPOCH — the ability to inject an arbitrary spatiotemporal amplitude
envelope (and optionally a separate phase profile) from a `.dat` file, rather
than relying on EPOCH's built-in analytical `gauss()` profile. The three tests
are milestones for my PhD progress report:

1. **Test 1 (`test1_gaussian_2d/`)** — 2D Gaussian beam, profile-only
   comparison. The analytical run uses EPOCH's native `gauss(y,0,w_bnd)` +
   `t_profile`; the numerical run injects the identical Gaussian envelope from
   `temporal_spatial_profile.dat`. Both share the same analytical paraxial phase
   formula, so any difference isolates the injection-path fidelity (file I/O +
   bilinear interpolation). Resolution: 50 cells/λ, nx=1750, ny=1250.

2. **Test 2 (`test2_lasy_2x2/`)** — 2D LASY 2×2 controlled experiment. Five
   runs sharing the same beam parameters, forming a 2×2 matrix (amplitude
   source × phase source) plus a sanity baseline:

   | Run subdir               | Amplitude source | Phase source   |
   |--------------------------|-----------------|----------------|
   | `analytical`             | deck `gauss()`  | deck formula   |
   | `numerical` (REFERENCE)  | numerical `.dat`| deck formula   |
   | `lasy_amp_deck_phase`    | LASY `.dat`     | deck formula   |
   | `numerical_amp_lasy_phase`| numerical `.dat`| LASY `.dat`    |
   | `lasy`                   | LASY `.dat`     | LASY `.dat`    |

   The `numerical` run is the reference. The 2×2 isolates whether the error
   between LASY and numerical comes from amplitude or phase. The key result
   (from the laptop runs): amplitude contributes ~94% of the variance; the
   residual is the physical gap between LASY's angular-spectrum propagation and
   the paraxial Gaussian.
   Resolution: 40 cells/λ, nx=2800, ny=2000.

3. **Test 3 (`test3_gaussian_3d/`)** — 3D Gaussian beam, profile-only. The
   3D counterpart of Test 1: native `gauss(y,0,w0)*gauss(z,0,w0)` vs a 2D
   (y,z) spatial profile from `spatial_profile.dat` (with `t_profile` handled
   analytically). Resolution: 30 cells/λ, nx=300, ny=nz=480.

## EPOCH binary requirements

The EPOCH binaries used **must** include the custom profile injection
modifications from `epoch_dev`:

- **epoch2d** (Tests 1 & 2): needs `use_custom_profile`, `use_spatiotemporal_profile`,
  `profile_data_file`, **and** `use_phase_from_file` + `phase_data_file` (the
  phase-from-file feature, needed by Test 2's `lasy` and `numerical_amp_lasy_phase`
  runs).
- **epoch3d** (Test 3): needs `use_custom_profile`, `use_spatiotemporal_profile=F`
  (spatial-only mode), and `profile_data_file`.

The unified job script `job_epoch` (repo root) selects the binary via
`--dim 2d|3d` and `--build dev|official`. These tests **must** use
`--build dev`, which points to:
- epoch2d: `/users/pnd531/epoch_dev/epoch2d` (GCC/14.3.0, OpenMPI/5.0.8)
- epoch3d: `/users/pnd531/epoch_dev/epoch3d` (GCC/14.3.0, OpenMPI/5.0.8)

**Before running**, confirm that these `epoch_dev` binaries have been compiled
with the custom injector modifications (including the phase-from-file
extension). If they haven't, the epoch_dev source must be rebuilt. Ask the user
rather than attempting to rebuild EPOCH yourself.

## Workflow

### Step 0 — Python environment

The analysis scripts need: `numpy`, `scipy`, `matplotlib`, `xarray`,
`sdf-xarray`. The LASY generator (`test2_lasy_2x2/generate_lasy_profiles.py`)
additionally needs `lasy >= 0.7.0`. Check these are available in the project
venv or install them.

### Step 1 — Generate `.dat` profile files

The `.dat` files are **not committed** (they're large, and regenerating them is
fast). Run the generators before EPOCH:

```bash
# Test 1
python lasers_in_epoch/Viking_tests/test1_gaussian_2d/numerical/generate_profile.py

# Test 2 — BOTH generators must run; they share a common (n_y=2000, n_t=1000) grid
python lasers_in_epoch/Viking_tests/test2_lasy_2x2/generate_numerical_profile.py
python lasers_in_epoch/Viking_tests/test2_lasy_2x2/generate_lasy_profiles.py

# Test 3
python lasers_in_epoch/Viking_tests/test3_gaussian_3d/numerical/generate_profile.py
```

These write `.dat` files directly into the run subdirectories. Verify the files
appeared (e.g. `temporal_spatial_profile.dat` in `test1.../numerical/`,
`amplitude_lasy.dat` and `phase_lasy.dat` in `test2.../lasy/`, etc.).

### Step 2 — Run EPOCH via Slurm

Submit each run directory as a separate Slurm job using the unified
`./job_epoch` wrapper (run it directly on the login node — it calls `sbatch`
internally). All tests use `--build dev`:

```bash
# Test 1 (2 runs, epoch2d)
./job_epoch --dir lasers_in_epoch/Viking_tests/test1_gaussian_2d/analytical --build dev
./job_epoch --dir lasers_in_epoch/Viking_tests/test1_gaussian_2d/numerical  --build dev

# Test 2 (5 runs, epoch2d)
for d in analytical numerical lasy lasy_amp_deck_phase numerical_amp_lasy_phase; do
  ./job_epoch --dir lasers_in_epoch/Viking_tests/test2_lasy_2x2/$d --build dev
done

# Test 3 (2 runs, epoch3d)
./job_epoch --dir lasers_in_epoch/Viking_tests/test3_gaussian_3d/analytical --dim 3d --build dev
./job_epoch --dir lasers_in_epoch/Viking_tests/test3_gaussian_3d/numerical  --dim 3d --build dev
```

EPOCH writes `.sdf` output files into the target directory (alongside the deck).
Wait for all jobs to complete before running analysis.

**Important notes:**
- Test 2's five runs are independent — submit them all in parallel.
- Test 3 is 3D and much more memory-intensive. If the 30 cells/λ resolution is
  too light, the decks have a comment noting `nx=400, ny=nz=640` for 40 cells/λ
  — but this significantly increases compute time. You can pass `--mem 256G` to
  cap memory if needed, but the default (full-node ~500 G) safely covers it.
- Do **NOT** run EPOCH on the login node — `./job_epoch` handles `sbatch`
  submission for you.

### Step 3 — Run analysis

Once all `.sdf` files are in place:

```bash
python lasers_in_epoch/Viking_tests/test1_gaussian_2d/analyse.py
python lasers_in_epoch/Viking_tests/test2_lasy_2x2/analyse.py
python lasers_in_epoch/Viking_tests/test3_gaussian_3d/analyse.py
```

Each script writes its output to a `results/` subdirectory under its test:

```
results/
  figures/     ← PNG field maps, difference panels, line-outs, bar charts
  summary.txt  ← narrative + headline metrics (replaces the notebook prose)
  metrics*.csv ← per-snapshot or per-case numerical data
```

The analysis scripts accept an optional `base_dir` argument (defaults to the
script's own directory). They import shared helpers from
`common/viking_analysis_lib.py` via a relative `../common` path, so **keep the
directory structure together**.

Analysis is lightweight (reads SDF, computes differences, writes PNGs) — it can
run on a login node or as a short Slurm job.

### Step 4 — Review results

Read `summary.txt` in each test's `results/` folder for the headline numbers.
Expected from the laptop runs:
- **Test 1**: max relative difference < 0.1% (confirms the file-injection path
  is faithful).
- **Test 2**: sanity ~0.02%, isolate-amplitude ~1.9%, isolate-phase ~0.5%,
  full-LASY ~2.0%, quadrature matches full-LASY; amplitude share of variance
  ~94%.
- **Test 3**: similar to Test 1 but in 3D — max relative difference should be
  small (< 1%).

The higher resolution on Viking may shift these slightly, but the qualitative
picture should be the same.

## `.dat` file format reference

- **2D spatiotemporal** (Tests 1 & 2 amplitude): line 1 = `n_t n_y`, line 2 =
  y-coordinates, line 3 = t-coordinates, then `n_t` rows of `n_y` values.
  Values are the normalised amplitude envelope (0 to 1).
- **2D phase** (Test 2 phase): same format, values in radians.
- **2D spatial-only** (Test 3): line 1 = `n1 n2` (n_y, n_z), line 2 =
  y-coordinates, line 3 = z-coordinates, then `n2` rows of `n1` values.

## Known gotchas

1. **`profile = custom` bug**: do NOT add `profile = custom` to any deck that
   uses `profile_data_file`. It triggers a premature load of the default
   filename during deck parsing, before `profile_data_file` is read. The decks
   here already avoid this.

2. **LASY carrier-offset (piston phase)**: the LASY generator applies a
   carrier-offset fix: `phase_epoch = -(phase_lasy - phi_ref) + (-gouy)`, where
   `phi_ref` is the LASY phase at (r=0, t_centre). Without this, a constant
   offset of ~π from LASY's propagator causes a ~137% error. The generator
   prints `phi_ref` — if it's near ±π, the fix is working as expected.

3. **Grid matching**: Test 2's numerical and LASY generators **must** share the
   same (n_y, n_t) grid (currently n_y=2000, n_t=1000, y ∈ [-25, 25] µm,
   t ∈ [0, 90 fs]), because the isolate-phase cell mixes the numerical
   amplitude with the LASY phase. If you change the grid in one, change it in
   both.

4. **`separate_times=True`**: the analysis library does NOT use this
   `sdf_xarray` option (it causes issues with coordinate rescaling); it loads
   files with plain `sdfxr.open_mfdataset(files, data_vars=[...])`. Do not add
   `separate_times=True` to the analysis scripts.
