# Viking laser-injection validation tests

High-resolution versions of the three laser-injection comparison milestones,
packaged for running **head-less on Viking** (or any HPC node). The laptop
versions used Jupyter notebooks; here the analysis is done by plain Python
scripts that write PNG figures and plain-text/CSV results into a `results/`
folder per test — no display, no animations.

## The three tests

| Dir | What it shows |
|-----|---------------|
| `test1_gaussian_2d/` | **Milestone 1** — 2D Gaussian beam, profile-only. Custom `.dat` amplitude injection vs EPOCH's native analytical `gauss()` (same analytical phase in both). Validates the injector. |
| `test2_lasy_2x2/` | **Milestone 2** — 2D LASY 2×2 (amplitude source × phase source). Isolates whether the LASY-vs-numerical error comes from amplitude or phase; includes the phase-from-file path. |
| `test3_gaussian_3d/` | **Milestone 3** — 3D counterpart of Test 1. Native `gauss(y)·gauss(z)` vs a 2D (y,z) `.dat` spatial profile. |

`common/viking_analysis_lib.py` holds shared head-less analysis helpers used by
all three `analyse.py` scripts. **Keep the folder together** when moving it — the
analysis scripts import the shared lib via a relative path (`../common`).

## Requirements

* EPOCH built with the custom-injector modifications:
  * epoch2d for tests 1 and 2 (amplitude **and** phase from file — the latter
    needs the `use_phase_from_file` / `phase_data_file` deck elements)
  * epoch3d for test 3 (spatial-only custom profile)
* Python: `numpy scipy matplotlib xarray sdf-xarray` for analysis;
  additionally `lasy` (>= 0.7.0) for the Test 2 LASY generator.

## Workflow on Viking

For each run directory that has a generator, **generate the `.dat` files first**
(they are not committed — regenerate them on the target machine):

```bash
# Test 1
python test1_gaussian_2d/numerical/generate_profile.py

# Test 2 (writes into the relevant run subdirs; both share one (y,t) grid)
python test2_lasy_2x2/generate_numerical_profile.py
python test2_lasy_2x2/generate_lasy_profiles.py

# Test 3
python test3_gaussian_3d/numerical/generate_profile.py
```

Then run EPOCH in each run subdirectory (output `.sdf` files land alongside the
deck — the analysis also accepts a `sdf_files/` subdir). For example with a
Slurm job step:

```bash
cd test1_gaussian_2d/analytical && mpirun epoch2d <<< "."
cd ../numerical               && mpirun epoch2d <<< "."
```

Finally run the analysis (writes `results/` under each test):

```bash
python test1_gaussian_2d/analyse.py
python test2_lasy_2x2/analyse.py
python test3_gaussian_3d/analyse.py
```

Each `analyse.py` also accepts an explicit base directory as its first argument,
in case you run from elsewhere.

## Resolution knobs

These are deliberately higher-resolution than the laptop runs. To rescale,
edit `nx`/`ny`(/`nz`) in the decks **and** the matching `n_y`/`n_t` (2D) or
`n1`/`n2` (3D) in the generators so the `.dat` grid still resolves the boundary
at least as finely as the simulation grid.

* Test 1: 50 cells/λ → `nx=1750, ny=1250` on a 35 µm × 25 µm box.
* Test 2: 40 cells/λ → `nx=2800, ny=2000` on a 70 µm × 50 µm box; `.dat` grid
  `n_y=2000, n_t=1000` (the numerical and LASY generators **must** share this
  grid — the isolate-phase cell mixes the analytical amplitude with the LASY
  phase).
* Test 3: 30 cells/λ → `nx=300, ny=nz=480` on a 10 µm × 16 µm × 16 µm box
  (the expensive one; `nx=400, ny=nz=640` gives 40 cells/λ).

## Outputs

Each test writes `results/`:
* `figures/*.png` — field maps (reference | variant | difference), relative-
  difference panels, on-axis line-outs; plus the 2×2 bar chart (Test 2) and a
  z-slice montage (Test 3).
* `summary.txt` — narrative + headline metrics (the notebook's prose, head-less).
* `metrics*.csv` — per-snapshot (or per-case) max absolute / relative difference.
