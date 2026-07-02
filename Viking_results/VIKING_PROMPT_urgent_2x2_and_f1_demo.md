# Viking agent prompt — urgent progress-report runs (injector 2x2 + f/1 LASY-vs-paraxial demo)

**Priority: run this before anything else.** Two test campaigns, 10 EPOCH
runs total, needed for a progress report under deadline. A third, larger
campaign (moderate-NA 3D field benchmark + 3D two-channel polarisation
test) was designed the same session but is **postponed** — see
`3D/HANDOVER.md` for that one; do not start it unless explicitly asked.

Both campaigns here use `--build dev` (the `epoch_dev` fork with the
custom file-injector mods — `use_custom_profile`, `use_spatiotemporal_
profile`, `use_phase_from_file`). **Before submitting anything: confirm
the Viking-side `epoch_dev` checkout is at or past commit `aa605cca`**
("Fix per-rank slab window vs default startup load balancing") — this
fixed a real MPI abort in the spatiotemporal file-injection path that
Test B below depends on. If the checkout predates it, `git pull` and
rebuild both `epoch_dev/epoch2d` and `epoch_dev/epoch3d` first.

## Campaign A — `injector_2x2_validation/` (amplitude/phase injector isolation)

**What it tests**: does the file-injection pipeline (binary format,
interpolation) reproduce EPOCH's own native deck-expression evaluation,
for the SAME closed-form Gaussian beam formula delivered two different
ways? Moderate NA (~0.265) deliberately chosen so native-deck output is
close to ground truth — this isolates injector-pipeline error, not
paraxial-approximation error (that's Campaign B's job). "Cast together"
across 2D and 3D (same beam, same box in both, side by side).

Originally specified as a full 2x2 (amplitude source x phase source, each
{file, deck}) but **one cell is architecturally unreachable** with the
current code — `use_phase_from_file` is only consulted inside the static-
path loader that itself unconditionally loads amplitude from file too, so
"native amplitude + file phase" can't be expressed. Reduced to a 3-cell
chain that gives equivalent diagnostic power (see `injector_2x2_
validation/3d/amp_deck_phase_file/README.md` for the full reasoning):

| Cell | Amplitude | Phase | Purpose |
|---|---|---|---|
| `amp_deck_phase_deck` | native | native | baseline / ground truth |
| `amp_file_phase_deck` | **file** | native | isolates amplitude-injector error |
| `amp_file_phase_file` | **file** | **file** | isolates phase-injector error (vs the row above) |

Both `2d/` and `3d/` have this same 3-cell layout.

### Steps
```bash
cd injector_2x2_validation/3d && python generate_profile.py   # writes .dat into amp_file_* dirs
cd ../2d          && python generate_profile.py
```
(Files are NOT committed — regenerate on Viking. `generate_profile.py`
imports `../physics_params.py`; run it from inside `3d/` or `2d/`.)

Then run each of the 6 decks (`echo . | mpirun -n N epoch{2,3}d` inside
each run directory, or via `job_epoch --dir <dir> --dim {2d,3d} --build dev`).

Analysis (after all 3 cells in a dimensionality have run):
```bash
python injector_2x2_validation/3d/analyse.py
python injector_2x2_validation/2d/analyse.py
```
Writes `results/{summary.txt,metrics.csv,figures/waist_vs_x.png}` in each.

### Resource estimates (Campaign A)

All 6 are vacuum (field-only), short (`t_end` 350 fs), small-to-moderate
grids.

| Run dir | Dim | Cells | `--ntasks` | `--time` |
|---|---|---|---|---|
| `injector_2x2_validation/3d/amp_deck_phase_deck` | 3D | 434x451x451 = 88.3M | 48 | 4:00:00 |
| `injector_2x2_validation/3d/amp_file_phase_deck` | 3D | 88.3M | 48 | 4:00:00 |
| `injector_2x2_validation/3d/amp_file_phase_file` | 3D | 88.3M | 48 | 4:00:00 |
| `injector_2x2_validation/2d/amp_deck_phase_deck` | 2D | 434x451 = 0.196M | 8 | 2:00:00 |
| `injector_2x2_validation/2d/amp_file_phase_deck` | 2D | 0.196M | 8 | 2:00:00 |
| `injector_2x2_validation/2d/amp_file_phase_file` | 2D | 0.196M | 8 | 2:00:00 |

(88.3M cells sits in the "3D vacuum 50-200M" band per the sizing table ->
32-64 tasks; picked 48 as the middle. The 2D runs are trivial in absolute
terms — 8 tasks and the walltime is generous padding, not a real estimate
of need.)

## Campaign B — `tight_focus_f1/` (LASY vs paraxial-formula, f/1 demo)

**What it tests**: at NA=0.5 (f/1), does the paraxial closed-form formula
still describe the beam, or does genuine non-paraxial (angular-spectrum)
propagation via LASY differ measurably? **Injector-only** — both `lasy/`
and `paraxial/` inject amplitude+phase from binary files (spatiotemporal
path); there is no native-deck comparison cell here (unlike Campaign A,
only one physically-meaningful formula exists per source: LASY's actual
propagated field, or the closed-form paraxial approximation — nothing to
isolate against a third "native" baseline).

Builds on the 30 June session's `dev_test/laser_profile_injection/
tight_focusing_f1/` (same LASY pipeline, same base parameters, same
skimage 2D-phase-unwrap fix) but: (1) new raw-binary file format (old text
format no longer read), (2) actually injects into EPOCH and runs (the
30 June work only ever compared raw numpy arrays, never ran EPOCH), (3)
adds a genuine 3D (y,z) version (30 June was 2D/y-only), (4) uses the
dimension-correct Gouy phase (HALF for 2D-slab, FULL for 3D
circularly-symmetric — a deliberate correction vs the 30 June 2D script;
Gouy phase is a spatially-uniform constant so this doesn't change that
session's amplitude/wing-structure conclusions, only the phase reference).

**A real sign bug was caught and fixed while building this** (see
`tight_focus_f1/2d/generate_paraxial.py`'s inline comment): a naive
transcription of the 30 June script's `phase = k0*y^2/(2*RC) - gouy`
formula into this session's `R_BND`-based convention (where `R_BND =
-RC`) needs the sign flipped (`-k0*y^2/(2*R_BND)`), not copied verbatim.
Missing the flip made the beam DIVERGE instead of converge — caught by a
local smoke test (`w(x)` monotonically increasing across the box instead
of dipping to a minimum near the focus), fixed, and re-verified
(`x_focus` and `w0` now land within ~10% of the target f/1 focus at
coarse smoke-test resolution). **The 3D paraxial generator was written
with the corrected sign from the start** (verified by inspection against
the fixed 2D formula) but has NOT been smoke-tested end-to-end in EPOCH
(only the 2D cells were locally run-checked this session, given time
constraints) — recommend a quick low-resolution local check of `3d/
paraxial/` before committing the full allocation, using the same waist-
scan method as the 2D check (Hilbert envelope along x, 2D Gaussian fit
per x-plane, minimum should land near `x_focus = 5.09 um`, `w0 = 0.637
um`). The 3D LASY generator (`3d/generate_lasy.py`, resampling LASY's
`rt`-mode output onto a Cartesian (y,z) grid via `rho=sqrt(y^2+z^2)`) is
new code, not previously used anywhere in this project — same
recommendation applies, with even less prior confidence since it's not
just a formula transcription.

### Steps
```bash
cd tight_focus_f1/2d && python generate_lasy.py && python generate_paraxial.py
cd ../3d              && python generate_lasy.py && python generate_paraxial.py
```
Requires `lasy>=0.7.0` and `scikit-image` (both already project
dependencies, see `pyproject.toml`). LASY propagation + phase unwrap takes
a little while (rt-mode, 800 radial points x 500 time points for 2D; 3D
resamples the same rt output onto a 400x400 Cartesian grid so the LASY
call itself is no more expensive, just the resampling loop).

Then run all 4 decks, and analyse (compare the two resulting field
evolutions — no analyse.py was written for Campaign B this session, given
time constraints; adapt `injector_2x2_validation/*/analyse.py`'s waist-
scan approach, or the existing `dev_test/laser_profile_injection/
tight_focusing_f1/compare_f1_vs_paraxial.py` for the metrics/figure style,
noting that script compares raw arrays, not EPOCH output, so it needs
adapting to read SDF fields instead).

### Resource estimates (Campaign B)

| Run dir | Dim | Cells | `--ntasks` | `--time` |
|---|---|---|---|---|
| `tight_focus_f1/3d/paraxial` | 3D | 326x672x672 = 147.2M | 64 | 6:00:00 |
| `tight_focus_f1/3d/lasy` | 3D | 147.2M | 64 | 6:00:00 |
| `tight_focus_f1/2d/paraxial` | 2D | 326x672 = 0.219M | 8 | 2:00:00 |
| `tight_focus_f1/2d/lasy` | 2D | 0.219M | 8 | 2:00:00 |

(147.2M cells is in the "3D vacuum 50-200M" band -> 32-64 tasks; picked
64 given it's toward the top of that band and this is the priority
result. Walltime padded per the project's own "don't truncate below
t_end" rule — `t_end=150 fs` at this resolution should be much faster
than 6h, but pad generously rather than trim, per CLAUDE.md.)

## After both campaigns: calibrate, don't just trust the table

Per this project's standing Viking guidance: launch one run at the
planned core count, check real cost with `sacct`, and size the rest of
the batch off that measured data point rather than the generic table
alone.

## Known caveats to carry into the progress report

1. Campaign A is a 3-cell chain, not a true 2x2 — see the architectural
   reason above. Worth mentioning explicitly rather than presenting it as
   a full 2x2 that happens to have one empty cell.
2. Campaign B's 3D cells are less validated than the 2D cells (see
   above) — flag results as preliminary until a smoke test confirms the
   3D generator, or note this explicitly if reporting 3D numbers.
3. The per-rank-slab MPI bug that Campaign B's spatiotemporal path
   depends on being fixed is now resolved (`aa605cca`) — but this was
   discovered and fixed mid-session; if Viking's checkout is stale,
   Campaign B's 3D runs specifically are at risk of the original abort.
