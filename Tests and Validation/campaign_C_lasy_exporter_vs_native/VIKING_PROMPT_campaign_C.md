# Viking prompt — Campaign C: lasy-fork exporter vs native injector

Run the four EPOCH simulations in `Tests and Validation/
campaign_C_lasy_exporter_vs_native/` (this folder) and the analysis
scripts. Read the folder's `README.md` first for the design and the
calibrated pass/fail expectations (short version: expect a smooth ~1–2 %
residual between the cells — that is physics, not a bug; ~200 % or a
uniform large offset IS a bug, stop and report).

## Prerequisites

1. EPOCH builds: current `epoch_dev` epoch2d and epoch3d (the
   spatiotemporal custom-laser paths must be present).
2. The 2D binary pair (`2d/lasy_file/laser_{amplitude,phase}.dat`) is in
   the repo. The 3D pair (~91 MB each) is not: either rsync it from the
   local machine, or generate on Viking —
   `pip install --user git+https://github.com/ZheyuanChen/lasy.git`
   then run `3d/generate_lasy_exporter.py` **via sbatch** (numpy FFTs on
   a ~180 MB envelope; not a login-node job). Verify against
   `lasy_file/laser_metadata.txt` and the deck's declared shape — EPOCH's
   only load-time check is file size = n_tr1 × n_tr2 × n_t × 8 bytes.

## Order of operations

1. **2D pair first** (cheap): `2d/native/`, `2d/lasy_file/`.
2. **3D MPI smoke test before the full 3D pair**: the epoch3d
   spatiotemporal path previously aborted under multi-rank domain
   decomposition (pre-slab-loading builds, commit 9d50a6a6 era). Run a
   small spatiotemporal deck (e.g. the 3D lasy_file deck cut to
   nx=64, ny=nz=64, t_end = 20 femto, with a matching small file) at the
   SAME rank count as the production run, and confirm it initialises and
   steps. Only then submit the full pair.
3. **3D pair**: `3d/native/`, `3d/lasy_file/`.
4. Analysis (all via sbatch — SDF reads are login-node-prohibited):
   `2d/scripts/analyse.py`, `2d/scripts/field_comparison.py`,
   `3d/scripts/analyse.py`. For a 3D pointwise field comparison, adapt
   campaign A's `3d/scripts/field_comparison.py` (REF/VAR → native /
   lasy_file).

## Resource estimates

Per CLAUDE.md sizing (verify against the actual decks before submitting):

- `2d/native`: 2D vacuum, 434×451 = 0.196M cells → `--ntasks 8, --time 2:00:00`
- `2d/lasy_file`: 2D vacuum, 434×451 = 0.196M cells → `--ntasks 8, --time 2:00:00`
- `3d/native`: 3D vacuum, 434×451×451 = 88.3M cells → `--ntasks 64, --time 8:00:00`
- `3d/lasy_file`: 3D vacuum, 88.3M cells → `--ntasks 64, --time 8:00:00`
- 3D smoke test: 64³ spatiotemporal → `--ntasks 64, --time 0:30:00`

Caveats (from CLAUDE.md — read them, they are binding):

- **Partition granularity is unresolved** (two Viking reports disagreed on
  whether sub-node requests share nodes; check
  `scontrol show partition <name>` for BOTH `OverSubscribe` and
  `Exclusive` before trusting `--ntasks 64`). Either way: **calibrate,
  don't guess** — submit ONE 3D run first, measure real cost with
  `sacct`, then size the second off that data point.
- **Do not truncate walltime below the deck's natural `t_end`** — a
  truncated final snapshot silently corrupts "final state" analyses. The
  estimates above are already padded; pad more if the smoke test suggests
  it.

## Reporting

Collect `results/` (metrics.csv, summary.txt, figures) per dimensionality
back into this folder, following the campaign A layout. Leave a
`Daily-log relay` in `HANDOFF.md` for the next local session (this repo
has no daily_log tree — do not create one).
