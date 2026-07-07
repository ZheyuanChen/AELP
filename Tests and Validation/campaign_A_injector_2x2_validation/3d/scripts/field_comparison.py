"""
Test A (3D) field comparison -- analytical vs numerical E_y, animated +
screenshots, normalised to a0. Same design as the 2D version
(injector_2x2_validation/2d/field_comparison.py) -- reference =
amp_deck_phase_deck (native deck), variant = amp_file_phase_file (full
file-injection pipeline), both solved by the identical EPOCH Maxwell
solver so the difference isolates injector-pipeline fidelity.

3D-specific: sliced at z=0 (the beam's symmetry midplane) to get a 2D
(x,y) field comparable to the 2D triptych/animation machinery -- after
slicing, structurally identical to a native 2D run's data, so this
REUSES the 2D library's save_triptych/animate_triptych functions.

Memory note: SDF stores each field as one contiguous per-file block (see
the "chunks separate stored chunks" warning sdf_xarray prints), so even a
z=0 slice requires reading the FULL 3D block into memory first -- slicing
doesn't reduce the read. The waist-scan analysis (analyse.py) only ever
touched ONE snapshot per cell and still needed --mem=120G; this script
needs ALL 15 snapshots for TWO cells, so it processes one SDF FILE AT A
TIME (opening/closing individually, keeping only the small z=0 slice
across iterations) rather than bulk-loading the whole time series, to
avoid multiplying that per-file overhead across 15 snapshots at once.

Must be run via sbatch (SDF reads are login-node-prohibited on Viking).
Usage: python field_comparison.py [base_dir]
"""
import sys
import glob
from pathlib import Path

import numpy as np
import xarray as xr
import sdf_xarray as sdfxr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "2D" / "common"))
sys.path.insert(0, str(HERE.parent.parent / "3D" / "common"))
sys.path.insert(0, str(HERE.parent))
import viking_analysis_lib as val
import viking_analysis_lib_3d as lib3d
import physics_params as P

BASE = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
REF_NAME, VAR_NAME = "amp_deck_phase_deck", "amp_file_phase_file"
REF_LABEL, VAR_LABEL = "Analytical (native deck)", "Numerical (file-injected)"
A0 = 9.1093837015e-31 * 299792458.0 * (2 * np.pi * 299792458.0 / P.LAMBDA0) / 1.602176634e-19

results, figures = lib3d.make_results_dirs(BASE)
results, figures = Path(results), Path(figures)
animations = results / "animations"
animations.mkdir(parents=True, exist_ok=True)


def load_z0_timeseries(run_dir):
    """One SDF file at a time: extract the z=0 slice + its own timestamp,
    then close before moving to the next file, so peak memory is bounded
    by one file's block, not all files' blocks at once."""
    files = sorted(glob.glob(str(run_dir / "*.sdf")))
    slices, times = [], []
    x = y = None
    for f in files:
        ds = sdfxr.open_mfdataset([f], separate_times=True, data_vars=["Electric_Field_Ey"])
        ey = ds["Electric_Field_Ey"]
        time_dim = [d for d in ey.dims if d.startswith("time")][0]
        iz0 = int(np.argmin(np.abs(ey["Z_Grid_mid"].values)))
        slab = ey.isel(Z_Grid_mid=iz0, **{time_dim: 0})
        if x is None:
            x = slab["X_Grid_mid"].values
            y = slab["Y_Grid_mid"].values
        slices.append(np.asarray(slab.values) / A0)
        times.append(float(ey[time_dim].values[0]))
        ds.close()
    order = np.argsort(times)
    data = np.stack([slices[i] for i in order])
    times_arr = np.array(times)[order]
    return xr.DataArray(data, dims=("time", "X_Grid_mid", "Y_Grid_mid"),
                        coords={"time": times_arr, "X_Grid_mid": x, "Y_Grid_mid": y})


print(f"Loading {REF_NAME} (z=0 slice, one file at a time)...")
ref = load_z0_timeseries(BASE / REF_NAME)
print(f"Loading {VAR_NAME} (z=0 slice, one file at a time)...")
var = load_z0_timeseries(BASE / VAR_NAME)

n_t = min(ref.sizes["time"], var.sizes["time"])
ref, var = ref.isel(time=slice(0, n_t)), var.isel(time=slice(0, n_t))
diff = var - ref
peak = float(np.abs(ref).max())

env_vals = val.hilbert_envelope_along_x(ref)
times_fs = ref.coords["time"].values * 1e15

sel = sorted(set([0, n_t // 4, n_t // 2, 3 * n_t // 4, n_t - 1]))
for idx in sel:
    t = times_fs[idx]
    val.save_triptych(ref.isel(time=idx), var.isel(time=idx), diff.isel(time=idx),
                      t, figures / f"field_comparison_t{idx:02d}.png",
                      ref_label=REF_LABEL, var_label=VAR_LABEL)

anim_path = val.animate_triptych(ref, var, diff, times_fs,
                                 animations / "field_comparison.mp4",
                                 ref_label=REF_LABEL, var_label=VAR_LABEL)

rows = []
for idx in range(n_t):
    m = val.difference_metrics(diff.isel(time=idx).values, env_vals[idx])
    rows.append((f"{times_fs[idx]:.2f}", f"{m['max_abs']:.6e}",
                f"{m['max_abs']/peak*100:.5f}", f"{m['max_rel']*100:.5f}"))
val.write_metrics_csv(results / "field_comparison_metrics.csv", rows,
                      ("time_fs", "max_abs_a0", "max_abs_pct_peak", "max_rel_pct"))

overall = val.difference_metrics(diff.values, env_vals)
summary = f"""Test A (3D) -- field comparison: analytical (native deck) vs numerical (file-injected)
========================================================================================
Reference  : {REF_NAME} ({REF_LABEL}), z=0 midplane slice
Comparison : {VAR_NAME} ({VAR_LABEL}), z=0 midplane slice
Peak |Ey|  : {peak:.4f} a0
Snapshots  : {n_t}

Overall metrics (vs analytical reference):
  max |difference|         : {overall['max_abs']:.6e} a0
  max |difference| / peak  : {overall['max_abs']/peak*100:.6f} %
  max relative difference  : {overall['max_rel']*100:.6f} %
  mean relative difference : {overall['mean_rel']*100:.6f} %

Animation: {anim_path}
Screenshots: field_comparison_t*.png in results/figures/
"""
val.write_text(results / "field_comparison_summary.txt", summary)
print(summary)
