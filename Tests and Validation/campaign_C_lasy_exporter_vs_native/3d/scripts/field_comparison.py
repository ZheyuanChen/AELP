"""
Campaign C (3D) field comparison -- native-deck vs lasy-exporter E_y,
z=0 midplane slice, animated + screenshots, normalised to a0.

Adapted from campaign A's 3d/scripts/field_comparison.py with REF/VAR
set to native/lasy_file. Same design as campaign C's own 2D
field_comparison.py -- reference = native
(EPOCH's own gauss()/quadratic-phase deck expression), variant =
lasy_file (the lasy fork's write_to_file exporter, epoch3d
spatiotemporal path), both solved by the identical EPOCH Maxwell solver
so the difference isolates (exporter + lasy-vs-paraxial physics)
fidelity, not solver discretisation error.

INTERPRETATION (differs from campaign A): campaign A's file/deck cells
shared one closed-form source, so its pointwise residual (~1e-3 % of
peak) was the injector floor. Campaign C's sources genuinely differ
(angular-spectrum envelope vs paraxial formula), so expect a smooth
residual of order 1-2 % of peak, amplitude-dominated, consistent with
the 2D twin's measured ~2.5% w(x) RMS / ~3.6% peak-normalised field
residual. A residual near 200 % of peak means a sign-flipped carrier; a
spatially uniform few-10s-of-% residual means a CEP/carrier-pin mismatch
(check carrier_phase_ref = psi_bnd_3d in the generator against the
native deck's phase_const). See the 2D twin's module docstring for the
full failure-signature list.

Memory note (from campaign A's docstring, carried over unchanged): SDF
stores each field as one contiguous per-file block, so even a z=0 slice
requires reading the FULL 3D block into memory first -- slicing doesn't
reduce the read. Campaign A's waist-scan analysis touched only ONE
snapshot per cell and still needed --mem=120G; this script processes one
SDF FILE AT A TIME (opening/closing individually, keeping only the small
z=0 slice across iterations) to avoid multiplying that overhead across
all 15 snapshots x 2 cells at once, but still needs a generously-sized
sbatch job (see run_analysis.sbatch: --mem=150G).

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
sys.path.insert(0, str(HERE.parent.parent.parent / "shared_libraries"))
sys.path.insert(0, str(HERE.parent.parent))
import viking_analysis_lib as val
import physics_params as P

BASE = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parent
REF_NAME, VAR_NAME = "native", "lasy_file"
REF_LABEL = "Native deck (paraxial expression)"
VAR_LABEL = "lasy exporter (file-injected)"

results, figures = val.make_results_dirs(BASE)
results, figures = Path(results), Path(figures)
animations = Path(BASE) / "results" / "animations"
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
        slices.append(np.asarray(slab.values) / val.a0_norm(P.LAMBDA0))
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
val.write_metrics_csv(Path(results) / "field_comparison_metrics.csv", rows,
                      ("time_fs", "max_abs_a0", "max_abs_pct_peak", "max_rel_pct"))

overall = val.difference_metrics(diff.values, env_vals)
summary = f"""Campaign C (3D) -- field comparison: native deck vs lasy-exporter file injection
================================================================================
Reference  : {REF_NAME} ({REF_LABEL}), z=0 midplane slice
Comparison : {VAR_NAME} ({VAR_LABEL}), z=0 midplane slice
Peak |Ey|  : {peak:.4f} a0
Snapshots  : {n_t}

Overall metrics (vs native reference):
  max |difference|         : {overall['max_abs']:.6e} a0
  max |difference| / peak  : {overall['max_abs']/peak*100:.6f} %
  max relative difference  : {overall['max_rel']*100:.6f} %
  mean relative difference : {overall['mean_rel']*100:.6f} %

Expected: ~1-2 % of peak, smooth and amplitude-dominated (the known
lasy-vs-paraxial physics gap at NA~0.265, consistent with the 2D twin's
measured residual) -- NOT ~200 % (sign-flipped carrier) or a uniform
few-10s-of-% offset (carrier-pin/CEP mismatch); see module docstring.

Animation: {anim_path}
Screenshots: field_comparison_t*.png in results/figures/
"""
val.write_text(Path(results) / "field_comparison_summary.txt", summary)
print(summary)
