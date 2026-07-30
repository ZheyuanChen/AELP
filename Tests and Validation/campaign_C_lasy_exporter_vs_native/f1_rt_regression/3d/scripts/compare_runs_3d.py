"""
Campaign C part 3 (3D) -- RUN-LEVEL comparison: EPOCH E_y (z=0 midplane)
from the new rt-exporter injection (../lasy_exporter/) vs test B's
original 3d/lasy run.

Same expectation as the 2D version (see its docstring): the two runs
inject the same physical field into the same deck, so the residual
should sit at the injector floor + the constant CEP piston's imprint.

3D-specific: sliced at z=0, one SDF file at a time (147M cells / 5.9 GB
per snapshot -- bulk-loading all 16 would OOM; same pattern as test B's
own 3d/field_comparison.py).

Must be run via sbatch (SDF reads are login-node-prohibited on Viking).
"""
import sys
import glob
from pathlib import Path

import numpy as np
import xarray as xr
import sdf_xarray as sdfxr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))
import physics_params as P
import epoch_tools.analysis.viking_analysis_lib as val
import epoch_tools.analysis.analytical_field as af

NEW_RUN = HERE.parent / "lasy_exporter"
OLD_RUN = (HERE.parents[3] / "campaign_B_tight_focus_f1"
           / "3d_full_resolution" / "lasy")

results = HERE.parent / "results"
figures = results / "figures"
figures.mkdir(parents=True, exist_ok=True)


def load_z0_timeseries(run_dir):
    """One SDF file at a time -- see test B's 3d/field_comparison.py."""
    files = sorted(glob.glob(str(run_dir / "*.sdf")))
    slices, times = [], []
    x = y = None
    for f in files:
        ds = sdfxr.open_mfdataset([f], separate_times=True,
                                  data_vars=["Electric_Field_Ey"])
        ey = ds["Electric_Field_Ey"]
        time_dim = [d for d in ey.dims if d.startswith("time")][0]
        iz0 = int(np.argmin(np.abs(ey["Z_Grid_mid"].values)))
        slab = ey.isel(Z_Grid_mid=iz0, **{time_dim: 0})
        if x is None:
            x = slab["X_Grid_mid"].values
            y = slab["Y_Grid_mid"].values
        slices.append(np.asarray(slab.values) / af.a0_norm(P.LAMBDA0))
        times.append(float(ey[time_dim].values[0]))
        ds.close()
    order = np.argsort(times)
    data = np.stack([slices[i] for i in order])
    times_arr = np.array(times)[order]
    return xr.DataArray(data, dims=("time", "X_Grid_mid", "Y_Grid_mid"),
                        coords={"time": times_arr,
                                "X_Grid_mid": x, "Y_Grid_mid": y})


new = load_z0_timeseries(NEW_RUN)
old = load_z0_timeseries(OLD_RUN)
t_new = new.coords["time"].values
t_old = old.coords["time"].values
assert len(t_new) == len(t_old), (len(t_new), len(t_old))
assert np.allclose(t_new, t_old, rtol=0, atol=1e-18), "snapshot times differ"

old_aligned = xr.DataArray(old.values, coords=new.coords, dims=new.dims)
diff = new - old_aligned
peak = float(np.abs(old_aligned).max())
env_vals = val.hilbert_envelope_along_x(old_aligned)
times_fs = t_new * 1e15
n_t = new.sizes["time"]

sel = sorted(set([0, n_t // 4, n_t // 2, 3 * n_t // 4, n_t - 1]))
for idx in sel:
    val.save_triptych(old_aligned.isel(time=idx), new.isel(time=idx),
                      diff.isel(time=idx), times_fs[idx],
                      figures / f"run_comparison_3d_t{idx:02d}.png",
                      ref_label="Original pipeline (test B lasy run, z=0)",
                      var_label="New rt-exporter run (z=0)")

rows = []
for idx in range(n_t):
    m = val.difference_metrics(diff.isel(time=idx).values, env_vals[idx])
    rows.append((f"{times_fs[idx]:.2f}", f"{m['max_abs']:.6e}",
                 f"{m['max_abs']/peak*100:.5f}", f"{m['max_rel']*100:.5f}"))
val.write_metrics_csv(results / "run_comparison_3d_metrics.csv", rows,
                      ("time_fs", "max_abs_a0", "max_abs_pct_peak",
                       "max_rel_pct"))

overall = val.difference_metrics(diff.values, env_vals)
summary = f"""Campaign C part 3 (3D) -- run-level comparison (z=0 midplane)
=================================================================
New run : {NEW_RUN}
Old run : {OLD_RUN}
Peak |Ey| (old run, z=0) : {peak:.4f} a0
Snapshots                : {n_t}

Overall metrics (new-exporter run vs original-pipeline run):
  max |difference|         : {overall['max_abs']:.6e} a0
  max |difference| / peak  : {overall['max_abs']/peak*100:.4f} %
  max relative difference  : {overall['max_rel']*100:.4f} %
  mean relative difference : {overall['mean_rel']*100:.4f} %

Expectation: injector floor (~1e-3 % band) + <=0.14 % piston imprint.
"""
val.write_text(results / "run_comparison_3d_summary.txt", summary)
print(summary)
