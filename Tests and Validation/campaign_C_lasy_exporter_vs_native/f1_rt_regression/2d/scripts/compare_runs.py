"""
Campaign C part 3 (2D) -- RUN-LEVEL comparison: EPOCH E_y from the new
rt-exporter injection (../lasy_exporter/) vs test B's original LASY run
(`campaign_B_tight_focus_f1/2d/lasy/`).

Both runs inject the SAME physical field (file-level crosscheck:
amplitude bit-identical, phase 1e-14-identical where amp>1%, constant
-1.4 mrad CEP piston) into the SAME deck/grid, and the injector code is
unchanged between the two epoch_dev binaries (checked: only QED-channel
commits separate them). So the expected residual is the injector/
discretisation floor plus the 1.4 mrad piston's field imprint (~0.14%
of peak bound) -- i.e. campaign A's ~1e-3-of-peak band, NOT parts 1-2's
1-2% physics-difference band. Anything larger flags an exporter bug.

Must be run via sbatch (SDF reads are login-node-prohibited on Viking).
"""
import sys
from pathlib import Path

import numpy as np
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))
import physics_params as P
import epoch_tools.analysis.viking_analysis_lib as val

NEW_RUN = HERE.parent / "lasy_exporter"
OLD_RUN = (HERE.parents[3] / "campaign_B_tight_focus_f1" / "2d" / "lasy")

results = HERE.parent / "results"
figures = results / "figures"
figures.mkdir(parents=True, exist_ok=True)

new = val.load_ey(str(NEW_RUN), P.LAMBDA0).load()
old = val.load_ey(str(OLD_RUN), P.LAMBDA0).load()

# identical decks -> identical grids and snapshot times; align by index
# and assert the times really do match before differencing.
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
                      figures / f"run_comparison_t{idx:02d}.png",
                      ref_label="Original pipeline (test B lasy run)",
                      var_label="New rt-exporter run")

rows = []
for idx in range(n_t):
    m = val.difference_metrics(diff.isel(time=idx).values, env_vals[idx])
    rows.append((f"{times_fs[idx]:.2f}", f"{m['max_abs']:.6e}",
                 f"{m['max_abs']/peak*100:.5f}", f"{m['max_rel']*100:.5f}"))
val.write_metrics_csv(results / "run_comparison_metrics.csv", rows,
                      ("time_fs", "max_abs_a0", "max_abs_pct_peak",
                       "max_rel_pct"))

overall = val.difference_metrics(diff.values, env_vals)
summary = f"""Campaign C part 3 (2D) -- run-level comparison
=================================================================
New run : {NEW_RUN}
Old run : {OLD_RUN}
Peak |Ey| (old run) : {peak:.4f} a0
Snapshots           : {n_t}

Overall metrics (new-exporter run vs original-pipeline run):
  max |difference|         : {overall['max_abs']:.6e} a0
  max |difference| / peak  : {overall['max_abs']/peak*100:.4f} %
  max relative difference  : {overall['max_rel']*100:.4f} %
  mean relative difference : {overall['mean_rel']*100:.4f} %

Expectation: injector floor (~1e-3 % band) + <=0.14 % piston imprint.
"""
val.write_text(results / "run_comparison_summary.txt", summary)
print(summary)
