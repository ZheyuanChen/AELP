"""
Test B (3D) field comparison -- analytical (closed-form paraxial, dim="3d"
i.e. full Gouy phase + circularly-symmetric amplitude scaling) vs
numerical (EPOCH-simulated) E_y, animated + screenshots, normalised to
a0. Run against BOTH paraxial/ and lasy/ (see the 2D version's docstring
for why -- paraxial should closely track the formula, lasy should
genuinely diverge at this f/1).

3D-specific: sliced at z=0 (the beam's symmetry midplane), then reuses
the 2D library's plotting/animation functions on that 2D slice.

Memory note: same one-file-at-a-time loading as
injector_2x2_validation/3d/field_comparison.py -- SDF's per-file block
storage means even a z=0 slice requires reading the whole 3D block first,
and this script needs all 16 snapshots (this campaign's tight_focus_f1
cells are 147.2M cells/5.9 GB per snapshot, bigger than Campaign A's
88.3M/3.5GB), so bulk-loading the full time series risked a severe OOM.

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
sys.path.insert(0, str(HERE.parent.parent))
sys.path.insert(0, str(HERE.parent))
import viking_analysis_lib as val
import analytical_field as af
import physics_params as P

BASE = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
AMP = 3.2e12  # V/m, matches `amp = 3.2e12` in paraxial/input.deck and lasy/input.deck


def make_results_dirs(base_dir):
    results = base_dir / "results"
    figures = results / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    return results, figures


results, figures = make_results_dirs(BASE)
animations = results / "animations"
animations.mkdir(parents=True, exist_ok=True)


def load_z0_timeseries(run_dir):
    """One SDF file at a time -- see injector_2x2_validation/3d/
    field_comparison.py's docstring for why."""
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
        slices.append(np.asarray(slab.values) / af.a0_norm(P.LAMBDA0))
        times.append(float(ey[time_dim].values[0]))
        ds.close()
    order = np.argsort(times)
    data = np.stack([slices[i] for i in order])
    times_arr = np.array(times)[order]
    return xr.DataArray(data, dims=("time", "X_Grid_mid", "Y_Grid_mid"),
                        coords={"time": times_arr, "X_Grid_mid": x, "Y_Grid_mid": y})


def compare(var_name, var_label):
    print(f"\n--- {var_name} ---")
    var = load_z0_timeseries(BASE / var_name)
    n_t = var.sizes["time"]
    times_fs = var.coords["time"].values * 1e15
    times_s = var.coords["time"].values
    x = var.coords["X_Grid_mid"].values
    y = var.coords["Y_Grid_mid"].values
    x_first = var.dims.index("X_Grid_mid") < var.dims.index("Y_Grid_mid")
    Xg, Yg = np.meshgrid(x, y, indexing="ij") if x_first else np.meshgrid(y, x, indexing="ij")

    ref_vals = np.empty(var.shape)
    for it, t in enumerate(times_s):
        field = af.analytical_paraxial_ey(
            Xg, Yg, np.zeros_like(Xg), t, w0=P.W0, x_r=P.X_R, x_spot=P.X_SPOT,
            tau=P.PULSE_TAU, t_centre=P.T_CENTRE, k0=P.K0, amp=AMP, dim="3d")
        ref_vals[it] = field / af.a0_norm(P.LAMBDA0)
    ref = xr.DataArray(ref_vals, coords=var.coords, dims=var.dims)

    diff = var - ref
    peak = float(np.abs(ref).max())
    env_vals = val.hilbert_envelope_along_x(ref)

    ref_label = "Analytical (closed-form paraxial, z=0)"
    sel = sorted(set([0, n_t // 4, n_t // 2, 3 * n_t // 4, n_t - 1]))
    for idx in sel:
        t = times_fs[idx]
        val.save_triptych(ref.isel(time=idx), var.isel(time=idx), diff.isel(time=idx),
                          t, figures / f"field_comparison_{var_name}_t{idx:02d}.png",
                          ref_label=ref_label, var_label=var_label)

    anim_path = val.animate_triptych(ref, var, diff, times_fs,
                                     animations / f"field_comparison_{var_name}.mp4",
                                     ref_label=ref_label, var_label=var_label)

    rows = []
    for idx in range(n_t):
        m = val.difference_metrics(diff.isel(time=idx).values, env_vals[idx])
        rows.append((f"{times_fs[idx]:.2f}", f"{m['max_abs']:.6e}",
                    f"{m['max_abs']/peak*100:.5f}", f"{m['max_rel']*100:.5f}"))
    val.write_metrics_csv(results / f"field_comparison_{var_name}_metrics.csv", rows,
                          ("time_fs", "max_abs_a0", "max_abs_pct_peak", "max_rel_pct"))

    overall = val.difference_metrics(diff.values, env_vals)
    summary = f"""Test B (3D) -- field comparison: analytical (closed-form paraxial) vs {var_label}
====================================================================================================
Comparison : {var_name} ({var_label}), z=0 midplane slice
Peak |Ey|  : {peak:.4f} a0
Snapshots  : {n_t}

Overall metrics (vs analytical closed-form reference):
  max |difference|         : {overall['max_abs']:.6e} a0
  max |difference| / peak  : {overall['max_abs']/peak*100:.4f} %
  max relative difference  : {overall['max_rel']*100:.4f} %
  mean relative difference : {overall['mean_rel']*100:.4f} %

Animation: {anim_path}
"""
    val.write_text(results / f"field_comparison_{var_name}_summary.txt", summary)
    print(summary)


for cell, label in [("paraxial", "Numerical (paraxial-injected)"),
                    ("lasy", "Numerical (LASY-injected)")]:
    compare(cell, label)
