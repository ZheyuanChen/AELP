"""
Test 3 (3D) analysis — Analytical vs Numerical Gaussian beam (profile only)
==========================================================================
Head-less replacement for the 3D notebook. Because the data is 3D, we analyse
the z = 0 mid-plane (the direct analogue of the 2D test) plus a montage of
z-slices at a fixed time, and an on-axis line cut. Slices are loaded lazily to
keep memory bounded.

Outputs (results/):
  figures/field_z0_t*.png    - analytical | numerical | difference at z=0
  figures/reldiff_z0_t*.png  - envelope | relative difference at z=0
  figures/lineout_axis_t*.png - on-axis (y=z=0) Ey overlay + residual
  figures/zmontage.png       - analytical & difference across several z (fixed t)
  summary.txt, metrics.csv

Run after both EPOCH runs have produced .sdf files:
  python analyse.py [base_dir]
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import viking_analysis_lib as val  # noqa: E402

LAMBDA0 = 1.0e-6
REF_LABEL, VAR_LABEL = "Analytical", "Numerical"


def main(base_dir):
    results, figures = val.make_results_dirs(base_dir)
    ref = val.load_ey(os.path.join(base_dir, "analytical"), LAMBDA0)   # lazy
    var = val.load_ey(os.path.join(base_dir, "numerical"), LAMBDA0)    # lazy
    n_t = min(ref.sizes["time"], var.sizes["time"])
    ref, var = ref.isel(time=slice(0, n_t)), var.isel(time=slice(0, n_t))

    # Overall peak and max abs difference via streaming (dask) reductions
    peak = float(np.abs(ref).max())
    max_abs_3d = float(np.abs(var - ref).max())

    # --- z = 0 mid-plane (materialise only this slice) ---
    ref_z0 = ref.sel(Z_Grid_mid=0, method="nearest").load()
    var_z0 = var.sel(Z_Grid_mid=0, method="nearest").load()
    diff_z0 = var_z0 - ref_z0
    env_vals = val.hilbert_envelope_along_x(ref_z0)
    env = xr.DataArray(env_vals, coords=ref_z0.coords, dims=ref_z0.dims)
    env_mask = env_vals > 0.01 * env_vals.max()
    rel = xr.DataArray(np.abs(diff_z0.values) / np.where(env_mask, env_vals, np.nan),
                       coords=ref_z0.coords, dims=ref_z0.dims)

    times_fs = ref_z0.coords["time"].values * 1e15
    sel = sorted(set([0, n_t // 4, n_t // 2, 3 * n_t // 4, n_t - 1]))
    for idx in sel:
        t = times_fs[idx]
        val.save_triptych(ref_z0.isel(time=idx), var_z0.isel(time=idx), diff_z0.isel(time=idx),
                          t, os.path.join(figures, f"field_z0_t{idx:02d}.png"),
                          ref_label=REF_LABEL, var_label=VAR_LABEL)
        val.save_relative_diff_panel(env.isel(time=idx), rel.isel(time=idx), t,
                                     os.path.join(figures, f"reldiff_z0_t{idx:02d}.png"),
                                     ref_label=REF_LABEL)

    # --- on-axis (y = z = 0) line cut ---
    ref_axis = ref.sel(Y_Grid_mid=0, Z_Grid_mid=0, method="nearest").load()
    var_axis = var.sel(Y_Grid_mid=0, Z_Grid_mid=0, method="nearest").load()
    xum = ref_axis.coords["X_Grid_mid"].values * 1e6
    for idx in sel:
        t = times_fs[idx]
        val.save_lineout(xum,
                         [(REF_LABEL, ref_axis.isel(time=idx).values, "b-"),
                          (VAR_LABEL, var_axis.isel(time=idx).values, "r--")],
                         t, os.path.join(figures, f"lineout_axis_t{idx:02d}.png"),
                         title=f"On-axis (y=z=0) E_y  |  t = {t:.1f} fs")

    # --- z-slice montage at the mid-time snapshot (single 3D snapshot in RAM) ---
    mid = n_t // 2
    t_mid = times_fs[mid]
    ref_t = ref.isel(time=mid).load()      # (X, Y, Z)
    var_t = var.isel(time=mid).load()
    diff_t = var_t - ref_t
    z_vals = ref_t.coords["Z_Grid_mid"].values
    z_pick = [int(round(f * (len(z_vals) - 1))) for f in (0.5, 0.6, 0.7, 0.8)]  # z=0 outward
    xu = ref_t.coords["X_Grid_mid"].values * 1e6
    yu = ref_t.coords["Y_Grid_mid"].values * 1e6
    ext = [xu.min(), xu.max(), yu.min(), yu.max()]
    vmax_f = float(np.abs(ref_t).max())
    vmax_d = float(np.abs(diff_t).max()) or 1e-30

    fig, ax = plt.subplots(2, len(z_pick), figsize=(4.2 * len(z_pick), 8))
    for col, iz in enumerate(z_pick):
        a = val.orient_xy(ref_t.isel(Z_Grid_mid=iz)).values
        d = val.orient_xy(diff_t.isel(Z_Grid_mid=iz)).values
        im0 = ax[0, col].imshow(a, origin="lower", aspect="auto", extent=ext,
                                cmap="RdBu_r", vmin=-vmax_f, vmax=vmax_f)
        ax[0, col].set(title=f"Analytical  z={z_vals[iz]*1e6:.1f} µm",
                       xlabel="x (µm)", ylabel="y (µm)"); fig.colorbar(im0, ax=ax[0, col])
        im1 = ax[1, col].imshow(d, origin="lower", aspect="auto", extent=ext,
                                cmap="RdBu_r", vmin=-vmax_d, vmax=vmax_d)
        ax[1, col].set(title=f"Difference  z={z_vals[iz]*1e6:.1f} µm",
                       xlabel="x (µm)", ylabel="y (µm)"); fig.colorbar(im1, ax=ax[1, col])
    fig.suptitle(f"z-slice montage at t = {t_mid:.1f} fs")
    fig.tight_layout(); fig.savefig(os.path.join(figures, "zmontage.png"), dpi=140); plt.close(fig)

    # --- metrics (z=0 plane) ---
    rows = []
    for idx in range(n_t):
        m = val.difference_metrics(diff_z0.isel(time=idx).values, env_vals[idx])
        rows.append((f"{times_fs[idx]:.2f}", f"{m['max_abs']:.6e}", f"{m['max_rel']*100:.5f}"))
    val.write_metrics_csv(os.path.join(results, "metrics.csv"), rows,
                          ("time_fs", "max_abs_a0_z0", "max_rel_pct_z0"))
    overall_z0 = val.difference_metrics(diff_z0.values, env_vals)

    summary = f"""Test 3 (3D) - Analytical vs Numerical Gaussian beam (profile only)
==================================================================

3D counterpart of Test 1. Native gauss(y,0,w0)*gauss(z,0,w0) (w0 = 3 um) vs a
2D (y,z) spatial profile injected from spatial_profile.dat. phase = 0.

Reference  : analytical (native separable gauss)
Comparison : numerical  (.dat 2D spatial injection)
Peak |Ey|  : {peak:.4f} a0
Snapshots  : {n_t}

Full-3D max |difference|      : {max_abs_3d:.6e} a0  ({max_abs_3d/peak*100:.4f} % of peak)
z=0 plane max relative diff   : {overall_z0['max_rel']*100:.4f} %
z=0 plane mean relative diff  : {overall_z0['mean_rel']*100:.4f} %

Interpretation:
  As in 2D, a small relative difference confirms the 3D custom injector
  reproduces the native analytical Gaussian. The z-slice montage shows the
  difference stays small across z (the beam's Gaussian z-structure is captured),
  and the on-axis line cuts overlay pointwise.

Figures in results/figures/: field_z0_t*.png, reldiff_z0_t*.png,
lineout_axis_t*.png, zmontage.png
"""
    val.write_text(os.path.join(results, "summary.txt"), summary)
    print(summary)
    print(f"Results written to {results}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("base_dir", nargs="?",
                   default=os.path.dirname(os.path.abspath(__file__)))
    main(p.parse_args().base_dir)
