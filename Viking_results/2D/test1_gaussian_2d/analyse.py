"""
Test 1 (2D) analysis — Analytical vs Numerical Gaussian beam (profile only)
==========================================================================
Head-less replacement for the laptop notebook. Writes figures and text results
to ``results/`` under this test directory.

Run (after both EPOCH runs have produced .sdf files):
    python analyse.py                 # uses this directory
    python analyse.py /path/to/test   # or an explicit base directory
"""

import argparse
import os
import sys

import numpy as np
import xarray as xr

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import viking_analysis_lib as val  # noqa: E402

LAMBDA0 = 1.0e-6
REF_LABEL, VAR_LABEL = "Analytical", "Numerical"


def main(base_dir):
    results, figures = val.make_results_dirs(base_dir)
    ref = val.load_ey(os.path.join(base_dir, "analytical"), LAMBDA0).load()
    var = val.load_ey(os.path.join(base_dir, "numerical"), LAMBDA0).load()

    # Align on the common time steps (both decks share dt_snapshot)
    n_t = min(ref.sizes["time"], var.sizes["time"])
    ref, var = ref.isel(time=slice(0, n_t)), var.isel(time=slice(0, n_t))

    diff = var - ref
    peak = float(np.abs(ref).max())

    # Reference envelope (Hilbert along x) wrapped back into a DataArray
    env_vals = val.hilbert_envelope_along_x(ref)
    env = xr.DataArray(env_vals, coords=ref.coords, dims=ref.dims)
    env_mask = env_vals > 0.01 * env_vals.max()
    rel = xr.DataArray(np.abs(diff.values) / np.where(env_mask, env_vals, np.nan),
                       coords=ref.coords, dims=ref.dims)

    times_fs = ref.coords["time"].values * 1e15
    sel = sorted(set([0, n_t // 4, n_t // 2, 3 * n_t // 4, n_t - 1]))

    # ---- Figures ----
    for idx in sel:
        t = times_fs[idx]
        val.save_triptych(ref.isel(time=idx), var.isel(time=idx), diff.isel(time=idx),
                          t, os.path.join(figures, f"field_t{idx:02d}.png"),
                          ref_label=REF_LABEL, var_label=VAR_LABEL)
        val.save_relative_diff_panel(env.isel(time=idx), rel.isel(time=idx), t,
                                     os.path.join(figures, f"reldiff_t{idx:02d}.png"),
                                     ref_label=REF_LABEL)
        iy0 = int(np.argmin(np.abs(ref.coords["Y_Grid_mid"].values)))
        xum = ref.coords["X_Grid_mid"].values * 1e6
        val.save_lineout(
            xum,
            [(REF_LABEL, ref.isel(time=idx, Y_Grid_mid=iy0).values, "b-"),
             (VAR_LABEL, var.isel(time=idx, Y_Grid_mid=iy0).values, "r--")],
            t, os.path.join(figures, f"lineout_y0_t{idx:02d}.png"),
            title=f"On-axis (y=0) E_y  |  t = {t:.1f} fs")

    # ---- Per-snapshot metrics ----
    rows = []
    for idx in range(n_t):
        m = val.difference_metrics(diff.isel(time=idx).values, env_vals[idx])
        rows.append((f"{times_fs[idx]:.2f}", f"{m['max_abs']:.6e}",
                     f"{m['max_abs']/peak*100:.5f}", f"{m['max_rel']*100:.5f}"))
    val.write_metrics_csv(os.path.join(results, "metrics.csv"), rows,
                          ("time_fs", "max_abs_a0", "max_abs_pct_peak", "max_rel_pct"))

    overall = val.difference_metrics(diff.values, env_vals)
    summary = f"""Test 1 (2D) - Analytical vs Numerical Gaussian beam (profile only)
==================================================================

Beam: lambda0 = 1 um, pulse_FWHM = 2 um (w0 ~ 1.2 um), focus 10 um downstream
of the x_min boundary. Both runs use the SAME analytical paraxial phase; the
numerical run injects the amplitude envelope from a .dat file. Any difference
therefore measures the fidelity of the custom injection path, not the physics.

Reference  : analytical (native gauss profile)
Comparison : numerical  (.dat amplitude injection)
Peak |Ey|  : {peak:.4f} a0
Snapshots  : {n_t}

Overall metrics (vs analytical reference):
  max |difference|        : {overall['max_abs']:.6e} a0
  max |difference| / peak  : {overall['max_abs']/peak*100:.4f} %
  max relative difference  : {overall['max_rel']*100:.4f} %
  mean relative difference : {overall['mean_rel']*100:.4f} %

Interpretation:
  A small (<~0.1 %) relative difference confirms the custom .dat injector
  reproduces EPOCH's native analytical Gaussian to the discretisation floor.
  The residual concentrates at the beam edges (steepest gradients), where
  bilinear interpolation of the .dat profile is least accurate.

Figures (results/figures/):
  field_t*.png    - analytical | numerical | difference (full 2D Ey)
  reldiff_t*.png  - reference envelope | relative difference
  lineout_y0_t*.png - on-axis Ey overlay + residual
metrics.csv        - per-snapshot max abs / relative difference
"""
    val.write_text(os.path.join(results, "summary.txt"), summary)
    print(summary)
    print(f"Results written to {results}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("base_dir", nargs="?",
                   default=os.path.dirname(os.path.abspath(__file__)))
    main(p.parse_args().base_dir)
