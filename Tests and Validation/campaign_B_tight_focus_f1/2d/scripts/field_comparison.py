"""
Test B (2D) field comparison -- analytical (closed-form paraxial) vs
numerical (EPOCH-simulated) E_y, animated + screenshots, normalised to
a0. Run against BOTH numerical cells:
  paraxial : expected to closely track the analytical formula (both
      encode the same physics; residual should be discretisation-only).
  lasy     : expected to genuinely diverge from the analytical paraxial
      formula at this f/1 (NA=0.5) -- that divergence IS Test B's
      finding (see analyse.py's RMS numbers), so seeing it appear here
      too, in the raw field rather than just the waist-scan metric,
      is the point.

Unlike Campaign A, there is no native-deck EPOCH cell to use as the
analytical reference here (both paraxial/ and lasy/ inject via file
only), so the reference is a genuine closed-form evaluation
(analytical_field.analytical_paraxial_ey), not another EPOCH run.

Must be run via sbatch (SDF reads are login-node-prohibited on Viking).
Usage: python field_comparison.py [base_dir]
"""
import sys
from pathlib import Path

import numpy as np
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "2D" / "common"))
sys.path.insert(0, str(HERE.parent.parent))
sys.path.insert(0, str(HERE.parent))
import viking_analysis_lib as val
import analytical_field as af
import physics_params as P

BASE = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
AMP = 3.2e12  # V/m, matches `amp = 3.2e12` in paraxial/input.deck and lasy/input.deck

results, figures = val.make_results_dirs(BASE)
results, figures = Path(results), Path(figures)
animations = Path(BASE) / "results" / "animations"
animations.mkdir(parents=True, exist_ok=True)


def compare(var_name, var_label):
    print(f"\n--- {var_name} ---")
    var = val.load_ey(str(BASE / var_name), P.LAMBDA0).load()
    n_t = var.sizes["time"]
    times_fs = var.coords["time"].values * 1e15
    times_s = var.coords["time"].values
    x = var.coords["X_Grid_mid"].values
    y = var.coords["Y_Grid_mid"].values
    Xg, Yg = np.meshgrid(x, y, indexing="ij") if var.dims.index("X_Grid_mid") < var.dims.index("Y_Grid_mid") \
        else np.meshgrid(y, x, indexing="ij")

    ref_vals = np.empty(var.shape)
    for it, t in enumerate(times_s):
        field = af.analytical_paraxial_ey(
            Xg, Yg, None, t, w0=P.W0, x_r=P.X_R, x_spot=P.X_SPOT,
            tau=P.PULSE_TAU, t_centre=P.T_CENTRE, k0=P.K0, amp=AMP, dim="2d")
        ref_vals[it] = field / af.a0_norm(P.LAMBDA0)
    ref = xr.DataArray(ref_vals, coords=var.coords, dims=var.dims)

    diff = var - ref
    peak = float(np.abs(ref).max())
    env_vals = val.hilbert_envelope_along_x(ref)

    ref_label = "Analytical (closed-form paraxial)"
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
    val.write_metrics_csv(Path(results) / f"field_comparison_{var_name}_metrics.csv", rows,
                          ("time_fs", "max_abs_a0", "max_abs_pct_peak", "max_rel_pct"))

    overall = val.difference_metrics(diff.values, env_vals)
    summary = f"""Test B (2D) -- field comparison: analytical (closed-form paraxial) vs {var_label}
====================================================================================================
Comparison : {var_name} ({var_label})
Peak |Ey|  : {peak:.4f} a0
Snapshots  : {n_t}

Overall metrics (vs analytical closed-form reference):
  max |difference|         : {overall['max_abs']:.6e} a0
  max |difference| / peak  : {overall['max_abs']/peak*100:.4f} %
  max relative difference  : {overall['max_rel']*100:.4f} %
  mean relative difference : {overall['mean_rel']*100:.4f} %

Animation: {anim_path}
"""
    val.write_text(Path(results) / f"field_comparison_{var_name}_summary.txt", summary)
    print(summary)


for cell, label in [("paraxial", "Numerical (paraxial-injected)"),
                    ("lasy", "Numerical (LASY-injected)")]:
    compare(cell, label)
