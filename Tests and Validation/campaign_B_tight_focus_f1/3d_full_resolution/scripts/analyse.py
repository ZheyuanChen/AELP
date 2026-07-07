"""
Test B (3D) analysis -- full-resolution paraxial vs LASY, f/1 tight-
focusing demo. Same waist-scan method as smoke_test_analyse.py (Hilbert
envelope along x at the snapshot nearest the beam's peak-at-focus time,
2D Gaussian fit per x-plane, parabolic interpolation to locate the waist
minimum), now on the full 326x672x672 (147.2M cell) runs in paraxial/ and
lasy/ directly, not the coarse smoke_test/ subdirs.

Must be run via sbatch (SDF reads are login-node-prohibited on Viking).

Usage: python analyse.py [base_dir]
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "3D" / "common"))
sys.path.insert(0, str(HERE.parent))
import viking_analysis_lib_3d as lib
import physics_params as P

BASE = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE

T_PEAK_AT_FOCUS = P.T_CENTRE + P.X_SPOT / P.C_LIGHT


def analyse_cell(name):
    run_dir = BASE / name
    ds = lib.load_fields(run_dir, ["Electric_Field_Ey"])
    ey = ds["Electric_Field_Ey"]
    x, y, z = (ey[c].values for c in
              ("X_Grid_mid", "Y_Grid_mid", "Z_Grid_mid"))
    times = ey["time0"].values

    # Bracket the ideal focal time between the two nearest AVAILABLE
    # snapshots and linearly interpolate the FITTED w(x) curve to that
    # exact time, instead of snapping to whichever snapshot happens to be
    # nearest. Interpolating the raw field itself would be unsound -- the
    # carrier oscillates several times per dt_snapshot, so linear
    # interpolation between distant-in-phase snapshots wouldn't
    # reconstruct the intermediate field correctly. The fitted width
    # w(x), by contrast, varies slowly over the pulse duration and
    # interpolates cleanly.
    it_lo = int(np.clip(np.searchsorted(times, T_PEAK_AT_FOCUS) - 1, 0, len(times) - 2))
    it_hi = it_lo + 1
    t_lo, t_hi = float(times[it_lo]), float(times[it_hi])
    frac = (T_PEAK_AT_FOCUS - t_lo) / (t_hi - t_lo)

    env_lo = lib.hilbert_envelope_along_x(ey.isel(time0=it_lo).values)
    env_hi = lib.hilbert_envelope_along_x(ey.isel(time0=it_hi).values)
    scan_lo = lib.waist_scan(env_lo, x, y, z)
    scan_hi = lib.waist_scan(env_hi, x, y, z)
    w_interp = (1.0 - frac) * scan_lo["w"] + frac * scan_hi["w"]
    x_focus, w0 = lib.find_waist_minimum(x, w_interp)

    # Nearest-snapshot result too, for comparison (the previous method's
    # number, kept so the size of the quantisation correction is visible).
    nearest_is_lo = abs(t_lo - T_PEAK_AT_FOCUS) <= abs(t_hi - T_PEAK_AT_FOCUS)
    t_nearest = t_lo if nearest_is_lo else t_hi
    x_focus_nearest, w0_nearest = lib.find_waist_minimum(
        x, scan_lo["w"] if nearest_is_lo else scan_hi["w"])

    return dict(x=x, w=w_interp, x_focus=x_focus, w0=w0, t=T_PEAK_AT_FOCUS,
               t_lo=t_lo, t_hi=t_hi, frac=frac,
               x_focus_nearest=x_focus_nearest, w0_nearest=w0_nearest,
               t_nearest=t_nearest)


cells = {name: analyse_cell(name) for name in ("paraxial", "lasy")}

x_theory = np.linspace(0, 2 * P.X_SPOT, 400)
w_theory = P.W0 * np.sqrt(1.0 + ((x_theory - P.X_SPOT) / P.X_R) ** 2)

results, figures = lib.make_results_dirs(BASE)
results, figures = Path(results), Path(figures)

rows = []
for name, c in cells.items():
    rows.append((name, c["x_focus"] * 1e6, c["w0"] * 1e6,
                abs(c["x_focus"] - P.X_SPOT) * 1e6,
                abs(c["w0"] - P.W0) / P.W0,
                c["x_focus_nearest"] * 1e6, c["w0_nearest"] * 1e6))
lib.write_metrics_csv(results / "metrics.csv", rows,
                      ("cell", "x_focus_um", "w0_um", "x_focus_err_um", "w0_rel_err",
                       "x_focus_nearest_snapshot_um", "w0_nearest_snapshot_um"))

lib.save_waist_comparison(
    x_theory * 1e6,
    [("theory", x_theory * 1e6, w_theory * 1e6, "k-")]
    + [(name, c["x"] * 1e6, c["w"] * 1e6, style) for name, style, c in
       zip(cells, ["C0o-", "C1s--"], cells.values())],
    P.X_SPOT * 1e6, P.W0 * 1e6, figures / "waist_vs_x.png",
    title="Test B (3D): paraxial vs LASY, full resolution",
)

paraxial_vs_lasy_valid = np.isfinite(cells["paraxial"]["w"]) & np.isfinite(cells["lasy"]["w"])
paraxial_vs_lasy_rms = float(np.sqrt(np.nanmean(
    ((cells["paraxial"]["w"][paraxial_vs_lasy_valid] - cells["lasy"]["w"][paraxial_vs_lasy_valid])
     / cells["lasy"]["w"][paraxial_vs_lasy_valid]) ** 2)))

# The RMS number alone hides WHERE along x the two models diverge -- this
# plot shows the relative difference directly as a function of x, which is
# the actual finding this test is after (unlike Campaign A, divergence from
# theory/each-other here is the EXPECTED result, not a bug signature).
fig, ax = plt.subplots(figsize=(9, 4.5))
x_um = cells["paraxial"]["x"][paraxial_vs_lasy_valid] * 1e6
diff_pct = ((cells["paraxial"]["w"][paraxial_vs_lasy_valid] - cells["lasy"]["w"][paraxial_vs_lasy_valid])
           / cells["lasy"]["w"][paraxial_vs_lasy_valid] * 100)
ax.plot(x_um, diff_pct, "C3-", lw=1.0, label="paraxial vs LASY")
ax.axhline(0, color="grey", ls=":", lw=1)
ax.set(xlabel="x (um)", ylabel="relative w(x) difference, paraxial vs LASY (%)",
      title="Test B (3D): paraxial-approximation breakdown vs x (full resolution)")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(figures / "paraxial_vs_lasy_diff_vs_x.png", dpi=140)
plt.close(fig)

any_cell = next(iter(cells.values()))
summary = f"""Test B (3D) -- LASY vs paraxial closed-form, f/1 tight-focusing demo (full resolution)
==========================================================================================
Theory: x_focus={P.X_SPOT*1e6:.4f}um  w0={P.W0*1e6:.4f}um  (f/{P.F_NUMBER:.1f}, NA={P.NA:.4f})

Snapshot timing: ideal t_peak@focus = {T_PEAK_AT_FOCUS*1e15:.3f} fs, but
output only exists every dt_snapshot -- bracketing snapshots at
{any_cell['t_lo']*1e15:.1f} fs and {any_cell['t_hi']*1e15:.1f} fs
(interpolation fraction {any_cell['frac']:.3f}). x_focus/w0 below use the
fitted w(x) curve LINEARLY INTERPOLATED between those two snapshots to
the exact ideal time (not the raw field -- see analyse_cell's comment).
Old method (snap to nearest snapshot, t={any_cell['t_nearest']*1e15:.1f} fs)
numbers included alongside for comparison.

"""
for name, c in cells.items():
    err_x = abs(c["x_focus"] - P.X_SPOT) * 1e6
    err_w = abs(c["w0"] - P.W0) / P.W0 * 100
    err_x_nearest = abs(c["x_focus_nearest"] - P.X_SPOT) * 1e6
    err_w_nearest = abs(c["w0_nearest"] - P.W0) / P.W0 * 100
    summary += (f"{name:10s}: x_focus={c['x_focus']*1e6:8.4f}um (vs theory {err_x:+.4f}um)  "
               f"w0={c['w0']*1e6:7.4f}um (vs theory {err_w:+.2f}%)\n"
               f"{'':10s}  nearest-snapshot was: x_focus={c['x_focus_nearest']*1e6:8.4f}um "
               f"(vs theory {err_x_nearest:+.4f}um)  w0={c['w0_nearest']*1e6:7.4f}um "
               f"(vs theory {err_w_nearest:+.2f}%)\n")

summary += f"""
w(x) RMS relative difference, paraxial vs LASY (over the interior valid
range) = {paraxial_vs_lasy_rms*100:.3f}%
This is the key Test B number: how much the closed-form paraxial formula
diverges from LASY's genuine non-paraxial propagation at NA=0.5, at full
(not coarse smoke-test) resolution.
"""
lib.write_text(results / "summary.txt", summary)
print(summary)
