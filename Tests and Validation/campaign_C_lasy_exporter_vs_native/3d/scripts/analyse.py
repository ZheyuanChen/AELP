"""
Campaign C (3D) analysis -- lasy-exporter cell vs native-injector cell.

Adapted from campaign A's 3d/scripts/analyse.py (same waist-scan library,
same snapshot-timing interpolation fix), reduced to the two campaign C
cells (native, lasy_file). See the 2D twin and ../../physics_params.py
for the interpretation: expected residual is the ~1 % lasy-vs-paraxial
physics difference, not campaign A's ~1e-3 % injector floor.

Must be run via sbatch on Viking (SDF reads are login-node-prohibited).
Usage: python analyse.py [base_dir]
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent.parent / "shared_libraries"))
sys.path.insert(0, str(HERE.parent.parent))
import viking_analysis_lib_3d as lib
import physics_params as P

BASE = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parent
results, figures = lib.make_results_dirs(BASE)
results, figures = Path(results), Path(figures)

CELLS = ("native", "lasy_file")


def analyse_cell(name):
    ds = lib.load_fields(BASE / name, ["Electric_Field_Ey"])
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
    it_lo = int(np.clip(np.searchsorted(times, P.T_PEAK_AT_FOCUS) - 1,
                        0, len(times) - 2))
    it_hi = it_lo + 1
    t_lo, t_hi = float(times[it_lo]), float(times[it_hi])
    frac = (P.T_PEAK_AT_FOCUS - t_lo) / (t_hi - t_lo)

    env_lo = lib.hilbert_envelope_along_x(ey.isel(time0=it_lo).values)
    env_hi = lib.hilbert_envelope_along_x(ey.isel(time0=it_hi).values)
    scan_lo = lib.waist_scan(env_lo, x, y, z)
    scan_hi = lib.waist_scan(env_hi, x, y, z)
    w_interp = (1.0 - frac) * scan_lo["w"] + frac * scan_hi["w"]
    x_focus, w0 = lib.find_waist_minimum(x, w_interp)

    nearest_is_lo = abs(t_lo - P.T_PEAK_AT_FOCUS) <= abs(t_hi - P.T_PEAK_AT_FOCUS)
    t_nearest = t_lo if nearest_is_lo else t_hi
    x_focus_nearest, w0_nearest = lib.find_waist_minimum(
        x, scan_lo["w"] if nearest_is_lo else scan_hi["w"])

    return dict(x=x, w=w_interp, x_focus=x_focus, w0=w0,
                t_lo=t_lo, t_hi=t_hi, frac=frac,
                x_focus_nearest=x_focus_nearest, w0_nearest=w0_nearest,
                t_nearest=t_nearest)


cells = {name: analyse_cell(name) for name in CELLS}

x_theory = np.linspace(0, 2 * P.X_SPOT, 400)
w_theory = P.w_of_xi(x_theory - P.X_SPOT)

rows = []
for name, c in cells.items():
    rows.append((name, c["x_focus"] * 1e6, c["w0"] * 1e6,
                 abs(c["x_focus"] - P.X_SPOT) * 1e6,
                 abs(c["w0"] - P.W0) / P.W0,
                 c["x_focus_nearest"] * 1e6, c["w0_nearest"] * 1e6))
lib.write_metrics_csv(Path(results) / "metrics.csv", rows,
                      ("cell", "x_focus_um", "w0_um", "x_focus_err_um",
                       "w0_rel_err", "x_focus_nearest_snapshot_um",
                       "w0_nearest_snapshot_um"))

lib.save_waist_comparison(
    x_theory * 1e6,
    [("theory", x_theory * 1e6, w_theory * 1e6, "k-")]
    + [(name, c["x"] * 1e6, c["w"] * 1e6, style) for name, style, c in
       zip(cells, ["C0o-", "C1s--"], cells.values())],
    P.X_SPOT * 1e6, P.W0 * 1e6, figures / "waist_vs_x.png",
    title="Campaign C (3D): lasy exporter vs native injector",
)

base, test = cells["native"], cells["lasy_file"]
valid = np.isfinite(base["w"]) & np.isfinite(test["w"])
rel = (test["w"][valid] - base["w"][valid]) / base["w"][valid]
w_rms = float(np.sqrt(np.mean(rel ** 2)))

# The pass/fail metric is the cells matching EACH OTHER (finite-grid
# effects offset the whole chain from theory equally); expected scale is
# the ~1 % lasy-vs-paraxial physics gap, amplitude-dominated -- see the
# 2D twin's comment for the failure signatures that would instead
# indicate an exporter bug.
fig, ax = plt.subplots(figsize=(9, 4.5))
ax.plot(base["x"][valid] * 1e6, rel * 100, "C1-", lw=1.0,
        label="lasy_file vs native")
ax.axhline(0, color="grey", ls=":", lw=1)
ax.set(xlabel="x (um)", ylabel="relative w(x) difference vs native (%)",
       title="Campaign C (3D): exporter+physics residual "
             "(expect ~1 %, smooth; NOT a pure injector floor)")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(figures / "exporter_error_vs_x.png", dpi=140)
plt.close(fig)

any_cell = next(iter(cells.values()))
summary = f"""Campaign C (3D) -- lasy-fork exporter vs native injector
=========================================================
Beam: lambda0={P.LAMBDA0*1e6:.3f}um w0={P.W0*1e6:.3f}um NA~{P.NA:.4f}
(campaign A's beam; circularly-symmetric r_yz, FULL Gouy phase).

Snapshot timing: ideal t_peak@focus = {P.T_PEAK_AT_FOCUS*1e15:.3f} fs;
bracketing snapshots at {any_cell['t_lo']*1e15:.1f} fs and
{any_cell['t_hi']*1e15:.1f} fs (interpolation fraction
{any_cell['frac']:.3f}); x_focus/w0 use the fitted w(x) interpolated to
the exact ideal time. Nearest-snapshot numbers alongside for comparison.

x_focus / w0 per cell (theory: x_focus={P.X_SPOT*1e6:.4f}um, w0={P.W0*1e6:.4f}um):
"""
for name, c in cells.items():
    summary += (f"  {name:10s}: x_focus={c['x_focus']*1e6:10.6f}um  "
                f"w0={c['w0']*1e6:9.6f}um   (nearest-snapshot: "
                f"x_focus={c['x_focus_nearest']*1e6:10.6f}um  "
                f"w0={c['w0_nearest']*1e6:9.6f}um)\n")
summary += f"""
w(x) RMS relative difference, lasy_file vs native = {w_rms*100:.3f}%
Expected scale: ~1 % (lasy angular-spectrum vs paraxial at NA~0.265).
For a pointwise field comparison, adapt campaign A's
3d/scripts/field_comparison.py with REF/VAR set to native/lasy_file --
the utilities are shared, only the cell names and expectation text
change (see the 2D twin for the calibrated failure signatures).
"""
lib.write_text(Path(results) / "summary.txt", summary)
print(summary)
