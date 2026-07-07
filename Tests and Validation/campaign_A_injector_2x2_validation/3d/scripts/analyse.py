"""
Test A (3D) analysis -- 3-cell injector isolation chain.
  amp_deck_phase_deck  -> baseline (ground truth)
  amp_file_phase_deck  -> isolates amplitude-injector error
  amp_file_phase_file  -> total (amp+phase)-injector error

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
results, figures = lib.make_results_dirs(BASE)
results, figures = Path(results), Path(figures)


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
    it_lo = int(np.clip(np.searchsorted(times, P.T_PEAK_AT_FOCUS) - 1, 0, len(times) - 2))
    it_hi = it_lo + 1
    t_lo, t_hi = float(times[it_lo]), float(times[it_hi])
    frac = (P.T_PEAK_AT_FOCUS - t_lo) / (t_hi - t_lo)

    env_lo = lib.hilbert_envelope_along_x(ey.isel(time0=it_lo).values)
    env_hi = lib.hilbert_envelope_along_x(ey.isel(time0=it_hi).values)
    scan_lo = lib.waist_scan(env_lo, x, y, z)
    scan_hi = lib.waist_scan(env_hi, x, y, z)
    w_interp = (1.0 - frac) * scan_lo["w"] + frac * scan_hi["w"]
    x_focus, w0 = lib.find_waist_minimum(x, w_interp)

    # Nearest-snapshot result too, for comparison (the previous method's
    # number, kept so the size of the quantisation correction is visible).
    nearest_is_lo = abs(t_lo - P.T_PEAK_AT_FOCUS) <= abs(t_hi - P.T_PEAK_AT_FOCUS)
    t_nearest = t_lo if nearest_is_lo else t_hi
    x_focus_nearest, w0_nearest = lib.find_waist_minimum(
        x, scan_lo["w"] if nearest_is_lo else scan_hi["w"])

    return dict(x=x, w=w_interp, x_focus=x_focus, w0=w0, t=P.T_PEAK_AT_FOCUS,
               t_lo=t_lo, t_hi=t_hi, frac=frac,
               x_focus_nearest=x_focus_nearest, w0_nearest=w0_nearest,
               t_nearest=t_nearest)


cells = {name: analyse_cell(name) for name in
        ("amp_deck_phase_deck", "amp_file_phase_deck", "amp_file_phase_file")}

x_theory = np.linspace(0, 2 * P.X_SPOT, 400)
w_theory = P.w_of_xi(x_theory - P.X_SPOT)

rows = []
for name, c in cells.items():
    rows.append((name, c["x_focus"] * 1e6, c["w0"] * 1e6,
                abs(c["x_focus"] - P.X_SPOT) * 1e6,
                abs(c["w0"] - P.W0) / P.W0,
                c["x_focus_nearest"] * 1e6, c["w0_nearest"] * 1e6))
lib.write_metrics_csv(Path(results) / "metrics.csv", rows,
                      ("cell", "x_focus_um", "w0_um", "x_focus_err_um", "w0_rel_err",
                       "x_focus_nearest_snapshot_um", "w0_nearest_snapshot_um"))

lib.save_waist_comparison(
    x_theory * 1e6,
    [("theory", x_theory * 1e6, w_theory * 1e6, "k-")]
    + [(name, c["x"] * 1e6, c["w"] * 1e6, style) for name, style, c in
       zip(cells, ["C0o-", "C1s--", "C2^:"], cells.values())],
    P.X_SPOT * 1e6, P.W0 * 1e6, figures / "waist_vs_x.png",
    title="Test A (3D): injector isolation chain",
)

base = cells["amp_deck_phase_deck"]
amp_only = cells["amp_file_phase_deck"]
both = cells["amp_file_phase_file"]
valid = np.isfinite(base["w"]) & np.isfinite(amp_only["w"]) & np.isfinite(both["w"])
amp_injector_rms = float(np.sqrt(np.nanmean(
    ((amp_only["w"][valid] - base["w"][valid]) / base["w"][valid]) ** 2)))
total_injector_rms = float(np.sqrt(np.nanmean(
    ((both["w"][valid] - base["w"][valid]) / base["w"][valid]) ** 2)))
phase_injector_rms = float(np.sqrt(np.nanmean(
    ((both["w"][valid] - amp_only["w"][valid]) / amp_only["w"][valid]) ** 2)))

# waist_vs_x.png compares each cell against the IDEALISED THEORY curve, which
# every cell (including the native-deck baseline) sits ~0.4um off of -- a
# finite-resolution/numerical-dispersion effect common to the whole chain,
# not an injector-pipeline error. That's not what this test actually checks
# for (injector fidelity = cells matching EACH OTHER, not matching theory
# exactly, which no finite-grid PIC run does). This second plot makes the
# real pass/fail metric visible directly: relative w(x) difference between
# each file-injected cell and the native-deck baseline, which the RMS
# numbers below show is ~1e-4 to 1e-3 % -- i.e. the curves overlap almost
# exactly once you're looking at the right quantity.
fig, ax = plt.subplots(figsize=(9, 4.5))
x_um = base["x"][valid] * 1e6
ax.plot(x_um, (amp_only["w"][valid] - base["w"][valid]) / base["w"][valid] * 100,
       "C1-", lw=1.0, label="amp_file_phase_deck vs baseline")
ax.plot(x_um, (both["w"][valid] - base["w"][valid]) / base["w"][valid] * 100,
       "C2-", lw=1.0, label="amp_file_phase_file vs baseline")
ax.axhline(0, color="grey", ls=":", lw=1)
ax.set(xlabel="x (um)", ylabel="relative w(x) difference vs baseline (%)",
      title="Test A (3D): injector-pipeline error (the actual pass/fail metric)")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(figures / "injector_error_vs_x.png", dpi=140)
plt.close(fig)

any_cell = next(iter(cells.values()))
summary = f"""Test A (3D) -- injector 2x2 (3-cell) isolation
================================================
Beam: lambda0={P.LAMBDA0*1e6:.3f}um w0={P.W0*1e6:.3f}um NA~{P.NA:.4f} (moderate;
validates the pipeline, not paraxial theory).

Snapshot timing: ideal t_peak@focus = {P.T_PEAK_AT_FOCUS*1e15:.3f} fs, but
output only exists every dt_snapshot -- bracketing snapshots at
{any_cell['t_lo']*1e15:.1f} fs and {any_cell['t_hi']*1e15:.1f} fs
(interpolation fraction {any_cell['frac']:.3f}). x_focus/w0 below use the
fitted w(x) curve LINEARLY INTERPOLATED between those two snapshots to
the exact ideal time (not the raw field -- see analyse_cell's comment).
Old method (snap to nearest snapshot, t={any_cell['t_nearest']*1e15:.1f} fs)
numbers included alongside for comparison.

x_focus / w0 per cell (theory: x_focus={P.X_SPOT*1e6:.4f}um, w0={P.W0*1e6:.4f}um):
"""
for name, c in cells.items():
    summary += (f"  {name:22s}: x_focus={c['x_focus']*1e6:10.6f}um  "
               f"w0={c['w0']*1e6:9.6f}um   (nearest-snapshot was: "
               f"x_focus={c['x_focus_nearest']*1e6:10.6f}um  "
               f"w0={c['w0_nearest']*1e6:9.6f}um)\n")

summary += f"""
w(x) RMS relative error, over the interior valid range:
  amplitude-injector-only (amp_file_phase_deck vs baseline)   = {amp_injector_rms*100:.3e}%
  phase-injector-additional (amp_file_phase_file vs amp-only) = {phase_injector_rms*100:.3e}%
  total (amp_file_phase_file vs baseline)                     = {total_injector_rms*100:.3e}%

Note: amp_deck_phase_file (native amplitude + file phase) is architecturally
unreachable with the current epoch3d code -- see amp_deck_phase_file/README.md.
The 3-cell chain above still isolates amplitude-injector error and phase-
injector error, just via a chain rather than a symmetric 2x2 grid.
"""
lib.write_text(Path(results) / "summary.txt", summary)
print(summary)
