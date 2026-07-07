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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "3D" / "common"))
sys.path.insert(0, str(HERE.parent))
import viking_analysis_lib_3d as lib
import physics_params as P

BASE = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
results, figures = lib.make_results_dirs(BASE)


def analyse_cell(name):
    ds = lib.load_fields(BASE / name, ["Electric_Field_Ey"])
    ey = ds["Electric_Field_Ey"]
    x, y, z = (ey[c].values for c in
              ("X_Grid_mid", "Y_Grid_mid", "Z_Grid_mid"))
    times = ey["time0"].values
    it = int(np.argmin(np.abs(times - P.T_PEAK_AT_FOCUS)))
    field = ey.isel(time0=it).values
    env = lib.hilbert_envelope_along_x(field)
    scan = lib.waist_scan(env, x, y, z)
    x_focus, w0 = lib.find_waist_minimum(x, scan["w"])
    return dict(x=x, w=scan["w"], x_focus=x_focus, w0=w0, t=float(times[it]))


cells = {name: analyse_cell(name) for name in
        ("amp_deck_phase_deck", "amp_file_phase_deck", "amp_file_phase_file")}

x_theory = np.linspace(0, 2 * P.X_SPOT, 400)
w_theory = P.w_of_xi(x_theory - P.X_SPOT)

lib.save_waist_comparison(
    x_theory * 1e6,
    [("theory", x_theory * 1e6, w_theory * 1e6, "k-")]
    + [(name, c["x"] * 1e6, c["w"] * 1e6, style) for name, style, c in
       zip(cells, ["C0o-", "C1s--", "C2^:"], cells.values())],
    P.X_SPOT * 1e6, P.W0 * 1e6, figures / "waist_vs_x.png",
    title="Test A (3D): injector isolation chain",
)

rows = []
for name, c in cells.items():
    rows.append((name, c["x_focus"] * 1e6, c["w0"] * 1e6,
                abs(c["x_focus"] - P.X_SPOT) * 1e6,
                abs(c["w0"] - P.W0) / P.W0))
lib.write_metrics_csv(Path(results) / "metrics.csv", rows,
                      ("cell", "x_focus_um", "w0_um", "x_focus_err_um", "w0_rel_err"))

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

summary = f"""Test A (3D) -- injector 2x2 (3-cell) isolation
================================================
Beam: lambda0={P.LAMBDA0*1e6:.3f}um w0={P.W0*1e6:.3f}um NA~{P.NA:.4f} (moderate;
validates the pipeline, not paraxial theory).

x_focus / w0 per cell (theory: x_focus={P.X_SPOT*1e6:.4f}um, w0={P.W0*1e6:.4f}um):
"""
for name, c in cells.items():
    summary += (f"  {name:22s}: x_focus={c['x_focus']*1e6:8.4f}um  "
               f"w0={c['w0']*1e6:7.4f}um\n")

summary += f"""
w(x) RMS relative error, over the interior valid range:
  amplitude-injector-only (amp_file_phase_deck vs baseline)   = {amp_injector_rms*100:.3f}%
  phase-injector-additional (amp_file_phase_file vs amp-only) = {phase_injector_rms*100:.3f}%
  total (amp_file_phase_file vs baseline)                     = {total_injector_rms*100:.3f}%

Note: amp_deck_phase_file (native amplitude + file phase) is architecturally
unreachable with the current epoch3d code -- see amp_deck_phase_file/README.md.
The 3-cell chain above still isolates amplitude-injector error and phase-
injector error, just via a chain rather than a symmetric 2x2 grid.
"""
lib.write_text(Path(results) / "summary.txt", summary)
print(summary)
