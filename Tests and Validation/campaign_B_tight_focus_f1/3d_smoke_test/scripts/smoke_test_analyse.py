"""
Test B (3D) smoke-test analysis -- coarse-resolution sanity check for the
paraxial/ and lasy/ generators before committing the full allocation (see
VIKING_PROMPT_urgent_2x2_and_f1_demo.md). Runs the same waist-scan method
as injector_2x2_validation/3d/analyse.py: Hilbert envelope along x at the
snapshot nearest the beam's peak-at-focus time, 2D Gaussian fit per
x-plane, parabolic interpolation to locate the waist minimum. Confirms the
beam actually converges (not diverges) and lands near the theoretical
x_focus/w0 -- this is exactly the kind of check that caught the 30 June
session's sign-transcription bug in the 2D paraxial generator.

Must be run via sbatch (SDF reads are login-node-prohibited on Viking).

Usage: python smoke_test_analyse.py [base_dir]
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

T_PEAK_AT_FOCUS = P.T_CENTRE + P.X_SPOT / P.C_LIGHT


def analyse_cell(name):
    run_dir = BASE / name / "smoke_test"
    ds = lib.load_fields(run_dir, ["Electric_Field_Ey"])
    ey = ds["Electric_Field_Ey"]
    x, y, z = (ey[c].values for c in
              ("X_Grid_mid", "Y_Grid_mid", "Z_Grid_mid"))
    times = ey["time0"].values
    it = int(np.argmin(np.abs(times - T_PEAK_AT_FOCUS)))
    field = ey.isel(time0=it).values
    env = lib.hilbert_envelope_along_x(field)
    scan = lib.waist_scan(env, x, y, z)
    x_focus, w0 = lib.find_waist_minimum(x, scan["w"])
    return dict(x=x, w=scan["w"], x_focus=x_focus, w0=w0, t=float(times[it]))


cells = {name: analyse_cell(name) for name in ("paraxial", "lasy")}

x_theory = np.linspace(0, 2 * P.X_SPOT, 400)
w_theory = P.W0 * np.sqrt(1.0 + ((x_theory - P.X_SPOT) / P.X_R) ** 2)

results, figures = lib.make_results_dirs(BASE)
results, figures = Path(results), Path(figures)

rows = []
for name, c in cells.items():
    rows.append((name, c["x_focus"] * 1e6, c["w0"] * 1e6,
                abs(c["x_focus"] - P.X_SPOT) * 1e6,
                abs(c["w0"] - P.W0) / P.W0))
lib.write_metrics_csv(results / "smoke_test_metrics.csv", rows,
                      ("cell", "x_focus_um", "w0_um", "x_focus_err_um", "w0_rel_err"))

lib.save_waist_comparison(
    x_theory * 1e6,
    [("theory", x_theory * 1e6, w_theory * 1e6, "k-")]
    + [(name, c["x"] * 1e6, c["w"] * 1e6, style) for name, style, c in
       zip(cells, ["C0o-", "C1s--"], cells.values())],
    P.X_SPOT * 1e6, P.W0 * 1e6, figures / "smoke_test_waist_vs_x.png",
    title="Test B (3D) smoke test: paraxial vs LASY, coarse resolution",
)

summary = f"""Test B (3D) smoke test -- coarse-resolution sanity check
=========================================================
Theory: x_focus={P.X_SPOT*1e6:.4f}um  w0={P.W0*1e6:.4f}um  (f/{P.F_NUMBER:.1f}, NA={P.NA:.4f})
Snapshot analysed: nearest t_peak_at_focus={T_PEAK_AT_FOCUS*1e15:.2f}fs

"""
for name, c in cells.items():
    err_x = abs(c["x_focus"] - P.X_SPOT) * 1e6
    err_w = abs(c["w0"] - P.W0) / P.W0 * 100
    converges = "CONVERGES (minimum found)" if np.isfinite(c["x_focus"]) else "NO MINIMUM FOUND -- FAIL"
    summary += (f"{name:10s}: x_focus={c['x_focus']*1e6:8.4f}um (err {err_x:6.4f}um)  "
               f"w0={c['w0']*1e6:7.4f}um (err {err_w:6.2f}%)  t={c['t']*1e15:.2f}fs  {converges}\n")

summary += """
PASS criterion (coarse smoke test, not final precision): both cells show a
genuine w(x) minimum (not monotonic increase/divergence) landing within
~10-20% of theory x_focus/w0 at this coarse resolution.
"""
lib.write_text(Path(results) / "smoke_test_summary.txt", summary)
print(summary)
