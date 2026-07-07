"""
Test 1 analysis -- converging Gaussian beam, 3D field-level benchmark.

Compares three things:
  1. numerical/  -- amplitude+phase injected from binary files
                     (use_custom_profile, static path, use_phase_from_file)
  2. analytical/ -- EPOCH's native deck-expression profile/phase
                     (gauss(r_yz,...) + closed-form curvature/Gouy phase)
  3. theory      -- the closed-form paraxial Gaussian beam formula
                     (physics_params.w_theory)

Decomposition logic: if numerical deviates from analytical, that is an
injector/pipeline bug (interpolation, binary read, axis handling). If BOTH
deviate similarly from theory, that is paraxial-approximation error, not a
code bug (expected to be small here -- NA~0.265 is moderate, chosen
specifically to keep the beam within paraxial validity so this test can
cleanly validate the pipeline rather than the approximation).

Validated quantities: w(x) scan, x_focus, w0, and the temporal-envelope-
corrected transverse power-conservation ratio (see
common/viking_analysis_lib_3d.expected_envelope_shape for why raw P(x)
is not expected to be flat for a pulsed, non-CW beam).

Usage: python analyse.py [base_dir]   (default: this script's directory)
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(HERE))
import viking_analysis_lib_3d as lib
import physics_params as P

BASE = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
results, figures = lib.make_results_dirs(BASE)


def analyse_run(run_dir, label):
    ds = lib.load_fields(run_dir, ["Electric_Field_Ey", "Electric_Field_Ez",
                                   "Magnetic_Field_By", "Magnetic_Field_Bz"])
    ey = ds["Electric_Field_Ey"]
    x = ey["X_Grid_mid"].values
    y = ey["Y_Grid_mid"].values
    z = ey["Z_Grid_mid"].values
    times = ey["time0"].values

    it = int(np.argmin(np.abs(times - P.T_PEAK_AT_FOCUS)))
    t_snap = float(times[it])
    print(f"[{label}] using snapshot t={t_snap/1e-15:.2f} fs "
          f"(target {P.T_PEAK_AT_FOCUS/1e-15:.2f} fs, index {it}/{len(times)})")

    ey_f = ey.isel(time0=it).values
    ez_f = ds["Electric_Field_Ez"].isel(time0=it).values
    by_f = ds["Magnetic_Field_By"].isel(time0=it).values
    bz_f = ds["Magnetic_Field_Bz"].isel(time0=it).values

    env = lib.hilbert_envelope_along_x(ey_f)
    scan = lib.waist_scan(env, x, y, z)
    x_focus, w0 = lib.find_waist_minimum(x, scan["w"])

    P_env = lib.transverse_intensity_power(env, y, z)
    ratio = lib.energy_conservation_ratio(P_env, x, t_snap, P.T_CENTRE, P.TAU)

    P_poynting = lib.poynting_x_power(ey_f, ez_f, by_f, bz_f, y, z)

    return dict(x=x, y=y, z=z, t_snap=t_snap, scan=scan,
                x_focus=x_focus, w0=w0, P_env=P_env, ratio=ratio,
                P_poynting=P_poynting,
                ey_slice=ey_f[:, :, len(z) // 2])


num = analyse_run(BASE / "numerical", "numerical")
ana = analyse_run(BASE / "analytical", "analytical")

x_theory_um = np.linspace(0, P.BOX_X, 400) * 1e6
w_theory_um = P.w_theory(np.linspace(0, P.BOX_X, 400) - P.X_SPOT) * 1e6

lib.save_waist_comparison(
    x_theory_um,
    [
        ("theory", x_theory_um, w_theory_um, "k-"),
        ("numerical (file-injected)", num["x"] * 1e6, num["scan"]["w"] * 1e6, "C0o-"),
        ("analytical (native deck)", ana["x"] * 1e6, ana["scan"]["w"] * 1e6, "C1s--"),
    ],
    P.X_SPOT * 1e6, P.W0 * 1e6,
    figures / "waist_vs_x.png",
    title=f"Test 1: beam waist vs x (NA~{P.NA:.3f})",
)

lib.save_power_conservation(
    num["x"] * 1e6,
    [
        ("numerical", num["x"] * 1e6, num["ratio"], "C0-"),
        ("analytical", ana["x"] * 1e6, ana["ratio"], "C1--"),
    ],
    figures / "energy_conservation_ratio.png",
    title="Test 1: envelope-corrected transverse power ratio (should be ~1)",
)

lib.save_field_slice(num["ey_slice"], num["x"] * 1e6, num["y"] * 1e6,
                     "numerical: Ey(x,y,z=0)", figures / "ey_slice_numerical.png")
lib.save_field_slice(ana["ey_slice"], ana["x"] * 1e6, ana["y"] * 1e6,
                     "analytical: Ey(x,y,z=0)", figures / "ey_slice_analytical.png")

# --- Metrics -----------------------------------------------------------
valid_mask = lambda arr: np.isfinite(arr)


def rel_err(a, b):
    return abs(a - b) / abs(b) if (b and np.isfinite(a) and np.isfinite(b)) else np.nan


w_theory_at_x = lambda xarr: P.w_theory(xarr - P.X_SPOT)

num_w_theory = w_theory_at_x(num["x"])
ana_w_theory = w_theory_at_x(ana["x"])
num_valid = valid_mask(num["scan"]["w"])
ana_valid = valid_mask(ana["scan"]["w"])

num_w_rms_rel = float(np.sqrt(np.nanmean(
    ((num["scan"]["w"][num_valid] - num_w_theory[num_valid])
     / num_w_theory[num_valid]) ** 2)))
ana_w_rms_rel = float(np.sqrt(np.nanmean(
    ((ana["scan"]["w"][ana_valid] - ana_w_theory[ana_valid])
     / ana_w_theory[ana_valid]) ** 2)))
num_vs_ana_w_rms_rel = float(np.sqrt(np.nanmean(
    ((num["scan"]["w"][num_valid & ana_valid]
      - ana["scan"]["w"][num_valid & ana_valid])
     / ana["scan"]["w"][num_valid & ana_valid]) ** 2))) \
    if np.any(num_valid & ana_valid) else np.nan

# energy-conservation metric: interior only (exclude 5% of box at each end,
# where near-boundary transients / imperfect absorption dominate -- see
# lib.energy_conservation_ratio docstring)
interior = (num["x"] > 0.05 * P.BOX_X) & (num["x"] < 0.95 * P.BOX_X)
num_energy_dev = float(np.nanmax(np.abs(num["ratio"][interior] - 1.0)))
ana_energy_dev = float(np.nanmax(np.abs(ana["ratio"][interior] - 1.0)))

rows = [
    ("x_focus_theory_um", P.X_SPOT * 1e6),
    ("w0_theory_um", P.W0 * 1e6),
    ("x_focus_numerical_um", num["x_focus"] * 1e6),
    ("w0_numerical_um", num["w0"] * 1e6),
    ("x_focus_analytical_um", ana["x_focus"] * 1e6),
    ("w0_analytical_um", ana["w0"] * 1e6),
    ("x_focus_numerical_err_um", abs(num["x_focus"] - P.X_SPOT) * 1e6),
    ("w0_numerical_rel_err", rel_err(num["w0"], P.W0)),
    ("x_focus_analytical_err_um", abs(ana["x_focus"] - P.X_SPOT) * 1e6),
    ("w0_analytical_rel_err", rel_err(ana["w0"], P.W0)),
    ("w_scan_rms_rel_err_numerical_vs_theory", num_w_rms_rel),
    ("w_scan_rms_rel_err_analytical_vs_theory", ana_w_rms_rel),
    ("w_scan_rms_rel_err_numerical_vs_analytical", num_vs_ana_w_rms_rel),
    ("energy_conservation_max_dev_numerical_interior", num_energy_dev),
    ("energy_conservation_max_dev_analytical_interior", ana_energy_dev),
]
lib.write_metrics_csv(Path(results) / "metrics.csv", rows, ("metric", "value"))

summary = f"""Test 1 -- converging Gaussian beam, 3D field-level benchmark
=============================================================

Beam: lambda0={P.LAMBDA0*1e6:.3f} um, w0={P.W0*1e6:.3f} um, NA~{P.NA:.4f}
(moderate NA, chosen to stay within paraxial validity -- this validates
the FILE-INJECTION PIPELINE, not the paraxial approximation, and is
deliberately NOT the f/1 (NA=0.5) regime that motivated this feature).
Focus at x_spot={P.X_SPOT*1e6:.4f} um (box centre). Snapshot analysed:
t={num['t_snap']/1e-15:.2f} fs (target t_peak_at_focus={P.T_PEAK_AT_FOCUS/1e-15:.2f} fs).

Waist scan (2D Gaussian fit at each x-plane, see
common/viking_analysis_lib_3d.fit_transverse_gaussian):

  x_focus: theory={P.X_SPOT*1e6:.4f} um
           numerical={num['x_focus']*1e6:.4f} um (err {abs(num['x_focus']-P.X_SPOT)*1e6:.4f} um)
           analytical={ana['x_focus']*1e6:.4f} um (err {abs(ana['x_focus']-P.X_SPOT)*1e6:.4f} um)

  w0:      theory={P.W0*1e6:.4f} um
           numerical={num['w0']*1e6:.4f} um (rel err {rel_err(num['w0'],P.W0)*100:.3f}%)
           analytical={ana['w0']*1e6:.4f} um (rel err {rel_err(ana['w0'],P.W0)*100:.3f}%)

  w(x) RMS relative error vs theory across the full scan:
           numerical  = {num_w_rms_rel*100:.3f}%
           analytical = {ana_w_rms_rel*100:.3f}%
           numerical vs analytical (isolates injector-pipeline error) = {num_vs_ana_w_rms_rel*100:.3f}%

Energy conservation (transverse power P(x), corrected for the imposed
temporal-envelope shape -- see lib.energy_conservation_ratio; a raw,
uncorrected P(x) is NOT expected to be flat for a pulsed non-CW beam and
should not be used as the check):
  max |ratio-1| over the interior 90% of the box (excludes near-boundary
  transient and far-boundary absorption artefacts):
           numerical  = {num_energy_dev:.4f}
           analytical = {ana_energy_dev:.4f}

Interpretation:
  * If "numerical vs analytical" RMS error is small (comparable to or
    smaller than "numerical vs theory" and "analytical vs theory"), the
    file-injection pipeline (interpolation, binary format, static
    phase-from-file path) is behaving correctly -- any residual error vs
    theory is paraxial-approximation error common to both, not a pipeline
    bug.
  * If "numerical vs analytical" is markedly larger than either vs-theory
    number, that specifically implicates the file-injection pipeline.
  * Energy-conservation ratio should sit close to 1.0 (interior); a
    systematic trend (not just noise) indicates the injected wavefront is
    not self-consistent (e.g. amplitude/phase mismatch), independent of
    whether the paraxial formula itself is exactly right.

Figures: results/figures/waist_vs_x.png, energy_conservation_ratio.png,
ey_slice_numerical.png, ey_slice_analytical.png
Metrics: results/metrics.csv
"""
lib.write_text(Path(results) / "summary.txt", summary)
print(summary)
