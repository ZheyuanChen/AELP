"""
Test 2 analysis -- 3D two-channel polarisation independence.

Fits a 2D Gaussian to Ey (expected: laser A's shape) and Ez (expected:
laser B's shape) a few cells inside the x_min boundary, and checks each
against its OWN laser's encoded input -- and explicitly against the OTHER
laser's input, to make a regression (per-laser storage silently broken
again) structurally loud rather than just "doesn't match theory".

Usage: python analyse.py [base_dir]   (default: this script's directory)
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
import viking_analysis_lib_3d as lib

BASE = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
RUN = BASE / "run"
results, figures = lib.make_results_dirs(BASE)

micron = 1e-6

# Ground truth (MUST match generate_profiles.py)
laser_A = dict(y0=2.0 * micron, z0=1.0 * micron, w=0.5 * micron, amp=1.0)
laser_B = dict(y0=-2.0 * micron, z0=-1.5 * micron, w=1.5 * micron, amp=0.4)

ds = lib.load_fields(RUN, ["Electric_Field_Ey", "Electric_Field_Ez"])
ey = ds["Electric_Field_Ey"]
ez = ds["Electric_Field_Ez"]
x = ey["X_Grid_mid"].values
y = ey["Y_Grid_mid"].values
z = ey["Z_Grid_mid"].values
times = ey["time0"].values

# Sample a few cells inside the boundary (not at the boundary cell
# itself), at a snapshot once the field has locally established --
# per the project's own lesson (daily_log 2026-06-30): always sample
# right at/near the injection boundary, not an x-aggregate further in.
IX = 3
t_target = 15e-15
it = int(np.argmin(np.abs(times - t_target)))
print(f"Sampling at x={x[IX]*1e6:.4f} um (index {IX}), "
      f"t={times[it]/1e-15:.2f} fs (target {t_target/1e-15:.1f} fs)")

ey_slab = ey.isel(time0=it, X_Grid_mid=IX).values
ez_slab = ez.isel(time0=it, X_Grid_mid=IX).values

fit_ey = lib.fit_transverse_gaussian(np.abs(ey_slab), y, z)
fit_ez = lib.fit_transverse_gaussian(np.abs(ez_slab), y, z)

lib.save_field_slice(ey_slab, y * 1e6, z * 1e6,
                     f"Ey(y,z) at x={x[IX]*1e6:.2f}um, t={times[it]/1e-15:.1f}fs",
                     figures / "ey_slice.png")
lib.save_field_slice(ez_slab, y * 1e6, z * 1e6,
                     f"Ez(y,z) at x={x[IX]*1e6:.2f}um, t={times[it]/1e-15:.1f}fs",
                     figures / "ez_slice.png")


def fwhm_of(w):
    return 1.665 * w


def compare(fit, truth, label):
    fwhm_fit = fwhm_of(fit["w"])
    fwhm_truth = fwhm_of(truth["w"])
    dist_um = np.hypot(fit["y0"] - truth["y0"], fit["z0"] - truth["z0"]) / micron
    return dict(
        label=label,
        y0_fit_um=fit["y0"] / micron, z0_fit_um=fit["z0"] / micron,
        y0_truth_um=truth["y0"] / micron, z0_truth_um=truth["z0"] / micron,
        position_err_um=dist_um,
        fwhm_fit_um=fwhm_fit / micron, fwhm_truth_um=fwhm_truth / micron,
        fwhm_rel_err=abs(fwhm_fit - fwhm_truth) / fwhm_truth,
        amp_fit=fit["amp"],
    )


ey_vs_A = compare(fit_ey, laser_A, "Ey vs laser A (own input, expected match)")
ey_vs_B = compare(fit_ey, laser_B, "Ey vs laser B (cross-check, expected mismatch)")
ez_vs_B = compare(fit_ez, laser_B, "Ez vs laser B (own input, expected match)")
ez_vs_A = compare(fit_ez, laser_A, "Ez vs laser A (cross-check -- the exact epoch2d bug signature)")

amp_ratio_fit = fit_ez["amp"] / fit_ey["amp"] if fit_ey["amp"] else np.nan
amp_ratio_truth = laser_B["amp"] / laser_A["amp"]

rows = []
for c in (ey_vs_A, ey_vs_B, ez_vs_B, ez_vs_A):
    rows.append((c["label"], c["position_err_um"], c["fwhm_fit_um"],
                c["fwhm_truth_um"], c["fwhm_rel_err"]))
rows.append(("amplitude_ratio_Ez_over_Ey_fit", amp_ratio_fit, "", "", ""))
rows.append(("amplitude_ratio_truth_B_over_A", amp_ratio_truth, "", "", ""))
lib.write_metrics_csv(Path(results) / "metrics.csv", rows,
                      ("comparison", "position_err_um", "fwhm_fit_um",
                       "fwhm_truth_um", "fwhm_rel_err"))

pass_ey = ey_vs_A["position_err_um"] < 0.3 and ey_vs_A["fwhm_rel_err"] < 0.1
pass_ez = ez_vs_B["position_err_um"] < 0.3 and ez_vs_B["fwhm_rel_err"] < 0.1
regression_signature = (ez_vs_A["position_err_um"] < 0.3
                        and ez_vs_A["fwhm_rel_err"] < 0.1)

summary = f"""Test 2 -- 3D two-channel polarisation independence
====================================================

Sampled at x={x[IX]*1e6:.4f} um (index {IX} from x_min), t={times[it]/1e-15:.2f} fs.

Laser A (pol_angle=0, expected channel: Ey):
  input:  peak=({laser_A['y0']/micron:.2f},{laser_A['z0']/micron:.2f})um  w={laser_A['w']/micron:.2f}um  amp_scale={laser_A['amp']}
  Ey fit: peak=({fit_ey['y0']/micron:.4f},{fit_ey['z0']/micron:.4f})um  w={fit_ey['w']/micron:.4f}um  amp={fit_ey['amp']:.4e}
  Ey vs own input (laser A):   position err={ey_vs_A['position_err_um']:.4f} um, FWHM rel err={ey_vs_A['fwhm_rel_err']*100:.2f}%
  Ey vs OTHER input (laser B): position err={ey_vs_B['position_err_um']:.4f} um, FWHM rel err={ey_vs_B['fwhm_rel_err']*100:.2f}%
  -> {"PASS" if pass_ey else "FAIL"}: Ey matches laser A, not laser B.

Laser B (pol_angle=pi/2, expected channel: Ez):
  input:  peak=({laser_B['y0']/micron:.2f},{laser_B['z0']/micron:.2f})um  w={laser_B['w']/micron:.2f}um  amp_scale={laser_B['amp']}
  Ez fit: peak=({fit_ez['y0']/micron:.4f},{fit_ez['z0']/micron:.4f})um  w={fit_ez['w']/micron:.4f}um  amp={fit_ez['amp']:.4e}
  Ez vs own input (laser B):   position err={ez_vs_B['position_err_um']:.4f} um, FWHM rel err={ez_vs_B['fwhm_rel_err']*100:.2f}%
  Ez vs OTHER input (laser A -- the exact epoch2d bug signature):
                                position err={ez_vs_A['position_err_um']:.4f} um, FWHM rel err={ez_vs_A['fwhm_rel_err']*100:.2f}%
  -> {"PASS" if pass_ez else "FAIL"}: Ez matches laser B, not laser A.
  -> Regression signature present (Ez looks like laser A): {"YES -- BUG REPRODUCED" if regression_signature else "no"}

Amplitude ratio (Ez/Ey): fit={amp_ratio_fit:.4f}, truth (B/A)={amp_ratio_truth:.4f}
  rel err = {abs(amp_ratio_fit-amp_ratio_truth)/amp_ratio_truth*100:.2f}%

Overall: {"PASS -- per-laser storage independence confirmed at 3D/Viking scale" if (pass_ey and pass_ez and not regression_signature) else "FAIL -- investigate immediately, see regression signature above"}

Figures: results/figures/ey_slice.png, ez_slice.png
Metrics: results/metrics.csv
"""
lib.write_text(Path(results) / "summary.txt", summary)
print(summary)
