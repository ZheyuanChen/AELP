"""
Test 2 -- 3D two-channel polarisation independence (rigorous / Viking-scale).

This is the 3D, high-resolution repeat of the exact recipe that caught a
REAL bug in epoch2d (module-level SAVE storage for custom_laser.f90's
profile/phase matrices -- a global singleton silently shared across every
laser block on a boundary). The epoch3d port was built per-laser from the
start (fields moved onto laser_block in shared_data.F90, not module-level),
so the bug should not reproduce -- but that is a code-reading inference,
not a measurement, and this exact class of bug is specifically silent
(wrong answer, no crash, no error) rather than loud. Costs minutes to
actually check.

Design (mirrors the validated 2D recipe exactly, see daily_log/2026/06/30.md
"Validated the epoch_dev per-laser-storage fix" / "...new raw-binary
profile/phase format"): two x_min laser blocks, pol_angle 0 and pi/2,
pointing at spatiotemporal profile files with the SAME grid size (so the
file-size-mismatch abort guard does NOT fire) but DIFFERENT content --
deliberately isolating the dangerous failure mode (silent wrong answer,
laser B's data replaced by a stale copy of laser A's) rather than the
already-covered abort case (mismatched grid sizes).

  Laser A: pol_angle=0      -> should appear in Ey. Narrow, off-centre,
           peak amp scale 1.0.
  Laser B: pol_angle=pi/2   -> should appear in Ez. Wide, off-centre
           (different location), peak amp scale 0.4.

Both files are time-independent (same 2D shape at every t-slice) -- this
test isolates per-laser storage/data independence, not time-axis handling
(that is covered by Test 1's y_min boundary array-ordering exercise
elsewhere in this project's local dev_test battery).

File convention: spatiotemporal, Fortran (n_tr1, n_tr2, n_t), written as
numpy (n_t, n_tr2, n_tr1) via .tofile(), no transpose (see custom_laser.f90).
"""
import numpy as np
from pathlib import Path

micron = 1e-6

N_TR1, N_TR2, N_T = 320, 320, 3          # SAME size for both lasers
TR_MIN, TR_MAX = -8.0 * micron, 8.0 * micron   # matches the sim box exactly

y = np.linspace(TR_MIN, TR_MAX, N_TR1)
z = np.linspace(TR_MIN, TR_MAX, N_TR2)
Y, Z = np.meshgrid(y, z)  # shape (n_tr2, n_tr1)

# Laser A: narrow, peak at (y,z) = (+2, +1) um, w=0.5 um, amp scale 1.0
yA, zA, wA, ampA = 2.0 * micron, 1.0 * micron, 0.5 * micron, 1.0
profile_A = ampA * np.exp(-(((Y - yA) ** 2 + (Z - zA) ** 2)) / wA ** 2)

# Laser B: wide, peak at (y,z) = (-2, -1.5) um, w=1.5 um, amp scale 0.4
yB, zB, wB, ampB = -2.0 * micron, -1.5 * micron, 1.5 * micron, 0.4
profile_B = ampB * np.exp(-(((Y - yB) ** 2 + (Z - zB) ** 2)) / wB ** 2)

# Time-independent: replicate the same 2D plane at every t-slice
array_A = np.repeat(profile_A[np.newaxis, :, :], N_T, axis=0)
array_B = np.repeat(profile_B[np.newaxis, :, :], N_T, axis=0)

HERE = Path(__file__).parent
out = HERE / "run"
out.mkdir(exist_ok=True)
array_A.astype(np.float64).tofile(out / "profile_A.dat")
array_B.astype(np.float64).tofile(out / "profile_B.dat")

expected_fwhm_A = 1.665 * wA
expected_fwhm_B = 1.665 * wB
expected_ratio = ampB / ampA

print(f"Wrote {out/'profile_A.dat'}  grid {N_TR1}x{N_TR2}x{N_T}  "
      f"({array_A.nbytes/1e6:.1f} MB)")
print(f"Wrote {out/'profile_B.dat'}  grid {N_TR1}x{N_TR2}x{N_T}  "
      f"({array_B.nbytes/1e6:.1f} MB)")
print()
print(f"Laser A (pol_angle=0 -> Ey): peak=({yA/micron:.2f},{zA/micron:.2f})um  "
      f"w={wA/micron:.2f}um  FWHM={expected_fwhm_A/micron:.4f}um  amp_scale={ampA}")
print(f"Laser B (pol_angle=pi/2 -> Ez): peak=({yB/micron:.2f},{zB/micron:.2f})um  "
      f"w={wB/micron:.2f}um  FWHM={expected_fwhm_B/micron:.4f}um  amp_scale={ampB}")
print(f"Expected amplitude ratio (B/A) = {expected_ratio:.4f}")
print()
print("Failure signature to watch for (the exact epoch2d bug this repeats):")
print("  if the per-laser fix has regressed, Ez would show laser A's shape")
print("  (peak near (+2,+1)um, FWHM matching w=0.5um) instead of laser B's.")