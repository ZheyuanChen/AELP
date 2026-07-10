"""
Campaign C (2D) -- LOCAL cross-check, no EPOCH run needed: the NEW in-lasy
exporter vs the LEGACY hand-rolled conversion pipeline (campaign B's
generate_lasy.py recipe: rt-mode lasy propagate, RegularGridInterpolator
r -> |y|, skimage 2D phase unwrap, manual piston pin and .tofile()).

Both pipelines model the SAME physical beam (campaign A's moderate-NA
Gaussian) and are pinned to the SAME carrier reference (pi/2 +
psi_bnd_2d), so their binary outputs should agree closely; the residual
is rt-vs-xyt propagation/grid differences, not conventions. The legacy
recipe was validated in real EPOCH runs (campaign B), so agreement here
transfers that validation to the new exporter before any Viking hours are
spent. Run AFTER generate_lasy_exporter.py (reads its .dat output):

    python generate_lasy_exporter.py
    python crosscheck_legacy_pipeline.py

Pass metrics (with the measured 10 July 2026 values for the record):
  - amplitude max |diff| < 2e-2 of peak        (measured: 5.8e-4)
  - amplitude-weighted field error, i.e. 2*amp*|sin(dphi/2)| -- the actual
    contribution a phase difference makes to the injected field --
    max < 2e-3 of peak                          (measured: 5.9e-4)
  - raw phase max |diff| (mod 2pi) < 5e-2 rad where amp > 5e-2
                                                (measured: 3.8e-3)
The raw phase difference is deliberately NOT asserted at the amp > 1e-3
level: out at ~3 w_bnd the envelope phase disagrees by up to ~0.15 rad
between the rt and xyt propagations (different grids/windows resolve the
wings differently), but the amplitude there is ~1e-3 so the field-level
effect is negligible -- hence the amplitude-weighted metric above.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from skimage.restoration import unwrap_phase
from lasy.laser import Laser
from lasy.profiles.gaussian_profile import GaussianProfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import physics_params as P

N_Y, N_T = P.N_TRANS_FILE_2D, P.N_T_FILE
SHAPE = (N_T, N_Y)  # C-order (t, transverse), both pipelines

# --- 1. Load the NEW exporter's output (written by generate_lasy_exporter)
lasy_dir = HERE / "lasy_file"
amp_new = np.fromfile(lasy_dir / "laser_amplitude.dat",
                      dtype=np.float64).reshape(SHAPE)
phase_new = np.fromfile(lasy_dir / "laser_phase.dat",
                        dtype=np.float64).reshape(SHAPE)

# --- 2. Rebuild the same beam via the LEGACY pipeline (campaign B recipe,
# campaign A physics, same target grid as the new exporter's file).
profile = GaussianProfile(P.LAMBDA0, (1, 0), P.LASER_ENERGY, P.W0,
                          P.TAU, P.T_CENTRE)
laser = Laser(dim="rt", lo=(0.0, P.T_START_FILE),
              hi=(1.2 * P.Y_HALF, P.T_END_FILE),
              npoints=(512, N_T), profile=profile)
print(f"Legacy pipeline: rt back-propagate {P.X_SPOT*1e6:.3f} um...")
laser.propagate(-P.X_SPOT)

env = laser.grid.get_temporal_field()[0, :, :]
r_axis = np.real(laser.grid.axes[0])
t_axis = np.real(laser.grid.axes[1])

amp_interp = RegularGridInterpolator((r_axis, t_axis), np.abs(env),
                                     bounds_error=False, fill_value=0.0)
phase_unwrapped = unwrap_phase(np.angle(env))
phase_interp = RegularGridInterpolator((r_axis, t_axis), phase_unwrapped,
                                       bounds_error=False, fill_value=0.0)

y_arr = np.linspace(-P.Y_HALF, P.Y_HALF, N_Y)
t_arr = np.real(t_axis)
Tg, Yg = np.meshgrid(t_arr, y_arr, indexing="ij")
Rg = np.abs(Yg)

amp_legacy = amp_interp((Rg, Tg))
amp_legacy /= amp_legacy.max()

phase_lasy = phase_interp((Rg, Tg))
it_ref = int(np.argmin(np.abs(t_arr - P.T_CENTRE)))
iy_ref = int(np.argmin(np.abs(y_arr - 0.0)))
phi_ref = phase_lasy[it_ref, iy_ref]
# Same pin as the new exporter's carrier_phase_ref = PSI_BND_2D (campaign
# B's cell pinned to -gouy WITHOUT the +pi/2 because its paraxial partner
# deck was built that way; campaign C follows campaign A's native-deck
# convention instead, so both pipelines here carry the +pi/2 explicitly).
phase_legacy = -(phase_lasy - phi_ref) + np.pi / 2 + P.PSI_BND_2D

# --- 3. Compare where the amplitude is significant -----------------------
significant = np.logical_and(amp_new > 1e-3, amp_legacy > 1e-3)
core = np.logical_and(amp_new > 5e-2, amp_legacy > 5e-2)
amp_diff = np.abs(amp_new - amp_legacy)
dphi = phase_new - phase_legacy
dphi = np.abs(dphi - 2 * np.pi * np.round(dphi / (2 * np.pi)))
# Contribution the phase difference actually makes to the injected field
# amp*sin(w0 t + phase): |e^{i dphi} - 1| * amp = 2 amp |sin(dphi/2)|.
field_err = 2.0 * amp_new * np.abs(np.sin(dphi / 2.0))

amp_max = float(amp_diff[significant].max())
amp_rms = float(np.sqrt(np.mean(amp_diff[significant] ** 2)))
phi_max_sig = float(dphi[significant].max())
phi_max_core = float(dphi[core].max())
ferr_max = float(field_err[significant].max())
ferr_rms = float(np.sqrt(np.mean(field_err[significant] ** 2)))

print(f"\nNew exporter vs legacy pipeline "
      f"({int(significant.sum())} significant samples, "
      f"{int(core.sum())} core samples):")
print(f"  amplitude:  max |diff| = {amp_max:.3e} of peak, rms = {amp_rms:.3e}")
print(f"  field error (amp-weighted phase): max = {ferr_max:.3e} of peak, "
      f"rms = {ferr_rms:.3e}")
print(f"  raw phase (mod 2pi): max = {phi_max_core:.3e} rad in the core "
      f"(amp > 5e-2), {phi_max_sig:.3e} rad out to amp > 1e-3 (far wings, "
      "see docstring)")

assert amp_max < 2e-2, f"amplitude mismatch {amp_max:.3e} >= 2e-2"
assert ferr_max < 2e-3, f"field-level mismatch {ferr_max:.3e} >= 2e-3"
assert phi_max_core < 5e-2, f"core phase mismatch {phi_max_core:.3e} >= 5e-2"
print("\nPASS: the new exporter reproduces the validated legacy pipeline.")
