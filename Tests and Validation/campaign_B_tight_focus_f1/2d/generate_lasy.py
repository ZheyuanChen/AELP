"""
Test B (2D) -- LASY (non-paraxial) f/1 amplitude+phase generator.

Physics and pipeline are unchanged from the 30 June session's
`dev_test/laser_profile_injection/tight_focusing_f1/generate_f1_lasy_
profile.py` (same LASY rt-mode propagate, same skimage 2D phase unwrap
fix, same piston-removal convention) -- only the OUTPUT FORMAT changed
(raw binary, matching current epoch2d, instead of the old text format it
can no longer read) and physics_params.py is now the single shared source
of constants (was duplicated across scripts before).

File convention (epoch2d spatiotemporal path): Fortran (n_y, n_t), y
fastest-varying. Written as numpy (n_t, n_y) via .tofile(), no transpose.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from skimage.restoration import unwrap_phase
from lasy.laser import Laser
from lasy.profiles.gaussian_profile import GaussianProfile

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import physics_params as P

N_Y = 1600
Y_MIN, Y_MAX = -P.HALF_TRANSVERSE, P.HALF_TRANSVERSE

profile = GaussianProfile(P.LAMBDA0, (1, 0), 1.0, P.W0, P.PULSE_TAU, P.T_CENTRE)
laser = Laser(dim="rt", lo=(0, P.T_START), hi=(5 * P.W_BND, P.T_END),
              npoints=(N_Y // 2, P.N_T), profile=profile)
print("Propagating from focus to the injection boundary (2D)...")
laser.propagate(-P.X_SPOT)

env = laser.grid.get_temporal_field()[0, :, :]
r_axis = np.real(laser.grid.axes[0])
t_axis = np.real(laser.grid.axes[1])

amp_interp = RegularGridInterpolator((r_axis, t_axis), np.abs(env),
                                     bounds_error=False, fill_value=0.0)
phase_unwrapped = unwrap_phase(np.angle(env))
phase_interp = RegularGridInterpolator((r_axis, t_axis), phase_unwrapped,
                                       bounds_error=False, fill_value=0.0)

y_arr = np.linspace(Y_MIN, Y_MAX, N_Y)
t_arr = np.real(t_axis)
Tg, Yg = np.meshgrid(t_arr, y_arr, indexing="ij")
Rg = np.abs(Yg)

amp = amp_interp((Rg, Tg))
amp /= amp.max()

phase_lasy = phase_interp((Rg, Tg))
it_ref = int(np.argmin(np.abs(t_arr - P.T_CENTRE)))
iy_ref = int(np.argmin(np.abs(y_arr - 0.0)))
phi_ref = phase_lasy[it_ref, iy_ref]
# Pinned to the 2D (HALF) Gouy convention -- see physics_params.py docstring
# for why this deliberately differs from the 30 June script's full-Gouy pin.
phase_epoch = -(phase_lasy - phi_ref) + (-P.GOUY_2D)

print(f"amp peak (post-norm) = {amp.max():.4f}")
print(f"phi_ref (on-axis, t_centre) = {phi_ref:.4f} rad")
print(f"phase range = [{phase_epoch.min():.3f}, {phase_epoch.max():.3f}] rad")

HERE = Path(__file__).parent / "lasy"
HERE.mkdir(exist_ok=True)
amp.astype(np.float64).tofile(HERE / "spatial_profile.dat")
phase_epoch.astype(np.float64).tofile(HERE / "phase_profile.dat")
print(f"Wrote {HERE/'spatial_profile.dat'} and phase_profile.dat "
      f"(grid {N_Y} x {P.N_T})")
print(f"\nDeck fragment (n_y, n_t, y_min/max, t_start/t_end):")
print(f"  n_y = {N_Y}\n  n_t = {P.N_T}\n  y_min = {Y_MIN:.8e}\n  y_max = {Y_MAX:.8e}")
print(f"  t_start = {P.T_START:.8e}\n  t_end = {P.T_END:.8e}")
