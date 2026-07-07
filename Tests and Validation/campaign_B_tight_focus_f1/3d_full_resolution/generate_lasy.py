"""
Test B (3D) -- LASY (non-paraxial) f/1 amplitude+phase generator.

Same LASY rt-mode propagation as ../2d/generate_lasy.py (cylindrically
symmetric -> only needs a 1D radial profile), but resampled onto a genuine
Cartesian (y,z) grid via rho=sqrt(y^2+z^2) instead of a 1D y axis -- the
direct 3D extension, matching the resampling trick already used (and
validated) in the old Viking_results/2D/test2_lasy_2x2/generate_lasy_
profiles.py for the 2D case (r=|y| there; r=sqrt(y^2+z^2) here).

File convention (epoch3d spatiotemporal path): Fortran (n_tr1, n_tr2, n_t)
= (y,z,t), y fastest. Written as numpy (n_t, n_tr2, n_tr1) via .tofile(),
no transpose.
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

N_TR = 400
N_T = P.N_T
HALF_TRANSVERSE = P.HALF_TRANSVERSE

profile = GaussianProfile(P.LAMBDA0, (1, 0), 1.0, P.W0, P.PULSE_TAU, P.T_CENTRE)
laser = Laser(dim="rt", lo=(0, P.T_START), hi=(5 * P.W_BND, P.T_END),
              npoints=(N_TR // 2, N_T), profile=profile)
print("Propagating from focus to the injection boundary (3D)...")
laser.propagate(-P.X_SPOT)

env = laser.grid.get_temporal_field()[0, :, :]
r_axis = np.real(laser.grid.axes[0])
t_axis = np.real(laser.grid.axes[1])

amp_interp = RegularGridInterpolator((r_axis, t_axis), np.abs(env),
                                     bounds_error=False, fill_value=0.0)
phase_unwrapped = unwrap_phase(np.angle(env))
phase_interp = RegularGridInterpolator((r_axis, t_axis), phase_unwrapped,
                                       bounds_error=False, fill_value=0.0)

y = np.linspace(-HALF_TRANSVERSE, HALF_TRANSVERSE, N_TR)
z = np.linspace(-HALF_TRANSVERSE, HALF_TRANSVERSE, N_TR)
t_arr = np.real(t_axis)

Y2, Z2 = np.meshgrid(y, z)          # (n_tr2, n_tr1) = (z, y)
Rg2 = np.sqrt(Y2 ** 2 + Z2 ** 2)     # (n_tr2, n_tr1)

# Reference (piston removal) at r=0, t=t_centre, matching ../2d/generate_lasy.py
it_ref = int(np.argmin(np.abs(t_arr - P.T_CENTRE)))
phi_ref = float(phase_interp((0.0, t_arr[it_ref])))

amp_out = np.empty((N_T, N_TR, N_TR))
phase_out = np.empty((N_T, N_TR, N_TR))
for it, t in enumerate(t_arr):
    Tg2 = np.full_like(Rg2, t)
    amp_out[it] = amp_interp((Rg2, Tg2))
    phase_out[it] = -(phase_interp((Rg2, Tg2)) - phi_ref) + (-P.GOUY_3D)

amp_out /= amp_out.max()

print(f"amp peak (post-norm) = {amp_out.max():.4f}")
print(f"phi_ref (on-axis, t_centre) = {phi_ref:.4f} rad")
print(f"phase range = [{phase_out.min():.3f}, {phase_out.max():.3f}] rad")

HERE = Path(__file__).parent / "lasy"
HERE.mkdir(exist_ok=True)
amp_out.astype(np.float64).tofile(HERE / "spatial_profile.dat")
phase_out.astype(np.float64).tofile(HERE / "phase_profile.dat")
print(f"Wrote {HERE/'spatial_profile.dat'} and phase_profile.dat "
      f"(grid {N_TR}x{N_TR}x{N_T}, {amp_out.nbytes/1e6:.1f} MB each)")
