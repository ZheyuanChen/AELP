"""
Test B (3D) -- paraxial closed-form f/1 amplitude+phase generator.

3D (circularly-symmetric, rho=sqrt(y^2+z^2)) extension of ../2d/generate_
paraxial.py. Same corrected sign convention (empirically verified in both
injector_2x2_validation/ and ../2d/ -- the beam measurably converges to
the right focus; the naively-transcribed sign from the 30 June script's
RC-based formula was WRONG and caught by ../2d/'s local smoke test, see
that script's comment for the full story). Full (not half) Gouy phase,
matching EPOCH3D's genuinely 2D-transverse (circularly symmetric) beam.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import physics_params as P

N_TR = 400   # per transverse axis (y,z) -- 3D file is N_TR^2 * N_T, keep modest
N_T = P.N_T
HALF_TRANSVERSE = P.HALF_TRANSVERSE

y = np.linspace(-HALF_TRANSVERSE, HALF_TRANSVERSE, N_TR)
z = np.linspace(-HALF_TRANSVERSE, HALF_TRANSVERSE, N_TR)
t_arr = np.linspace(P.T_START, P.T_END, N_T)

Y, Z = np.meshgrid(y, z)   # (n_tr2, n_tr1) = (z, y)
rho2 = Y ** 2 + Z ** 2
spatial_amp = np.exp(-rho2 / P.W_BND ** 2)               # (n_tr2, n_tr1)
spatial_phase = -P.K0 * rho2 / (2.0 * P.R_BND) - P.GOUY_3D  # (n_tr2, n_tr1), no time dep

temporal_env = np.exp(-((t_arr - P.T_CENTRE) / P.PULSE_TAU) ** 2)  # (n_t,)

amp = spatial_amp[np.newaxis, :, :] * temporal_env[:, np.newaxis, np.newaxis]
phase = np.repeat(spatial_phase[np.newaxis, :, :], N_T, axis=0)

HERE = Path(__file__).parent / "paraxial"
HERE.mkdir(exist_ok=True)
amp.astype(np.float64).tofile(HERE / "spatial_profile.dat")
phase.astype(np.float64).tofile(HERE / "phase_profile.dat")
print(f"Wrote {HERE/'spatial_profile.dat'} and phase_profile.dat "
      f"(grid {N_TR}x{N_TR}x{N_T}, {amp.nbytes/1e6:.1f} MB each)")
print(f"w_bnd={P.W_BND*1e6:.4f}um  R_bnd={P.R_BND*1e6:.4f}um  gouy_3d={P.GOUY_3D:.4f}rad")
print(f"\nDeck fragment: n_y={N_TR} n_z={N_TR} n_t={N_T}  "
      f"y/z_min={-HALF_TRANSVERSE:.8e} y/z_max={HALF_TRANSVERSE:.8e}  "
      f"t_start={P.T_START:.8e} t_end={P.T_END:.8e}")
