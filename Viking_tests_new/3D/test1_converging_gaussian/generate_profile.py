"""
Test 1 profile generator -- static (y,z) amplitude + phase planes for the
converging-Gaussian boundary condition, encoding the exact same formula as
analytical/input.deck's native `gauss(r_yz,...)` / phase expression (see
physics_params.py for the derivation). Writes into numerical/.

File convention (epoch3d, v2.1.0, custom_laser.f90 -- static spatial path):
  Fortran array (n_tr1, n_tr2), tr1 (=y, for an x_min laser) fastest-varying.
  Written as a numpy array of shape (n_tr2, n_tr1) = (z, y) via .tofile()
  in C order -- NO transpose needed (numpy's slow axis lines up with
  Fortran's slow axis). This is the same no-transpose rule verified for the
  spatiotemporal (n_t, n_tr2, n_tr1) case, one axis fewer.
"""
import numpy as np
from pathlib import Path
import physics_params as P

y = np.linspace(P.FILE_TR_MIN, P.FILE_TR_MAX, P.N_FILE)
z = np.linspace(P.FILE_TR_MIN, P.FILE_TR_MAX, P.N_FILE)
Y, Z = np.meshgrid(y, z)  # shape (n_tr2, n_tr1) = (z, y)
rho2 = Y ** 2 + Z ** 2

amplitude = np.exp(-rho2 / P.W_BND ** 2)
phase = P.EPOCH_PHASE_CONST + P.EPOCH_PHASE_QUAD * rho2

out_dir = Path(__file__).parent / "numerical"
out_dir.mkdir(exist_ok=True)

amp_path = out_dir / "spatial_profile.dat"
phase_path = out_dir / "phase_profile.dat"
amplitude.astype(np.float64).tofile(amp_path)
phase.astype(np.float64).tofile(phase_path)

print(f"Wrote {amp_path}  ({amplitude.nbytes/1e6:.1f} MB)")
print(f"Wrote {phase_path}  ({phase.nbytes/1e6:.1f} MB)")
print(f"Grid: {P.N_FILE} x {P.N_FILE}, extent [{P.FILE_TR_MIN*1e6:.4f}, "
      f"{P.FILE_TR_MAX*1e6:.4f}] um")
print(f"amplitude: peak={amplitude.max():.6f} at centre, "
      f"edge={amplitude[0,0]:.3e} (should be ~0)")
print(f"phase: range [{phase.min():.4f}, {phase.max():.4f}] rad "
      f"(centre={P.EPOCH_PHASE_CONST:.4f})")
print()
print("--- Matching deck fragment (numerical/input.deck laser block) ---")
print(f"""
  n_y_points = {P.N_FILE}
  n_z_points = {P.N_FILE}
  y_min = {P.FILE_TR_MIN:.8e}
  y_max = {P.FILE_TR_MAX:.8e}
  z_min = {P.FILE_TR_MIN:.8e}
  z_max = {P.FILE_TR_MAX:.8e}
""")