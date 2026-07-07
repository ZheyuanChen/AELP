"""
Test A (3D) profile generator -- static (y,z) amplitude + phase planes.

Writes into amp_file_phase_file/ (needs both) and amp_file_phase_deck/
(needs amplitude only). amp_deck_phase_deck/ needs no files at all.
amp_deck_phase_file/ is unreachable with the current code -- see its
README.md.

Uses the EXACT SAME formula for the file content as the native deck
expressions in the *_deck decks (literal constants copied from
physics_params.py's printed values) -- the point of this test is to
isolate the injector's own behaviour (interpolation, binary I/O), so both
paths must start from identical numbers, not independently-computed ones
that could differ by floating-point evaluation order.

File convention: Fortran (n_tr1, n_tr2)=(y,z), tr1 fastest-varying,
written as numpy (n_tr2, n_tr1)=(z,y) via .tofile(), no transpose.
"""
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import physics_params as P

N_FILE = 1200
HALF_TRANSVERSE = 3.5 * P.W_BND

y = np.linspace(-HALF_TRANSVERSE, HALF_TRANSVERSE, N_FILE)
z = np.linspace(-HALF_TRANSVERSE, HALF_TRANSVERSE, N_FILE)
Y, Z = np.meshgrid(y, z)
rho2 = Y ** 2 + Z ** 2

amplitude = np.exp(-rho2 / P.W_BND ** 2)
phase = P.PHASE_CONST_3D + P.PHASE_QUAD * rho2

HERE = Path(__file__).parent
for d, need_phase in [("amp_file_phase_file", True), ("amp_file_phase_deck", False)]:
    out = HERE / d
    amplitude.astype(np.float64).tofile(out / "spatial_profile.dat")
    print(f"Wrote {out/'spatial_profile.dat'}")
    if need_phase:
        phase.astype(np.float64).tofile(out / "phase_profile.dat")
        print(f"Wrote {out/'phase_profile.dat'}")

print(f"\nGrid: {N_FILE}x{N_FILE}, extent +/-{HALF_TRANSVERSE*1e6:.4f} um "
      f"(== box transverse extent, no clamp/zero-pad edge case)")
print(f"amplitude peak={amplitude.max():.6f}, edge={amplitude[0,0]:.3e}")
print(f"phase range=[{phase.min():.4f}, {phase.max():.4f}] rad, "
      f"centre={P.PHASE_CONST_3D:.6f}")