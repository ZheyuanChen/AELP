"""
Test A (2D) profile generator -- static y amplitude + phase 1D arrays.

epoch2d's static path (converted to raw binary, commit `cb57b7c1`) is the
1D-transverse analogue of epoch3d's static (y,z) plane: a single array of
n_y values, no header, y fastest (only) axis.

Uses the 2D-specific (HALF Gouy phase, see physics_params.py) formula, with
the same lambda0/w0/x_spot as the 3D campaign so the two are directly
comparable ("cast together").
"""
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import physics_params as P

N_FILE = 1200
HALF_TRANSVERSE = 3.5 * P.W_BND

y = np.linspace(-HALF_TRANSVERSE, HALF_TRANSVERSE, N_FILE)
amplitude = np.exp(-(y ** 2) / P.W_BND ** 2)
phase = P.PHASE_CONST_2D + P.PHASE_QUAD * y ** 2

HERE = Path(__file__).parent
for d, need_phase in [("amp_file_phase_file", True), ("amp_file_phase_deck", False)]:
    out = HERE / d
    amplitude.astype(np.float64).tofile(out / "spatial_profile.dat")
    print(f"Wrote {out/'spatial_profile.dat'}")
    if need_phase:
        phase.astype(np.float64).tofile(out / "phase_profile.dat")
        print(f"Wrote {out/'phase_profile.dat'}")

print(f"\nGrid: {N_FILE}, extent +/-{HALF_TRANSVERSE*1e6:.4f} um")
print(f"amplitude peak={amplitude.max():.6f}, edge={amplitude[0]:.3e}")
print(f"phase range=[{phase.min():.4f}, {phase.max():.4f}] rad, "
      f"centre={P.PHASE_CONST_2D:.6f}")
