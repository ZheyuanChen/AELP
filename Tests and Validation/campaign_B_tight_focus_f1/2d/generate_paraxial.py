"""
Test B (2D) -- paraxial closed-form f/1 amplitude+phase generator.

Same formula as EPOCH's native deck laser block would use (profile =
gauss(y,0,w_bnd), phase = k0*y^2/(2*R_bnd) - gouy_2d), on the SAME (y,t)
grid as generate_lasy.py -- run that script first, this one reads the grid
straight out of lasy/spatial_profile.dat's shape (via the deck fragment
printed there) so the two are guaranteed to share identical sample points.

Deliberate correction vs the 30 June session's generate_paraxial_profile.py:
uses HALF Gouy phase (gouy_2d), the physically-correct convention for
EPOCH2D's 1D-transverse slab beam -- see physics_params.py docstring.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import physics_params as P

N_Y = 1600
Y_MIN, Y_MAX = -P.HALF_TRANSVERSE, P.HALF_TRANSVERSE
N_T = P.N_T

y_arr = np.linspace(Y_MIN, Y_MAX, N_Y)
t_arr = np.linspace(P.T_START, P.T_END, N_T)
Tg, Yg = np.meshgrid(t_arr, y_arr, indexing="ij")

amp = np.exp(-(Yg / P.W_BND) ** 2) * np.exp(-((Tg - P.T_CENTRE) / P.PULSE_TAU) ** 2)
# Sign convention matches injector_2x2_validation/physics_params.py (empirically
# verified there: the beam measurably converges to the right focus in EPOCH).
# NOTE: R_BND is negative here (boundary is before focus), so -k0/(2*R_BND) is
# positive -- this is NOT the same as the 30 June script's "+k0*y^2/(2*RC)" with
# RC defined as a positive foc_dist-based quantity; despite R_BND = -RC
# algebraically, transcribing the old formula naively (without flipping the
# sign) into this R_BND-based convention is an easy transcription error -- it
# was made once while building this script and caught by the local smoke test
# below (the beam diverged instead of converging). Kept as a comment as a
# specific instance of the sign-convention bug class this project keeps
# hitting (see CLAUDE.md's array-ordering/sign-convention emphasis).
phase = -P.K0 * Yg ** 2 / (2.0 * P.R_BND) - P.GOUY_2D  # no time dependence

HERE = Path(__file__).parent / "paraxial"
HERE.mkdir(exist_ok=True)
amp.astype(np.float64).tofile(HERE / "spatial_profile.dat")
phase.astype(np.float64).tofile(HERE / "phase_profile.dat")
print(f"Wrote {HERE/'spatial_profile.dat'} and phase_profile.dat "
      f"(grid {N_Y} x {N_T})")
print(f"w_bnd={P.W_BND*1e6:.4f}um  R_bnd={P.R_BND*1e6:.4f}um  "
      f"gouy_2d={P.GOUY_2D:.4f}rad")
