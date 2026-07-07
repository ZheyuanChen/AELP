"""
Test B -- f/1 tight-focusing demonstration: LASY (non-paraxial) vs the
paraxial closed-form formula, BOTH injected into EPOCH via the file
injector only (no native-deck comparison -- see injector_2x2_validation/
for that question; this test is about which PHYSICS model is closer to
reality at NA=0.5, not about injector correctness).

Builds directly on `dev_test/laser_profile_injection/tight_focusing_f1/`
(30 June session): same base parameters (lambda0, pulse_tau, t_centre,
f-number geometry rule), same LASY pipeline (rt-mode propagate, 2D
quality-guided phase unwrap via skimage, piston-removal pinned to -gouy
on-axis). That prior work only ever produced/compared raw numpy arrays
(amplitude_lasy.dat vs amplitude_paraxial.dat, never run through EPOCH).
This extends it to: (1) the new raw-binary file format (old text format no
longer read by current epoch2d/epoch3d), (2) actually injecting both into
EPOCH and running, (3) a genuine 3D (y,z) version (old work was 2D/y-only).

Deliberate correction vs the 30 June 2D paraxial script: this uses the
dimension-correct Gouy phase (HALF for 2D-slab, FULL for 3D
circularly-symmetric) -- see injector_2x2_validation/physics_params.py's
docstring for the derivation. Gouy phase is a spatially-uniform constant
at the boundary and does not affect focusing physics, so this does not
change the 30 June session's amplitude/wing-structure conclusions, only
the absolute phase reference.
"""
import numpy as np

LAMBDA0 = 1.0e-6
PULSE_TAU = 15.0e-15
T_CENTRE = 40.0e-15
F_NUMBER = 1.0
NA = 1.0 / (2.0 * F_NUMBER)
W0 = LAMBDA0 / (np.pi * NA)
X_R = np.pi * W0 ** 2 / LAMBDA0
X_SPOT = 4.0 * X_R              # boundary-to-focus distance (matches 30 June choice)
XI_BND = -X_SPOT
K0 = 2.0 * np.pi / LAMBDA0

W_BND = W0 * np.sqrt(1.0 + (XI_BND / X_R) ** 2)
R_BND = XI_BND * (1.0 + (X_R / XI_BND) ** 2)
GOUY_3D = np.arctan(XI_BND / X_R)       # full, EPOCH3D (circularly symmetric)
GOUY_2D = 0.5 * np.arctan(XI_BND / X_R)  # half, EPOCH2D (1D-transverse slab)

MARGIN = 4.0
HALF_TRANSVERSE = MARGIN * W_BND

C_LIGHT = 299792458.0
T_START, T_END, N_T = 0.0, 6.0 * PULSE_TAU, 500


def summary():
    return "\n".join([
        f"f/{F_NUMBER:.1f}  NA={NA:.4f}",
        f"w0={W0*1e6:.4f}um  x_R={X_R*1e6:.4f}um  x_spot={X_SPOT*1e6:.4f}um "
        f"({X_SPOT/X_R:.2f} x_R)",
        f"w_bnd={W_BND*1e6:.4f}um  R_bnd={R_BND*1e6:.4f}um",
        f"gouy_3d={GOUY_3D:.4f}rad  gouy_2d(half)={GOUY_2D:.4f}rad",
        f"half_transverse={HALF_TRANSVERSE*1e6:.4f}um (margin={MARGIN}x w_bnd)",
        f"pulse_tau={PULSE_TAU*1e15:.1f}fs  t_centre={T_CENTRE*1e15:.1f}fs",
    ])


if __name__ == "__main__":
    print(summary())
