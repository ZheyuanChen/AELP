"""
Test A -- injector 2x2 validation (amplitude x phase, file vs native deck).
Shared physics for both the 2D and 3D campaigns (same lambda0, w0, x_spot,
"cast together" so results are directly comparable across dimensionality).

Purpose: isolate whether the FILE-INJECTION PIPELINE (interpolation,
binary format, static/phase-from-file paths) reproduces EPOCH's own native
deck-expression evaluation -- NOT whether paraxial theory is accurate.
Moderate NA is used deliberately (NA~0.265, well inside paraxial validity)
so the native deck expression itself is close to "ground truth": any
disagreement between the file-injected and native-deck runs is then
attributable to the injector, not to a paraxial-approximation gap. The
f/1 (NA=0.5) paraxial-vs-non-paraxial question is Test B's job, not this
test's.

2x2 design (both 2D and 3D): amplitude source in {file, deck} independently
crossed with phase source in {file, deck} -- 4 runs per dimensionality.
This is a deliberate improvement over the isolate-the-source-not-the-path
idea in the old test2_lasy_2x2 package (LASY vs analytical): that
compared two different SOURCE GENERATORS (LASY vs closed-form) which
conflates generator error with injector error. Here both sides use the
IDENTICAL closed-form formula -- the only thing that varies is whether it
reaches EPOCH via a binary file or a native deck expression -- so any
residual is unambiguously a pipeline effect.

Dimensionality note (Gouy phase): the Gouy phase for a Gaussian beam scales
as (D/2)*atan(z/zR), where D is the number of transverse dimensions
(Siegman; Self, Appl. Opt. 22, 658 (1983)). EPOCH3D models a genuinely 2D
transverse (D=2, circularly symmetric r_yz) beam -- full Gouy phase,
atan(z/zR). EPOCH2D models a 1D transverse slab (D=1, y only, uniform/
infinite in z) -- HALF Gouy phase, 0.5*atan(z/zR). The 30 June
`dev_test/laser_profile_injection/tight_focusing_f1/generate_paraxial_
profile.py` script used the full (D=2) Gouy phase in an EPOCH2D (D=1)
context -- a likely inconsistency, NOT reproduced here. Since Gouy phase
is a spatially-uniform constant at the boundary, it does not affect
amplitude/convergence physics (only an overall carrier-phase reference),
so it did not invalidate that script's amplitude/wing-structure findings,
but the 2D formula here is deliberately corrected.

EPOCH phase-sign convention (see the lasy-epoch-field-convention memory
and physics_params.py in Viking_results/3D/test1_converging_gaussian/,
where this exact formula was empirically validated -- the beam measurably
converges to the right focus location and width in a real EPOCH3D run,
not just reasoned about): EPOCH's E_phys = amp*profile*sin(wt+phase);
physical field = Re[E_env*exp(-i*wt)]. The r-independent constant terms
(pi/2, Gouy phase) only shift the overall carrier-phase reference and do
NOT affect focusing -- confirmed empirically in that Test 1 run.
"""
import numpy as np

LAMBDA0 = 1.0e-6
W0 = 1.2e-6                    # same beam in 2D and 3D
K0 = 2.0 * np.pi / LAMBDA0
X_R = np.pi * W0 ** 2 / LAMBDA0
X_SPOT = 2.0 * X_R             # boundary-to-focus distance; focus at box centre
NA = LAMBDA0 / (np.pi * W0)

C_LIGHT = 299792458.0

# Temporal envelope (both dimensionalities): peak reaches focus near the
# run's midpoint, matching Viking_results/3D/test1_converging_gaussian.
TAU = 80e-15
T_CENTRE = 200e-15
T_PEAK_AT_FOCUS = T_CENTRE + X_SPOT / C_LIGHT
T_END = 350e-15
DT_SNAPSHOT = 25e-15

AMP = 3.2e12


def w_of_xi(xi):
    """Beam radius at beam-coordinate xi (xi=0 at focus), both dimensionalities."""
    return W0 * np.sqrt(1.0 + (xi / X_R) ** 2)


def R_of_xi(xi):
    return xi * (1.0 + (X_R / xi) ** 2)


def gouy_3d(xi):
    """Full (D=2) Gouy phase -- EPOCH3D, circularly-symmetric r_yz beam."""
    return np.arctan(xi / X_R)


def gouy_2d(xi):
    """Half (D=1) Gouy phase -- EPOCH2D, 1D-transverse slab beam."""
    return 0.5 * np.arctan(xi / X_R)


# --- Boundary-plane values (xi_bnd = -X_SPOT, before focus) ---
XI_BND = -X_SPOT
W_BND = w_of_xi(XI_BND)
R_BND = R_of_xi(XI_BND)
PSI_BND_3D = gouy_3d(XI_BND)
PSI_BND_2D = gouy_2d(XI_BND)

PHASE_CONST_3D = np.pi / 2.0 + PSI_BND_3D
PHASE_CONST_2D = np.pi / 2.0 + PSI_BND_2D
PHASE_QUAD = -K0 / (2.0 * R_BND)   # coefficient of rho^2 (3D) or y^2 (2D) -- same formula


def summary():
    return "\n".join([
        f"lambda0={LAMBDA0*1e6:.4f}um  w0={W0*1e6:.4f}um  NA={NA:.4f}",
        f"x_R={X_R*1e6:.4f}um  x_spot={X_SPOT*1e6:.4f}um  w_bnd={W_BND*1e6:.4f}um",
        f"R_bnd={R_BND*1e6:.4f}um",
        f"psi_bnd_3D(full)={PSI_BND_3D:.4f}rad  psi_bnd_2D(half)={PSI_BND_2D:.4f}rad",
        f"phase_const_3D={PHASE_CONST_3D:.4f}  phase_const_2D={PHASE_CONST_2D:.4f}",
        f"phase_quad(k0/2R, same both dims)={PHASE_QUAD:.6e} 1/m^2",
        f"amp={AMP:.4e} V/m  t_centre={T_CENTRE*1e15:.1f}fs  tau={TAU*1e15:.1f}fs",
        f"t_peak_at_focus={T_PEAK_AT_FOCUS*1e15:.2f}fs  t_end={T_END*1e15:.1f}fs",
    ])


if __name__ == "__main__":
    print(summary())