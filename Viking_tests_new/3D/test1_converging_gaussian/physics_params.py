"""
Test 1 — converging Gaussian beam: single source of truth for all physics
and grid constants.

Both `generate_profile.py` (writes the binary profile/phase files for the
numerical/ run) and the two decks (numerical/input.deck, analytical/
input.deck, hand-authored from this module's printed output) MUST agree on
every number here. Run this file directly (`python physics_params.py`) to
print everything needed to check or regenerate the decks.

Physics
-------
Standard paraxial Gaussian beam, propagating in +x, focus at x = x_focus
inside the box. Injection boundary is x_min = 0, i.e. at beam-coordinate
xi_bnd = -x_spot (before the focus -- the beam is converging as it enters
the box).

  w(xi)   = w0 * sqrt(1 + (xi/x_R)^2)             beam radius
  R(xi)   = xi * (1 + (x_R/xi)^2)                 wavefront radius of curvature
  psi(xi) = arctan(xi/x_R)                        Gouy phase

Field (paraxial, SVEA): E(rho,xi) ~ (w0/w(xi)) * exp(-rho^2/w(xi)^2)
    * exp(i*[k*xi - psi(xi) + k*rho^2/(2R(xi))])          [+ exp(-i*omega*t)]
with rho^2 = y^2 + z^2 (cylindrically symmetric -- moderate NA, scalar OK).

EPOCH-file convention (see TIGHT_FOCUSING_INJECTION_PROMPT.md / the
lasy-epoch-field-convention memory): EPOCH's E_phys = amp*profile*sin(wt+phase),
physical field = Re[E_env * exp(-i*omega*t)]. Matching:
  |E| cos(wt-phi) = |E| sin(wt + (pi/2-phi))  =>  phase = pi/2 - phi
The spatial phase at the boundary (excluding the internally-generated wt
carrier) is phi(rho) = k*rho^2/(2*R(xi_bnd)) - psi(xi_bnd), so:
  EPOCH phase(y,z) = pi/2 + psi(xi_bnd) - k*(y^2+z^2)/(2*R(xi_bnd))
This is a CONSTANT-plus-quadratic form -- directly expressible as a native
EPOCH deck expression (analytical/input.deck) using the built-in r_yz
variable, or as a binary (y,z) plane file (numerical/input.deck).

Validated quantities (see analyse.py):
  1. w(x) fit at several x-planes vs the w(xi) formula above.
  2. x_focus (location of minimum w) vs x_spot.
  3. w0_fit (minimum waist) vs w0.
  4. Transverse-plane power integral P(x) = integral |envelope|^2 dy dz,
     which should be flat (energy-conserving) vs x -- independent of
     whether the paraxial approximation itself is exact, this checks
     internal self-consistency of the injected wavefront.
"""
import numpy as np

# ---------------------------------------------------------------------
# Beam parameters (moderate NA -- paraxial approximation should hold to
# a few % here; deliberately NOT the f/1 (NA=0.5) regime, so this test
# validates the *pipeline*, not the paraxial approximation itself).
# ---------------------------------------------------------------------
LAMBDA0 = 1.0e-6          # wavelength [m]
W0 = 1.2e-6               # focal waist (1/e field radius) [m]
K0 = 2.0 * np.pi / LAMBDA0
X_R = np.pi * W0**2 / LAMBDA0             # Rayleigh range
X_SPOT = 2.0 * X_R                        # boundary-to-focus distance
XI_BND = -X_SPOT                          # boundary position rel. to focus
W_BND = W0 * np.sqrt(1.0 + (XI_BND / X_R) ** 2)
R_BND = XI_BND * (1.0 + (X_R / XI_BND) ** 2)
PSI_BND = np.arctan(XI_BND / X_R)
NA = LAMBDA0 / (np.pi * W0)               # small-angle NA estimate

EPOCH_PHASE_CONST = np.pi / 2.0 + PSI_BND    # constant term of phase(y,z)
EPOCH_PHASE_QUAD = -K0 / (2.0 * R_BND)       # coefficient of (y^2+z^2)

# ---------------------------------------------------------------------
# Simulation box: x in [0, 2*X_SPOT] (focus exactly at box centre, so the
# beam is seen converging over the first half and diverging over the
# second -- symmetric, self-checking). Transverse: +/- 3.5*W_BND, well
# beyond where the Gaussian tail is negligible (exp(-3.5^2) ~ 7.5e-6) so
# the static-path clamp behaviour at the file/box edge (see the logged
# "out-of-file-grid edge behaviour" issue) never triggers.
# ---------------------------------------------------------------------
CELLS_PER_LAMBDA = 24
DX = LAMBDA0 / CELLS_PER_LAMBDA

BOX_X = 2.0 * X_SPOT
HALF_TRANSVERSE = 3.5 * W_BND
BOX_TRANSVERSE = 2.0 * HALF_TRANSVERSE

NX = int(round(BOX_X / DX))
NY = int(round(BOX_TRANSVERSE / DX))
NZ = NY
TOTAL_CELLS = NX * NY * NZ

# ---------------------------------------------------------------------
# Binary file grid (static (y,z) plane, matches the whole simulation
# transverse extent exactly -- no clamp/zero-pad edge case at all -- at
# a resolution well above the simulation grid so interpolation error is
# negligible compared to the physics being tested).
# ---------------------------------------------------------------------
N_FILE = 1200               # per transverse axis
FILE_TR_MIN = -HALF_TRANSVERSE
FILE_TR_MAX = HALF_TRANSVERSE

# ---------------------------------------------------------------------
# Temporal envelope (deck t_profile, both numerical and analytical decks):
# chosen so the pulse PEAK reaches the focus (x = X_SPOT) close to the
# middle of the run, giving the best SNR exactly where w(x) is smallest
# and hardest to resolve.
# ---------------------------------------------------------------------
C_LIGHT = 299792458.0
TAU = 80e-15
T_CENTRE = 200e-15
T_PEAK_AT_FOCUS = T_CENTRE + X_SPOT / C_LIGHT
T_END = 350e-15
DT_SNAPSHOT = 25e-15

# Deck 'amp' (V/m) -- arbitrary convenient scale (~a0 order 1 at 1 um).
AMP = 3.2e12


def summary():
    lines = [
        f"lambda0        = {LAMBDA0*1e6:.6f} um",
        f"w0             = {W0*1e6:.6f} um",
        f"NA (approx)    = {NA:.6f}",
        f"x_R            = {X_R*1e6:.6f} um",
        f"x_spot         = {X_SPOT*1e6:.6f} um   (focus location, box centre)",
        f"w_bnd          = {W_BND*1e6:.6f} um",
        f"R_bnd          = {R_BND*1e6:.6f} um",
        f"psi_bnd        = {PSI_BND:.6f} rad",
        f"EPOCH phase(y,z) = {EPOCH_PHASE_CONST:.6f} + ({EPOCH_PHASE_QUAD:.6e})*(y^2+z^2)   [rad, y,z in m]",
        f"k0             = {K0:.6e} 1/m",
        "",
        f"dx = dy = dz   = {DX*1e9:.4f} nm  ({CELLS_PER_LAMBDA} cells/lambda)",
        f"box_x (full)   = {BOX_X*1e6:.6f} um   -> nx = {NX}",
        f"box_transverse = {BOX_TRANSVERSE*1e6:.6f} um   -> ny = nz = {NY}",
        f"total cells    = {TOTAL_CELLS:,}  ({TOTAL_CELLS/1e6:.2f} M)",
        "",
        f"file grid      = {N_FILE} x {N_FILE}, extent [{FILE_TR_MIN*1e6:.4f}, {FILE_TR_MAX*1e6:.4f}] um",
        f"file spacing   = {(FILE_TR_MAX-FILE_TR_MIN)/(N_FILE-1)*1e9:.4f} nm  (finer than sim dx: {DX*1e9:.4f} nm)",
        "",
        f"amp            = {AMP:.4e} V/m",
        f"t_centre       = {T_CENTRE*1e15:.4f} fs, tau = {TAU*1e15:.4f} fs",
        f"t_peak_at_focus= {T_PEAK_AT_FOCUS*1e15:.4f} fs   (best-SNR snapshot to analyse)",
        f"t_end          = {T_END*1e15:.4f} fs, dt_snapshot = {DT_SNAPSHOT*1e15:.4f} fs",
        f"  -> n_snapshots = {int(round(T_END/DT_SNAPSHOT))}",
    ]
    return "\n".join(lines)


def w_theory(xi):
    """Beam radius at beam-coordinate xi (xi=0 at focus)."""
    return W0 * np.sqrt(1.0 + (xi / X_R) ** 2)


if __name__ == "__main__":
    print(summary())