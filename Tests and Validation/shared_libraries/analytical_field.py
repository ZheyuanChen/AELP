"""
analytical_field.py
====================
Closed-form paraxial Gaussian-beam field generator, evaluable at ANY
(x, y[, z], t) point -- not just the injection boundary. Shared between
2D and 3D campaigns (dimension-agnostic via the ``dim`` argument), which
is why it lives at the Viking_tests_new root rather than inside either
dimension's ``common/`` folder.

Why this formula, not something derived fresh
-----------------------------------------------
This is the SAME functional form (w(xi), R(xi), gouy(xi), and the
amp/phase combination) already used at the injection BOUNDARY in every
campaign's ``generate_paraxial.py`` and empirically validated in this
project (see injector_2x2_validation/physics_params.py's docstring: "the
beam measurably converges to the right focus location and width in a
real EPOCH3D run, not just reasoned about"). w(xi)/R(xi)/gouy(xi) are, by
definition in Gaussian beam theory, exactly the functions that make this
same formula a self-consistent solution of the paraxial wave equation at
EVERY x, not just the one boundary plane it was originally evaluated at
-- so extending it to general x is not a new derivation, just evaluating
the same already-trusted formula at more points.

Two bugs found and fixed by comparing against real EPOCH output
(2 July 2026 session) before trusting this module:

1. AMPLITUDE NORMALISATION. ``amp`` here means what the deck literally
   injects: the peak field AT THE BOUNDARY (x=0), since
   generate_paraxial.py's boundary profile is exp(-(y/w_bnd)^2) with
   peak exactly 1.0, scaled by the deck's ``amp = 3.2e12``. Standard
   Gaussian-beam notation instead treats the amplitude prefactor
   w0/w(xi) [3D] or sqrt(w0/w(xi)) [2D] as normalised to 1 AT FOCUS (xi=0),
   not at the boundary -- so naively using ``amp`` as that formula's E0
   systematically UNDER-predicts the field everywhere except exactly at
   the boundary (worse near focus, where the true amplitude grows
   substantially above the boundary value). Fixed by rescaling: the
   effective focus-referenced E0 = amp / amp_prefactor(xi_bnd), so the
   formula reproduces exactly ``amp`` at the injection boundary and
   correctly grows toward focus from there. Caught because a triptych
   screenshot showed the numerical (real) field visibly brighter than
   the analytical prediction specifically near focus, not spread evenly
   across the whole beam -- an amplitude-normalisation signature, not a
   phase or shape error.

2. CAUSALITY. The temporal envelope exp(-((t-x/c-t_centre)/tau)^2) is an
   idealised INFINITE-duration Gaussian -- it has a small but nonzero
   value for retarded time arbitrarily far in the past, which is
   unphysical given the deck's injection genuinely starts at t_start (no
   field exists anywhere before that). Caught because a t=0 screenshot
   showed a small nonzero analytical field at x>0 while the real
   (numerical) field was correctly exactly zero everywhere (the pulse
   hadn't had time to arrive). Fixed with a hard mask: zero wherever the
   retarded time (t - x/c) is before t_start.

Carrier convention: E(x,y,t) = E0 * profile(x,y,t) * sin(omega0*t - k0*x
+ phase(x,y)), matching the project's documented EPOCH sign convention
(physics_params.py: "E_phys=amp*profile*sin(wt+phase)"), extended so the
carrier itself also propagates forward in +x (the omega0*t - k0*x
combination is the standard forward-travelling-wave phase; at x=0 it
reduces to the boundary's own sin(omega0*t+phase(y)) convention exactly).
"""
import numpy as np

M_E = 9.1093837015e-31      # electron mass [kg]
C_LIGHT = 299792458.0       # speed of light [m/s]
E_CHARGE = 1.602176634e-19  # elementary charge [C]


def a0_norm(lambda0, c_light=C_LIGHT):
    """E-field [V/m] corresponding to a0 = 1, for wavelength lambda0."""
    omega0 = 2.0 * np.pi * c_light / lambda0
    return M_E * c_light * omega0 / E_CHARGE


def _beam_shape(xi, x_r, dim):
    """w(xi), amp_prefactor(xi) [normalised to 1 AT FOCUS, xi=0], gouy(xi)."""
    w_xi = np.sqrt(1.0 + (xi / x_r) ** 2)  # in units of w0; multiply by w0 by caller
    if dim == "2d":
        return np.sqrt(1.0 / w_xi), 0.5 * np.arctan(xi / x_r)
    elif dim == "3d":
        return 1.0 / w_xi, np.arctan(xi / x_r)
    else:
        raise ValueError(f"dim must be '2d' or '3d', got {dim!r}")


def analytical_paraxial_ey(x, y, z, t, *, w0, x_r, x_spot, tau, t_centre,
                           k0, amp, dim, t_start=0.0, c_light=C_LIGHT):
    """
    Closed-form paraxial Gaussian beam E_y(x,y,[z],t), in V/m.

    Parameters
    ----------
    x, y, z : broadcastable arrays (metres). z is ignored for dim="2d"
        (pass zeros of the right shape, or None).
    t : scalar snapshot time (seconds).
    w0, x_r, x_spot, tau, t_centre, k0 : beam parameters, matching the
        campaign's own physics_params.py names exactly (W0, X_R, X_SPOT,
        PULSE_TAU, T_CENTRE, K0).
    amp : peak field AT THE INJECTION BOUNDARY (x=0), i.e. exactly the
        deck's `amp = ...` value -- NOT the standard Gaussian-beam-theory
        focus-plane E0. This function rescales internally (see module
        docstring, bug #1).
    t_start : deck's laser `t_start` (seconds); field is exactly zero for
        retarded time before this (module docstring, bug #2).
    dim : "2d" (1D-transverse slab beam, HALF Gouy phase, amplitude ~
        1/sqrt(w)) or "3d" (circularly-symmetric, FULL Gouy phase,
        amplitude ~ 1/w).
    """
    xi = x - x_spot
    xi_bnd = 0.0 - x_spot  # injection boundary is always at x=0 in these campaigns
    w_xi_norm, gouy_xi = _beam_shape(xi, x_r, dim)
    amp_prefactor_bnd, _ = _beam_shape(np.asarray(xi_bnd), x_r, dim)
    w_xi = w0 * w_xi_norm

    e0 = amp / float(amp_prefactor_bnd)  # rescale so `amp` is matched exactly at x=0

    with np.errstate(divide="ignore", invalid="ignore"):
        r_xi = np.where(xi != 0, xi * (1.0 + (x_r / np.where(xi != 0, xi, 1.0)) ** 2), np.inf)

    rho2 = y ** 2 if dim == "2d" else y ** 2 + (0.0 if z is None else z ** 2)

    retarded_t = t - x / c_light
    causal_mask = (retarded_t >= t_start).astype(float) if hasattr(retarded_t, "astype") \
        else float(retarded_t >= t_start)

    envelope = w_xi_norm * np.exp(-rho2 / w_xi ** 2) * \
        np.exp(-((retarded_t - t_centre) / tau) ** 2) * causal_mask
    quad_phase = -k0 * rho2 / (2.0 * r_xi) - gouy_xi
    omega0 = k0 * c_light
    carrier = np.sin(omega0 * t - k0 * x + quad_phase)
    return e0 * envelope * carrier
