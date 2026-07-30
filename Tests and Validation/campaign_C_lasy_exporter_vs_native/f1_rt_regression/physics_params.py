"""
Campaign C, part 3 (2D) -- tight-focus f/1 (test B) beam via the lasy-fork
exporter's NEW rt-geometry path, compared against test B's original
standalone `rt` generator result (`campaign_B_tight_focus_f1/2d/lasy/`).

ALL beam constants are imported from test B's own physics_params.py so
the two pipelines stay directly comparable and the constants remain
single-sourced (same pattern as campaign C's own physics_params.py
re-exporting campaign A's constants).

Why this is a separate part of campaign C rather than reusing
2d/physics_params.py: test B is a materially different, much tighter
beam (f/1, NA=0.5, w0=0.64um... actually w0=LAMBDA0/(pi*NA)=0.6366um,
tau=15fs) than campaign C's own beam (NA=0.265, w0=1.2um, tau=80fs) --
this project has already hit a transverse-resolution interpolation-
artifact bug on THIS specific beam (test B's 3D generate_lasy.py:
N_TR=400 was too coarse, fixed at 1200), so campaign C's own (coarser,
NA=0.265-calibrated) exporter resolution recipe cannot be assumed
adequate here without checking.

HISTORY NOTE (11 July 2026): this part was originally BLOCKED because
the exporter rejected dim="rt" outright and its suggested workaround
(to_cartesian()) did not exist -- the only route was a from-scratch
dim="xyt" reconstruction at 1200^2 transverse points, ~1800x campaign
C's own 2D grid. That block is now lifted: the fork's exporter gained
native rt support (fork commit adding _rt_to_cartesian; same r->|y| /
r->sqrt(y^2+z^2) resampling test B's own generators validated, with the
seam-free unwrap done in (r,t) and the piston reference taken exactly
on axis). This module therefore now drives the SAME rt construction
test B used -- npoints, radial extent, propagate distance all identical
-- through the new exporter path, making the comparison an isolated
test of the export machinery rather than of a reconstruction detour.
"""
import importlib.util
from pathlib import Path

_TEST_B = (Path(__file__).resolve().parent.parent.parent
           / "campaign_B_tight_focus_f1" / "physics_params.py")
_spec = importlib.util.spec_from_file_location("_test_b_physics", str(_TEST_B))
_B = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_B)
globals().update({k: v for k, v in vars(_B).items()
                  if not k.startswith("_") and k != "np"})

# --- rt construction, matching test B's own generators EXACTLY -------
# 2D: lasy npoints=(N_Y_OUT//2, N_T), r up to 5*W_BND, output N_Y_OUT
#     Cartesian points spanning +-HALF_TRANSVERSE (= 4*W_BND).
# 3D: lasy npoints=(N_TR_OUT//2, N_T), same radial extent, output
#     N_TR_OUT^2 -- N_TR_OUT=1200 is the value test B's own history
#     found necessary at this NA (400 gave interpolation artifacts).
N_Y_OUT_2D = 1600
N_TR_OUT_3D = 1200
R_MAX_FACTOR = 5.0        # lasy radial grid extent, in units of W_BND

LASER_ENERGY = 1.0  # lasy is more comfortable in normalised energy;
                    # amplitude is normalised out by the exporter, the
                    # physical scale re-enters via the deck's amp
                    # (test B's own 3.2e12, unchanged).

# Carrier pin: test B's generators wrote phase = -(phi - phi_ref) - GOUY;
# the exporter writes -(phi - phi_ref) + pi/2 + carrier_phase_ref. So to
# byte-match test B's convention the run must pass
#   carrier_phase_ref = -GOUY_{2D|3D} - pi/2.
import numpy as _np
CARRIER_REF_2D = -GOUY_2D - 0.5 * _np.pi   # noqa: F821 (from test B import)
CARRIER_REF_3D = -GOUY_3D - 0.5 * _np.pi   # noqa: F821


def summary():
    return "\n".join([
        _B.summary(),
        f"rt construction (2D): n_r={N_Y_OUT_2D // 2}, n_t={N_T}, "        # noqa: F821
        f"r_max={R_MAX_FACTOR}*w_bnd -> {N_Y_OUT_2D} Cartesian points",
        f"rt construction (3D): n_r={N_TR_OUT_3D // 2}, n_t={N_T}, "       # noqa: F821
        f"-> {N_TR_OUT_3D}^2 Cartesian points",
        f"carrier_phase_ref: 2D {CARRIER_REF_2D:+.6f}, 3D {CARRIER_REF_3D:+.6f}",
    ])


if __name__ == "__main__":
    print(summary())
