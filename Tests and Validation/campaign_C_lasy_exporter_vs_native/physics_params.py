"""
Campaign C -- the in-lasy EPOCH exporter vs EPOCH's native laser injector.

Shared physics for the 2D and 3D cells. ALL beam/timing constants are
imported from campaign A (campaign_A_injector_2x2_validation/
physics_params.py) so the two campaigns stay directly comparable and the
constants remain single-sourced -- campaign C deliberately reuses campaign
A's moderate-NA converging Gaussian (NA~0.265, well inside paraxial
validity) so that the native deck expression is close to ground truth and
any large disagreement is attributable to the exporter pipeline, not to a
paraxial-approximation gap.

What campaign C tests that campaign A did not: campaign A proved the
FILE-INJECTION path itself is faithful (file vs deck, identical closed-form
source on both sides -- residual ~1e-3 %). Campaign C now varies the
SOURCE GENERATOR: the binary files come from the lasy fork's new
`Laser.write_to_file(file_format="epoch"/"epoch2d")` exporter (angular-
spectrum-propagated envelope, exporter-internal normalisation/phase
referencing/unwrapping), while the reference run uses the native deck
Gaussian. The expected residual is therefore NOT campaign A's ~1e-3 %:
it is the known lasy-vs-paraxial physics difference at this NA
(~1.9 % amplitude-dominated, ~0.5 % phase, June 2025 2x2 validation),
plus whatever the exporter itself gets wrong. Pass = residuals at that
~2 % level with the same amplitude-dominated structure; fail = larger, or
structured like a carrier-phase/CEP error (which would indicate an
exporter convention bug).

Carrier-phase pin: the native decks use phase_const = pi/2 + psi_bnd
(campaign A convention, dimension-correct Gouy). The exporter pins the
stored phase at the amplitude peak to pi/2 + carrier_phase_ref, so the
generators pass carrier_phase_ref = PSI_BND_2D / PSI_BND_3D to match the
native runs' carrier-envelope reference exactly -- EPOCH ignores the
deck's own `phase` expression when use_phase_from_file = T, so this pin
can only live in the file.
"""
import importlib.util
from pathlib import Path

# Load campaign A's physics_params under a private module name (both files
# are called physics_params.py, so a plain sys.path import would collide
# with THIS module in sys.modules) and re-export everything public.
_CAMPAIGN_A = (Path(__file__).resolve().parent.parent
               / "campaign_A_injector_2x2_validation" / "physics_params.py")
_spec = importlib.util.spec_from_file_location("_campaign_a_physics",
                                               str(_CAMPAIGN_A))
_A = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_A)
globals().update({k: v for k, v in vars(_A).items()
                  if not k.startswith("_") and k != "np"})

import numpy as np  # noqa: E402  (after the re-export on purpose)

# --- Campaign C additions: file/lasy grid parameters ---------------------

# Transverse half-width of the EPOCH domain and of the exported files
# (identical to campaign A's deck y bounds: 3.5 x w_bnd).
Y_HALF = 3.5 * W_BND  # noqa: F821  (re-exported from campaign A)

# lasy is more comfortable in normalised energy; the amplitude is
# normalised out by the exporter and the physical scale re-enters via the
# deck's amp (= campaign A's AMP for both cells).
LASER_ENERGY = 1.0

# File time window = the whole run (matches the native cell, whose
# t_profile keeps injecting until t_end). lasy's propagate() works in the
# co-moving frame, so the envelope peak stays at T_CENTRE in the local
# time axis after back-propagation -- the same t axis the deck sees
# (validated by campaign B's lasy cell, which used exactly this timing).
T_START_FILE = 0.0
T_END_FILE = T_END  # noqa: F821   350 fs
N_T_FILE = 501      # 0.7 fs sampling of an 80 fs envelope -- ample for
                    # EPOCH's bilinear interpolation in t

# 2D cell: the exported file transverse axis is lasy's x (601 points over
# +/-Y_HALF ~ 31 nm spacing, oversampling the 451-cell EPOCH grid); lasy's
# y axis exists only so the xyt propagation of the sliced 3D beam is
# accurate (97 points ~ 196 nm ~ w0/6 at the focus -- adequate for the
# angular-spectrum step; the slice at y=0 is what gets exported).
N_TRANS_FILE_2D = 601
N_Y_LASY_2D = 97

# 3D cell: 151 points per transverse axis (~125 nm) keeps the file pair at
# ~91 MB each (n_tr^2 * n_t * 8 bytes). This is coarser than the 451-cell
# EPOCH grid, but bilinear interpolation of a w_bnd = 2.68 um Gaussian
# sampled every 125 nm has O((h/w)^2) ~ 0.2 % error -- an order below the
# ~2 % physics scale this campaign resolves.
N_TR_FILE_3D = 151


def summary():
    return "\n".join([
        _A.summary(),
        f"y_half={Y_HALF*1e6:.4f}um (=3.5 w_bnd, campaign A deck bounds)",
        f"file grid 2D: {N_TRANS_FILE_2D} (transverse) x {N_T_FILE} (t), "
        f"lasy y sampling {N_Y_LASY_2D}",
        f"file grid 3D: {N_TR_FILE_3D}^2 (transverse) x {N_T_FILE} (t)",
        f"file window: [{T_START_FILE*1e15:.0f}, {T_END_FILE*1e15:.0f}] fs",
    ])


if __name__ == "__main__":
    print(summary())
