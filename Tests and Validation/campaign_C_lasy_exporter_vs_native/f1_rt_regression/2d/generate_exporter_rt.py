"""
Campaign C part 3 (2D) -- generate the f/1 test-B beam's EPOCH injection
files through the lasy-fork exporter's NEW rt path.

Construction is identical to test B's own 2d/generate_lasy.py (same
GaussianProfile, same rt npoints=(800, 500), same radial extent 5*w_bnd,
same propagate(-X_SPOT)); the ONLY difference is that the resample-to-
Cartesian + unwrap + piston-removal + file-writing now happens inside
laser.write_to_file(file_format="epoch2d") instead of ~40 lines of
hand-rolled script. carrier_phase_ref reproduces test B's -GOUY_2D pin
(see physics_params.CARRIER_REF_2D).

Single small rt profile, seconds of numpy/lasy -- login-node OK.
"""
import sys
from pathlib import Path

from lasy.laser import Laser
from lasy.profiles.gaussian_profile import GaussianProfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import physics_params as P

profile = GaussianProfile(P.LAMBDA0, (1, 0), P.LASER_ENERGY, P.W0,
                          P.PULSE_TAU, P.T_CENTRE)
laser = Laser(dim="rt", lo=(0, P.T_START),
              hi=(P.R_MAX_FACTOR * P.W_BND, P.T_END),
              npoints=(P.N_Y_OUT_2D // 2, P.N_T), profile=profile)
print("Propagating from focus to the injection boundary (2D)...")
laser.propagate(-P.X_SPOT)

laser.write_to_file(
    file_prefix="laser",
    file_format="epoch2d",
    write_dir=str(HERE / "lasy_exporter"),
    carrier_phase_ref=P.CARRIER_REF_2D,
    transverse_points=P.N_Y_OUT_2D,
    transverse_extent=P.HALF_TRANSVERSE,
)
