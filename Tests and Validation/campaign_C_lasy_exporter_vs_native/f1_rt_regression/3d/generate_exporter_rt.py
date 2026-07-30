"""
Campaign C part 3 (3D) -- generate the f/1 test-B beam's epoch3d injection
files through the lasy-fork exporter's NEW rt path.

Construction is identical to test B's own 3d/generate_lasy.py (same
GaussianProfile, same rt npoints=(600, 500), same radial extent 5*w_bnd,
same propagate(-X_SPOT)); the resample onto the 1200^2 Cartesian (y,z)
grid via r=sqrt(y^2+z^2) now happens inside write_to_file(file_format=
"epoch") instead of the hand-rolled per-timestep interpolation loop.
carrier_phase_ref reproduces test B's -GOUY_3D (FULL Gouy) pin.

MUST run via sbatch (job_generate_3d.slurm): the two 1200x1200x500
float64 output arrays are ~5.8 GB each and the exporter transiently
holds ~4 such arrays (~25 GB).
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
              npoints=(P.N_TR_OUT_3D // 2, P.N_T), profile=profile)
print("Propagating from focus to the injection boundary (3D)...")
laser.propagate(-P.X_SPOT)

laser.write_to_file(
    file_prefix="laser",
    file_format="epoch",
    write_dir=str(HERE / "lasy_exporter"),
    carrier_phase_ref=P.CARRIER_REF_3D,
    transverse_points=P.N_TR_OUT_3D,
    transverse_extent=P.HALF_TRANSVERSE,
)
