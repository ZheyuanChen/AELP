"""
Campaign C (3D) -- generate the lasy_file cell's binary pair with the NEW
in-lasy exporter (`Laser.write_to_file(file_format="epoch")`, the epoch3d
spatiotemporal path). See ../2d/generate_lasy_exporter.py for the shared
pipeline commentary; the only differences here are:

  - the full (n_t, n_tr2, n_tr1) 3D array is exported (no y-slice), on a
    151 x 151 transverse grid (~91 MB per file -- deliberately kept out
    of git; regenerate with this script, locally or on Viking via sbatch);
  - the carrier pin is the FULL Gouy phase (PSI_BND_3D), matching the 3D
    native deck's phase_const = pi/2 + psi_bnd_3d (circularly-symmetric
    r_yz beam -- see campaign A's physics_params.py docstring for the
    half-vs-full Gouy derivation).

Run from this directory:  python generate_lasy_exporter.py
Writes lasy_file/laser_amplitude.dat, laser_phase.dat, laser_metadata.txt.
"""
import sys
from pathlib import Path

from lasy.laser import Laser
from lasy.profiles.gaussian_profile import GaussianProfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import physics_params as P

profile = GaussianProfile(P.LAMBDA0, (1, 0), P.LASER_ENERGY, P.W0,
                          P.TAU, P.T_CENTRE)
laser = Laser(
    dim="xyt",
    lo=(-P.Y_HALF, -P.Y_HALF, P.T_START_FILE),
    hi=(P.Y_HALF, P.Y_HALF, P.T_END_FILE),
    npoints=(P.N_TR_FILE_3D, P.N_TR_FILE_3D, P.N_T_FILE),
    profile=profile,
)
print(f"Back-propagating {P.X_SPOT*1e6:.3f} um from focus to the "
      "injection boundary (xyt)...")
laser.propagate(-P.X_SPOT)

laser.write_to_file(
    file_prefix="laser",
    file_format="epoch",
    write_dir=str(HERE / "lasy_file"),
    carrier_phase_ref=P.PSI_BND_3D,
)

print("\nDeck cross-check (lasy_file/input.deck laser block should read):")
print(f"  n_tr1 = {P.N_TR_FILE_3D}\n  n_tr2 = {P.N_TR_FILE_3D}"
      f"\n  n_t = {P.N_T_FILE}")
print(f"  tr1_min/tr2_min = {-P.Y_HALF:.15e}")
print(f"  tr1_max/tr2_max = {P.Y_HALF:.15e}")
print(f"  t_start = {P.T_START_FILE:.6e}\n  t_end = {P.T_END_FILE:.6e}")
print(f"  amp = {P.AMP:.4e}  (native cell's amp)")
print(f"  carrier pin used: psi_bnd_3d = {P.PSI_BND_3D:.15f} rad "
      f"(native phase_const = {P.PHASE_CONST_3D:.15f})")
