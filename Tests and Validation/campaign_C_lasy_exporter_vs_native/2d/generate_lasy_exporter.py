"""
Campaign C (2D) -- generate the lasy_file cell's binary pair with the NEW
in-lasy exporter (`Laser.write_to_file(file_format="epoch2d")` from the
ZheyuanChen/lasy fork).

This replaces the hand-rolled conversion scripts used by campaigns A/B
(interpolation + skimage unwrap + manual .tofile) with the exporter now
built into the fork -- that exporter is exactly what this campaign
validates end-to-end. Physics is campaign A's moderate-NA converging
Gaussian (see ../physics_params.py for the design rationale and the
expected ~1.9 % lasy-vs-paraxial residual).

Pipeline: define the Gaussian AT ITS FOCUS (waist w0, peak at t=T_CENTRE
-- lasy's GaussianProfile tau and EPOCH's gauss(time, t0, w) share the
same 1/e-of-field convention, so TAU carries over unchanged), then
back-propagate by -X_SPOT so the grid holds the converging boundary
field, then export. The exporter slices lasy-y = 0, normalises to the
slice peak, unwraps/references the phase, and writes the Fortran
(n_transverse, n_t) layout epoch2d reads.

Carrier pin: carrier_phase_ref = PSI_BND_2D (half Gouy at the boundary)
so the stored on-axis peak phase is pi/2 + psi_bnd_2d -- identical to the
native deck's phase_const. Without this the two runs would differ by a
constant CEP and the pointwise field comparison would report a spurious
2|sin(delta/2)| error (the campaign-A-era 137 % lesson).

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
    npoints=(P.N_TRANS_FILE_2D, P.N_Y_LASY_2D, P.N_T_FILE),
    profile=profile,
)
print(f"Back-propagating {P.X_SPOT*1e6:.3f} um from focus to the "
      "injection boundary (xyt)...")
laser.propagate(-P.X_SPOT)

laser.write_to_file(
    file_prefix="laser",
    file_format="epoch2d",
    write_dir=str(HERE / "lasy_file"),
    carrier_phase_ref=P.PSI_BND_2D,
)

print("\nDeck cross-check (lasy_file/input.deck laser block should read):")
print(f"  n_y = {P.N_TRANS_FILE_2D}\n  n_t = {P.N_T_FILE}")
print(f"  y_min = {-P.Y_HALF:.15e}\n  y_max = {P.Y_HALF:.15e}")
print(f"  t_start = {P.T_START_FILE:.6e}\n  t_end = {P.T_END_FILE:.6e}")
print(f"  amp = {P.AMP:.4e}  (native cell's amp -- NOT the lasy peak "
      "field, the 1 J lasy energy is arbitrary and normalised out)")
print(f"  carrier pin used: psi_bnd_2d = {P.PSI_BND_2D:.15f} rad "
      f"(native phase_const = {P.PHASE_CONST_2D:.15f})")
