"""
Campaign C part 3 (2D) -- FILE-LEVEL crosscheck: the new rt-exporter's
binary pair vs test B's original hand-rolled files (which drove the
existing, validated 2d/lasy EPOCH run).

Both pipelines resample the same rt construction with the same scheme,
so the amplitude should agree to numerical noise. The phase differs by
(at most) a constant piston: test B referenced phi at (r=0, t=T_CENTRE)
while the exporter references at the envelope peak (r=0, t_peak) --
identical up to the discrete-t position of the peak. Both raw and
piston-removed phase differences are reported; the amplitude-weighted
piston-removed number is the pass metric (target: interpolation floor,
<< the ~1-2%-of-peak band of campaign C parts 1-2).

Small 2D files (6.4 MB each) -- login-node OK, no SDF involved.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import physics_params as P

OLD = (HERE.parents[2] / "campaign_B_tight_focus_f1" / "2d" / "lasy")
NEW = HERE / "lasy_exporter"
SHAPE = (P.N_T, P.N_Y_OUT_2D)

amp_old = np.fromfile(OLD / "spatial_profile.dat").reshape(SHAPE)
ph_old = np.fromfile(OLD / "phase_profile.dat").reshape(SHAPE)
amp_new = np.fromfile(NEW / "laser_amplitude.dat").reshape(SHAPE)
ph_new = np.fromfile(NEW / "laser_phase.dat").reshape(SHAPE)

d_amp = np.abs(amp_new - amp_old)
raw_dphi = ph_new - ph_old
# piston estimate at the amplitude peak (constant across the whole
# amp>1% region -- verified 1e-14-flat, 11 July 2026: the two pipelines
# reference phi at t=T_CENTRE vs at the envelope peak, giving a -1.41
# mrad CEP offset). NOT an amplitude-weighted mean: that would be
# contaminated by the wings' 2*pi branch differences.
piston = float(raw_dphi[np.unravel_index(np.argmax(amp_old), SHAPE)])
# Wrap the piston-removed residual to (-pi, pi]: the two unwrap schemes
# (skimage quality-guided vs t-then-transverse seam-free) may pick
# different 2*pi branches in the <5%-amplitude wings, and a 2*pi branch
# is physically identical (the phase enters as sin(w0 t + phi)). What
# matters is the wrapped residual, amplitude-weighted.
res = raw_dphi - piston
res_wrapped = (res + np.pi) % (2.0 * np.pi) - np.pi

# Pass metric on the physically meaningful region only (amp > 1% of
# peak): below that the envelope phase is numerical noise in BOTH
# pipelines and the two unwrap schemes legitimately disagree (2*pi
# branches, half-branch boundary pixels). Wing behaviour is reported
# for information; its field-level impact is bounded by amp * |res| --
# also printed.
core = amp_old > 0.01
d_phi_core = np.abs(res_wrapped[core]).max() if core.any() else 0.0
d_phi_weighted = np.abs(res_wrapped) * amp_old

print(f"amplitude: max|diff| = {d_amp.max():.3e} of peak "
      f"(mean {d_amp.mean():.3e})")
print(f"phase: piston (at amp peak) = {piston:+.6e} rad "
      "(reference-convention CEP offset, constant)")
print(f"phase: wrapped piston-removed max|diff| where amp>1%: "
      f"{d_phi_core:.3e} rad")
print(f"phase (info): wing amp-weighted wrapped residual max = "
      f"{d_phi_weighted.max():.3e} (bounds field impact, of peak)")
n_branch = int(((np.abs(res) > 1.0) & (amp_old > 0.0)).sum())
print(f"phase (info): pixels on a different 2*pi branch (wings only): "
      f"{n_branch}")

ok = d_amp.max() < 1e-3 and d_phi_core < 1e-3
print("FILE-LEVEL CROSSCHECK:", "PASS" if ok else "FAIL",
      "(thresholds, amp>1% region: amp 1e-3 of peak / phase 1e-3 rad)")
