"""
Campaign C part 3 (3D) -- FILE-LEVEL crosscheck: the new rt-exporter's
binary pair vs test B's original hand-rolled 3D files (which drove the
existing 3d/lasy EPOCH run; those originals survive on scratch).

Same metric philosophy as the 2D crosscheck (see
../2d/crosscheck_original_files.py): amplitude difference of
peak; phase compared piston-removed and wrapped to (-pi, pi], with the
pass criterion evaluated where amp > 1% of peak only (below that the
envelope phase is numerical noise in both pipelines and the two unwrap
schemes legitimately pick different 2*pi branches).

The four files are ~5.8 GB each -- processed as memory-mapped arrays in
t-slabs, so the resident footprint stays ~small. Still a multi-GB I/O
scan: run via sbatch (job_crosscheck_3d.slurm), not on the login node.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import physics_params as P

OLD = (HERE.parents[2] / "campaign_B_tight_focus_f1"
       / "3d_full_resolution" / "lasy")
NEW = HERE / "lasy_exporter"
N = P.N_TR_OUT_3D
SHAPE = (P.N_T, N, N)

amp_old = np.memmap(OLD / "spatial_profile.dat", dtype=np.float64,
                    mode="r", shape=SHAPE)
ph_old = np.memmap(OLD / "phase_profile.dat", dtype=np.float64,
                   mode="r", shape=SHAPE)
amp_new = np.memmap(NEW / "laser_amplitude.dat", dtype=np.float64,
                    mode="r", shape=SHAPE)
ph_new = np.memmap(NEW / "laser_phase.dat", dtype=np.float64,
                   mode="r", shape=SHAPE)

# piston at the global amplitude peak: on-axis, near t_centre. Find it
# in a slab scan first.
peak_val, peak_idx = -1.0, None
for it in range(SHAPE[0]):
    m = float(amp_old[it].max())
    if m > peak_val:
        peak_val = m
        peak_idx = (it,) + tuple(
            np.unravel_index(np.argmax(amp_old[it]), (N, N)))
piston = float(ph_new[peak_idx] - ph_old[peak_idx])

max_damp = 0.0
d_phi_core = 0.0
d_phi_weighted = 0.0
n_branch = 0
two_pi = 2.0 * np.pi
for it in range(SHAPE[0]):
    a_o = np.asarray(amp_old[it])
    d_amp = np.abs(np.asarray(amp_new[it]) - a_o)
    max_damp = max(max_damp, float(d_amp.max()))
    res = np.asarray(ph_new[it]) - np.asarray(ph_old[it]) - piston
    res_w = (res + np.pi) % two_pi - np.pi
    core = a_o > 0.01
    if core.any():
        d_phi_core = max(d_phi_core, float(np.abs(res_w[core]).max()))
    d_phi_weighted = max(d_phi_weighted, float((np.abs(res_w) * a_o).max()))
    n_branch += int(((np.abs(res) > 1.0) & (a_o > 0.0)).sum())

print(f"amplitude: max|diff| = {max_damp:.3e} of peak")
print(f"phase: piston (at amp peak {peak_idx}) = {piston:+.6e} rad")
print(f"phase: wrapped piston-removed max|diff| where amp>1%: "
      f"{d_phi_core:.3e} rad")
print(f"phase (info): wing amp-weighted wrapped residual max = "
      f"{d_phi_weighted:.3e} (bounds field impact, of peak)")
print(f"phase (info): pixels on a different 2*pi branch (wings only): "
      f"{n_branch}")

ok = max_damp < 1e-3 and d_phi_core < 1e-3
print("FILE-LEVEL CROSSCHECK (3D):", "PASS" if ok else "FAIL",
      "(thresholds, amp>1% region: amp 1e-3 of peak / phase 1e-3 rad)")
