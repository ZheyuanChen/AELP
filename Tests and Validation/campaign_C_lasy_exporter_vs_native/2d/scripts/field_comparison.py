"""
Campaign C (2D) field comparison -- native-deck vs lasy-exporter E_y,
animated + screenshots, normalised to a0.

Adapted from campaign A's 2d/scripts/field_comparison.py (same
viking_analysis_lib utilities), pointed at the two campaign C cells:
  reference : native    -- EPOCH's own gauss()/quadratic-phase deck
      expression evaluation.
  variant   : lasy_file -- the lasy fork's write_to_file exporter output
      injected via the spatiotemporal file path.
Both are solved by the IDENTICAL EPOCH Maxwell solver, so the difference
measures (exporter + lasy-vs-paraxial physics) fidelity, not solver
discretisation error.

INTERPRETATION (differs from campaign A): campaign A's file/deck cells
shared one closed-form source, so its pointwise residual (~1e-3 % of
peak) was the injector floor. Campaign C's sources genuinely differ
(angular-spectrum envelope vs paraxial formula), so expect a smooth
residual of order 1-2 % of peak, amplitude-dominated, consistent with the
June 2025 lasy-vs-analytical numbers (~1.9 % amplitude, ~0.5 % phase). A
residual near 200 % of peak means a sign-flipped carrier (the -pi/2
carrier-offset bug the pytest suite guards against); a spatially uniform
few-10s-of-% residual means a CEP/carrier-pin mismatch (check
carrier_phase_ref = psi_bnd_2d in the generator against the native deck's
phase_const).

Must be run via sbatch (SDF reads are login-node-prohibited on Viking).
Usage: python field_comparison.py [base_dir]
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent.parent / "shared_libraries"))
sys.path.insert(0, str(HERE.parent.parent))
import viking_analysis_lib as val
import physics_params as P

BASE = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parent
REF_NAME, VAR_NAME = "native", "lasy_file"
REF_LABEL = "Native deck (paraxial expression)"
VAR_LABEL = "lasy exporter (file-injected)"

results, figures = val.make_results_dirs(BASE)
results, figures = Path(results), Path(figures)
animations = Path(BASE) / "results" / "animations"
animations.mkdir(parents=True, exist_ok=True)

ref = val.load_ey(str(BASE / REF_NAME), P.LAMBDA0).load()
var = val.load_ey(str(BASE / VAR_NAME), P.LAMBDA0).load()

n_t = min(ref.sizes["time"], var.sizes["time"])
ref, var = ref.isel(time=slice(0, n_t)), var.isel(time=slice(0, n_t))
diff = var - ref
peak = float(np.abs(ref).max())

env_vals = val.hilbert_envelope_along_x(ref)
times_fs = ref.coords["time"].values * 1e15

sel = sorted(set([0, n_t // 4, n_t // 2, 3 * n_t // 4, n_t - 1]))
for idx in sel:
    t = times_fs[idx]
    val.save_triptych(ref.isel(time=idx), var.isel(time=idx), diff.isel(time=idx),
                      t, figures / f"field_comparison_t{idx:02d}.png",
                      ref_label=REF_LABEL, var_label=VAR_LABEL)

anim_path = val.animate_triptych(ref, var, diff, times_fs,
                                 animations / "field_comparison.mp4",
                                 ref_label=REF_LABEL, var_label=VAR_LABEL)

rows = []
for idx in range(n_t):
    m = val.difference_metrics(diff.isel(time=idx).values, env_vals[idx])
    rows.append((f"{times_fs[idx]:.2f}", f"{m['max_abs']:.6e}",
                 f"{m['max_abs']/peak*100:.5f}", f"{m['max_rel']*100:.5f}"))
val.write_metrics_csv(Path(results) / "field_comparison_metrics.csv", rows,
                      ("time_fs", "max_abs_a0", "max_abs_pct_peak", "max_rel_pct"))

overall = val.difference_metrics(diff.values, env_vals)
summary = f"""Campaign C (2D) -- field comparison: native deck vs lasy-exporter file injection
================================================================================
Reference  : {REF_NAME} ({REF_LABEL})
Comparison : {VAR_NAME} ({VAR_LABEL})
Peak |Ey|  : {peak:.4f} a0
Snapshots  : {n_t}

Overall metrics (vs native reference):
  max |difference|         : {overall['max_abs']:.6e} a0
  max |difference| / peak  : {overall['max_abs']/peak*100:.6f} %
  max relative difference  : {overall['max_rel']*100:.6f} %
  mean relative difference : {overall['mean_rel']*100:.6f} %

Expected: ~1-2 % of peak, smooth and amplitude-dominated (the known
lasy-vs-paraxial physics gap at NA~0.265) -- NOT campaign A's ~1e-3 %
injector floor, and NOT ~200 % (sign-flipped carrier) or a uniform
few-10s-of-% offset (carrier-pin/CEP mismatch); see module docstring.

Animation: {anim_path}
Screenshots: field_comparison_t*.png in results/figures/
"""
val.write_text(Path(results) / "field_comparison_summary.txt", summary)
print(summary)
