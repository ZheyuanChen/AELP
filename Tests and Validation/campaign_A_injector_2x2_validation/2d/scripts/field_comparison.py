"""
Test A (2D) field comparison -- analytical vs numerical E_y, animated +
screenshots, normalised to a0.

Unlike the waist-scan analysis (analyse.py), which only looks at the
envelope at one snapshot near focus, this compares the FULL oscillating
E_y field across all snapshots between two EPOCH runs:
  reference ("Analytical") : amp_deck_phase_deck -- EPOCH's own native
      gauss()/quadratic-phase deck-expression evaluation.
  variant  ("Numerical")   : amp_file_phase_file -- the full file-
      injection pipeline (amplitude AND phase from binary files).
Both are solved by the IDENTICAL EPOCH Maxwell solver -- the only
difference between them is how the boundary condition was specified, so
any difference measures injector-pipeline fidelity directly, not solver
discretisation error (which affects both equally). This mirrors the
established pattern from Viking_tests_new/2D/test1_gaussian_2d/analyse.py
(same a0_norm/load_ey/save_triptych utilities), just pointed at this
campaign's actual 3-cell-chain directories instead of that test's own.

Must be run via sbatch (SDF reads are login-node-prohibited on Viking).
Usage: python field_comparison.py [base_dir]
"""
import sys
from pathlib import Path

import numpy as np
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "2D" / "common"))
sys.path.insert(0, str(HERE.parent))
import viking_analysis_lib as val
import physics_params as P

BASE = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
REF_NAME, VAR_NAME = "amp_deck_phase_deck", "amp_file_phase_file"
REF_LABEL, VAR_LABEL = "Analytical (native deck)", "Numerical (file-injected)"

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
summary = f"""Test A (2D) -- field comparison: analytical (native deck) vs numerical (file-injected)
========================================================================================
Reference  : {REF_NAME} ({REF_LABEL})
Comparison : {VAR_NAME} ({VAR_LABEL})
Peak |Ey|  : {peak:.4f} a0
Snapshots  : {n_t}

Overall metrics (vs analytical reference):
  max |difference|         : {overall['max_abs']:.6e} a0
  max |difference| / peak  : {overall['max_abs']/peak*100:.6f} %
  max relative difference  : {overall['max_rel']*100:.6f} %
  mean relative difference : {overall['mean_rel']*100:.6f} %

Animation: {anim_path}
Screenshots: field_comparison_t*.png in results/figures/
"""
val.write_text(Path(results) / "field_comparison_summary.txt", summary)
print(summary)
