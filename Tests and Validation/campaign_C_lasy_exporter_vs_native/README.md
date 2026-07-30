# Campaign C — LASY exporter vs native and legacy EPOCH injection

Campaign C validates the EPOCH exporter in the LASY fork end to end, from
profile generation and file export through EPOCH propagation. The returned
Viking runs completed on 10–11 July 2026.

**Overall verdict: PASS in 2D and 3D, with two interpretation caveats.** The
moderate-NA comparison against EPOCH's native paraxial source is not an
exporter-only test because LASY uses angular-spectrum propagation. Its 3D
peak-normalised field residual (8.52%) is also above the preliminary 1–2%
reference band. The more discriminating f/1 regression against Campaign B's
validated hand-written export is therefore the decisive test: the file cores
agree to numerical precision and the propagated fields agree within 0.47% in
2D and 0.28% on the 3D midplane.

## Results at a glance

| Comparison | Dimension | Waist-curve RMS | Max `|ΔE_y|` / reference peak | Mean relative difference | Verdict |
|---|---:|---:|---:|---:|---|
| LASY exporter vs native paraxial source | 2D | 2.500% | 3.558% | 5.739% | Pass with physics-model caveat |
| LASY exporter vs native paraxial source | 3D, `z=0` | 3.118% | 8.515% | 5.297% | Pass with physics-model caveat |
| New `rt` exporter vs Campaign B hand-written export | 2D | n/a | 0.474% | 0.156% | Pass |
| New `rt` exporter vs Campaign B hand-written export | 3D, `z=0` | n/a | 0.281% | 0.113% | Pass |

The mean-relative metric is reported for completeness but becomes unstable
where the reference envelope approaches zero. The peak-normalised absolute
difference is the robust headline field metric.

## 1. Moderate-NA baseline against EPOCH's native source

The 2D and 3D cells reuse Campaign A's `λ0 = 1 µm`, `w0 = 1.2 µm`
(`NA ≈ 0.265`) beam. `native/` uses EPOCH deck expressions; `lasy_file/`
uses the LASY exporter's normalised amplitude and phase files with a matched
carrier reference.

The fitted focus is closer to theory for the LASY-exported source in both
dimensions:

| Dimension | Cell | `x_focus` (µm) | `|x_focus − theory|` (µm) | `w0` (µm) | `w0` error |
|---|---|---:|---:|---:|---:|
| 2D | native | 8.7230 | 0.3247 | 1.2129 | 1.072% |
| 2D | LASY file | 9.0348 | 0.0130 | 1.2209 | 1.739% |
| 3D | native | 8.6440 | 0.4038 | 1.2015 | 0.125% |
| 3D | LASY file | 8.9945 | 0.0533 | 1.1997 | 0.026% |

The pointwise residuals are smooth and beam-shaped. The plots show none of
the expected exporter-failure signatures: no global sign reversal, uniform
carrier-envelope offset, or phase-unwrapping streaks. This supports the
exporter's conventions, but does not by itself isolate them from the known
angular-spectrum-versus-paraxial physics difference.

Evidence:

- [2D summary](2d/results/summary.txt), [metrics](2d/results/metrics.csv),
  [waist plot](2d/results/figures/waist_vs_x.png), and
  [field comparison](2d/results/figures/field_comparison_t07.png)
- [3D summary](3d/results/summary.txt), [metrics](3d/results/metrics.csv),
  [waist plot](3d/results/figures/waist_vs_x.png), and
  [field comparison](3d/results/figures/field_comparison_t07.png)

The pre-run local gate also passed: the exporter and Campaign B's earlier
hand-written pipeline differed by `5.8e-4` of peak in amplitude and `5.9e-4`
in amplitude-weighted field error.

## 2. f/1 `rt`-exporter regression

This test isolates the exporter by constructing the same Campaign B `rt`
beam and comparing the new LASY export path with the hand-written path that
produced Campaign B's validated EPOCH runs. Native `rt` support was added on
the LASY fork's `dev-zheyuan` branch at commit `aaacc68`.

File-level checks:

| Dimension | Amplitude max difference | Core phase difference (`amp > 1%`) | Constant piston |
|---|---:|---:|---:|
| 2D | bit-identical | `1.2e-14 rad` | `−1.412 mrad` |
| 3D | `3.1e-14` of peak | `4.5e-13 rad` | `−1.447 mrad` |

The piston comes from the two pipelines choosing adjacent temporal samples
for the phase reference. In the low-amplitude wings the two unwrapping
algorithms choose different, physically equivalent `2π` branches. Those
wings explain the short-lived propagated-field extrema above the nominal
`≈0.14%` piston imprint: 0.474% at 70 fs in 2D and 0.281% at 80 fs in 3D.
The fields otherwise track closely, and the residual falls as the pulse
leaves the domain.

Evidence:

- [2D run summary](f1_rt_regression/2d/results/run_comparison_summary.txt),
  [per-snapshot metrics](f1_rt_regression/2d/results/run_comparison_metrics.csv),
  and [80 fs comparison](f1_rt_regression/2d/results/figures/run_comparison_t08.png)
- [3D run summary](f1_rt_regression/3d/results/run_comparison_3d_summary.txt),
  [per-snapshot metrics](f1_rt_regression/3d/results/run_comparison_3d_metrics.csv),
  and [80 fs comparison](f1_rt_regression/3d/results/figures/run_comparison_3d_t08.png)

The 3D run-level comparison was Viking job `35615786`. A separate
64-rank smoke test confirmed that epoch3d's spatiotemporal path initialised
and stepped without the earlier MPI abort.

## Reproduction and retained files

- `2d/` and `3d/` retain the moderate-NA generators, decks, analysis scripts,
  machine-readable metrics, summaries, figures and animations. The committed
  2D binary pair is retained; regenerate the ~91 MB-per-file 3D pair with
  `3d/generate_lasy_exporter.py`.
- `f1_rt_regression/` retains the `rt` generators, file-level cross-checks,
  EPOCH decks, run-comparison scripts, Slurm wrappers and returned results.
  Its generated `.dat` files are deliberately ignored: the 3D pair is about
  5.8 GB per file.
- The original SDF snapshots remain on Viking and are not duplicated here.
  Field analysis must be submitted through Slurm; do not read the SDF files
  on a login node.

The copied CSVs, text summaries and figures were hash-checked against the
returned `laptop_results` bundle. The returned bundle's duplicated top-level
`results/2d` and `results/3d` export and the now-obsolete run prompt were not
retained.
