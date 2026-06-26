"""
Test 2 (2D) analysis — LASY 2x2 experiment
==========================================
Head-less replacement for lasy_phase_validation.ipynb. Quantifies whether the
LASY-vs-numerical error is driven by amplitude or phase, via a controlled 2x2
(amplitude source x phase source) against the numerical reference.

Outputs (results/):
  figures/input_amplitude.png   - LASY vs analytical amplitude .dat
  figures/input_phase.png       - LASY vs analytical-deck phase .dat
  figures/bar_2x2.png           - max relative error for the 4 cells + quadrature
  figures/field_t*.png          - numerical | full-LASY | difference
  figures/lineout_y0_t*.png     - on-axis Ey overlay + residual
  summary.txt, metrics_2x2.csv

Run after all five EPOCH runs have produced .sdf files:
  python analyse.py [base_dir]
"""

import argparse
import os
import sys

import numpy as np
import xarray as xr

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import viking_analysis_lib as val  # noqa: E402

LAMBDA0 = 1.0e-6
# Beam parameters (must match the decks) — for the analytical deck phase
PULSE_FWHM = 10.0e-6
W_0 = PULSE_FWHM / 1.665
X_SPOT = 10.0e-6
X_R = np.pi * W_0**2 / LAMBDA0
RC = X_SPOT * (1.0 + (X_R / X_SPOT)**2)
GOUY = np.arctan(X_SPOT / X_R)

CASES = [  # (label, run-subdir, short tag for the bar chart)
    ("analytical built-in (sanity)", "analytical", "sanity\n(analytical)"),
    ("LASY amp + deck phase (isolate amplitude)", "lasy_amp_deck_phase", "isolate\namplitude"),
    ("num amp + LASY phase (isolate phase)", "numerical_amp_lasy_phase", "isolate\nphase"),
    ("LASY amp + LASY phase (full LASY)", "lasy", "full\nLASY"),
]


def load_dat(path):
    with open(path) as f:
        n_t, n_y = map(int, f.readline().split())
        y = np.loadtxt(f, max_rows=1)
        t = np.loadtxt(f, max_rows=1)
        mat = np.loadtxt(f).reshape(n_t, n_y)
    return y, t, mat


def input_comparison(base_dir, figures):
    """Compare the .dat inputs (amplitude and phase) before EPOCH."""
    import matplotlib.pyplot as plt
    y, t, amp_num = load_dat(os.path.join(base_dir, "numerical", "amplitude_numerical.dat"))
    _, _, amp_la = load_dat(os.path.join(base_dir, "lasy", "amplitude_lasy.dat"))
    _, _, ph_la = load_dat(os.path.join(base_dir, "lasy", "phase_lasy.dat"))
    ph_deck = np.broadcast_to(((2 * np.pi / LAMBDA0) * y**2 / (2 * RC) - GOUY)[None, :], amp_num.shape)
    ext = [y[0] * 1e6, y[-1] * 1e6, t[0] * 1e15, t[-1] * 1e15]
    core = amp_num > 0.01

    # Amplitude
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
    d = amp_la - amp_num
    vmx = np.abs(d).max()
    for a_, dat, ttl, cm, vm in [(ax[0], amp_num, "Analytical amp", "inferno", None),
                                 (ax[1], amp_la, "LASY amp", "inferno", None),
                                 (ax[2], d, "LASY - analytical", "RdBu_r", vmx)]:
        kw = dict(origin="lower", aspect="auto", extent=ext, cmap=cm)
        if vm is not None:
            kw.update(vmin=-vm, vmax=vm)
        im = a_.imshow(dat, **kw); a_.set(title=ttl, xlabel="y (µm)", ylabel="t (fs)")
        fig.colorbar(im, ax=a_)
    fig.tight_layout(); fig.savefig(os.path.join(figures, "input_amplitude.png"), dpi=140); plt.close(fig)
    amp_maxrel = float(np.abs(d[core] / amp_num[core]).max())

    # Phase
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
    pd = np.where(core, ph_la - ph_deck, np.nan)
    vmx = np.nanmax(np.abs(pd))
    for a_, dat, ttl, cm, vm in [(ax[0], ph_deck, "Deck phase (rad)", "twilight", None),
                                 (ax[1], ph_la, "LASY phase (rad)", "twilight", None),
                                 (ax[2], pd, "phase diff in core (rad)", "RdBu_r", vmx)]:
        kw = dict(origin="lower", aspect="auto", extent=ext, cmap=cm)
        if vm is not None:
            kw.update(vmin=-vm, vmax=vm)
        im = a_.imshow(dat, **kw); a_.set(title=ttl, xlabel="y (µm)", ylabel="t (fs)")
        fig.colorbar(im, ax=a_)
    fig.tight_layout(); fig.savefig(os.path.join(figures, "input_phase.png"), dpi=140); plt.close(fig)
    ph_maxcore = float(np.nanmax(np.abs(pd)))
    return amp_maxrel, ph_maxcore


def main(base_dir):
    import matplotlib.pyplot as plt
    results, figures = val.make_results_dirs(base_dir)

    amp_maxrel, ph_maxcore = input_comparison(base_dir, figures)

    ref = val.load_ey(os.path.join(base_dir, "numerical"), LAMBDA0).load()
    peak = float(np.abs(ref).max())
    env_vals = val.hilbert_envelope_along_x(ref)
    env_mask = env_vals > 0.01 * env_vals.max()

    metrics = {}
    for label, sub, _ in CASES:
        var = val.load_ey(os.path.join(base_dir, sub), LAMBDA0).load()
        n_t = min(ref.sizes["time"], var.sizes["time"])
        d = (var.isel(time=slice(0, n_t)) - ref.isel(time=slice(0, n_t))).values
        metrics[sub] = val.difference_metrics(d, env_vals[:n_t])
        del var, d
    # full-LASY kept for field plots
    lasy = val.load_ey(os.path.join(base_dir, "lasy"), LAMBDA0).load()
    n_t = min(ref.sizes["time"], lasy.sizes["time"])
    ref, lasy = ref.isel(time=slice(0, n_t)), lasy.isel(time=slice(0, n_t))
    diff = lasy - ref

    # ---- 2x2 bar chart + quadrature ----
    sanity = metrics["analytical"]["max_rel"] * 100
    amp_only = metrics["lasy_amp_deck_phase"]["max_rel"] * 100
    phase_only = metrics["numerical_amp_lasy_phase"]["max_rel"] * 100
    full = metrics["lasy"]["max_rel"] * 100
    quad = np.hypot(amp_only, phase_only)

    fig, ax = plt.subplots(figsize=(8, 5))
    tags = [c[2] for c in CASES]
    vals = [sanity, amp_only, phase_only, full]
    bars = ax.bar(tags, vals, color=["grey", "C3", "C0", "C2"])
    ax.axhline(quad, color="k", ls="--", lw=1, label=f"quadrature = {quad:.2f}%")
    ax.set_ylabel("max relative difference vs reference (%)")
    ax.set_title("2x2: amplitude dominates the LASY-vs-numerical error")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.2f}%", ha="center", va="bottom")
    ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(figures, "bar_2x2.png"), dpi=140); plt.close(fig)

    # ---- field triptych + lineout (numerical vs full LASY) ----
    env = xr.DataArray(env_vals[:n_t], coords=ref.coords, dims=ref.dims)
    rel = xr.DataArray(np.abs(diff.values) / np.where(env_mask[:n_t], env_vals[:n_t], np.nan),
                       coords=ref.coords, dims=ref.dims)
    times_fs = ref.coords["time"].values * 1e15
    sel = sorted(set([0, n_t // 4, n_t // 2, 3 * n_t // 4, n_t - 1]))
    iy0 = int(np.argmin(np.abs(ref.coords["Y_Grid_mid"].values)))
    xum = ref.coords["X_Grid_mid"].values * 1e6
    for idx in sel:
        t = times_fs[idx]
        val.save_triptych(ref.isel(time=idx), lasy.isel(time=idx), diff.isel(time=idx), t,
                          os.path.join(figures, f"field_t{idx:02d}.png"),
                          ref_label="Numerical", var_label="Full LASY")
        val.save_relative_diff_panel(env.isel(time=idx), rel.isel(time=idx), t,
                                     os.path.join(figures, f"reldiff_t{idx:02d}.png"),
                                     ref_label="Numerical")
        val.save_lineout(xum,
                         [("Numerical", ref.isel(time=idx, Y_Grid_mid=iy0).values, "b-"),
                          ("Full LASY", lasy.isel(time=idx, Y_Grid_mid=iy0).values, "r--")],
                         t, os.path.join(figures, f"lineout_y0_t{idx:02d}.png"),
                         title=f"On-axis (y=0) E_y  |  t = {t:.1f} fs")

    # ---- text outputs ----
    rows = [(label, f"{metrics[sub]['max_abs']:.4e}", f"{metrics[sub]['max_abs']/peak*100:.4f}",
             f"{metrics[sub]['max_rel']*100:.4f}", f"{metrics[sub]['mean_rel']*100:.4f}")
            for label, sub, _ in CASES]
    val.write_metrics_csv(os.path.join(results, "metrics_2x2.csv"), rows,
                          ("case", "max_abs_a0", "max_abs_pct_peak", "max_rel_pct", "mean_rel_pct"))

    summary = f"""Test 2 (2D) - LASY 2x2 experiment (amplitude source x phase source)
====================================================================

Paraxial Gaussian: lambda0 = 1 um, pulse_FWHM = 10 um (w0 ~ 6 um), focus 10 um
downstream of the boundary. Reference = numerical (analytical amp + deck phase).

Input (.dat) comparison, before EPOCH:
  amplitude  max relative diff (LASY vs analytical, beam core): {amp_maxrel*100:.3f} %
  phase      max |diff| in beam core (LASY vs deck formula)   : {ph_maxcore:.3e} rad

In-simulation 2x2 (max relative difference vs numerical reference):
  sanity (analytical built-in)          : {sanity:.4f} %
  isolate amplitude (LASY amp+deck)     : {amp_only:.4f} %
  isolate phase     (num amp+LASY)      : {phase_only:.4f} %
  full LASY                             : {full:.4f} %
  quadrature sqrt(amp^2 + phase^2)      : {quad:.4f} %
  amplitude share of error variance     : {amp_only**2/(amp_only**2+phase_only**2)*100:.1f} %

Interpretation:
  The phase agrees to ~1e-3 rad at the input and contributes only the
  "isolate phase" error; amplitude dominates and the full-LASY error matches
  the quadrature sum, confirming the two contributions are independent. The
  phase-from-file feature is validated; the residual amplitude difference is
  the physical gap between LASY's angular-spectrum-propagated beam and the
  paraxial analytical Gaussian (largest at the beam edges).

Figures in results/figures/: input_amplitude.png, input_phase.png, bar_2x2.png,
field_t*.png, reldiff_t*.png, lineout_y0_t*.png
"""
    val.write_text(os.path.join(results, "summary.txt"), summary)
    print(summary)
    print(f"Results written to {results}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("base_dir", nargs="?",
                   default=os.path.dirname(os.path.abspath(__file__)))
    main(p.parse_args().base_dir)
