"""
Test A (2D) analysis -- 3-cell injector isolation chain (see 3d/analyse.py
for the full design rationale; this is the 1D-transverse analogue).

Usage: python analyse.py [base_dir]
"""
import sys
import glob
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sdf_xarray as sdfxr
from scipy.signal import hilbert
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import physics_params as P

BASE = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
RESULTS = BASE / "results"
FIGURES = RESULTS / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)


def gauss1d(y, amp, y0, w):
    return amp * np.exp(-((y - y0) / w) ** 2)


def waist_scan_1d(env_xy, x, y, snr_frac=0.05):
    peak = float(np.abs(env_xy).max())
    w = np.full(len(x), np.nan)
    for ix in range(len(x)):
        slab = env_xy[ix]
        if slab.max() < snr_frac * peak:
            continue
        try:
            popt, _ = curve_fit(gauss1d, y, slab,
                                p0=[slab.max(), 0.0, 1e-6], maxfev=5000)
            w[ix] = abs(popt[2])
        except Exception:
            pass
    return w


def find_minimum(x, w):
    valid = np.isfinite(w)
    if valid.sum() < 3:
        return np.nan, np.nan
    xi, wi = x[valid], w[valid]
    i = int(np.argmin(wi))
    if i == 0 or i == len(wi) - 1:
        return float(xi[i]), float(wi[i])
    x0, x1, x2 = xi[i-1], xi[i], xi[i+1]
    y0, y1, y2 = wi[i-1], wi[i], wi[i+1]
    denom = (x0-x1)*(x0-x2)*(x2-x1)
    if denom == 0:
        return float(x1), float(y1)
    a = (x2*(y1-y0) + x1*(y0-y2) + x0*(y2-y1)) / denom
    b = (x2**2*(y0-y1) + x1**2*(y2-y0) + x0**2*(y1-y2)) / denom
    if a == 0:
        return float(x1), float(y1)
    xf = -b/(2*a)
    c = y0 - a*x0**2 - b*x0
    return float(xf), float(a*xf**2 + b*xf + c)


def analyse_cell(name):
    d = BASE / name
    files = sorted(glob.glob(str(d / "*.sdf")))
    ds = sdfxr.open_mfdataset(files, separate_times=True,
                              data_vars=["Electric_Field_Ey"])
    ey = ds["Electric_Field_Ey"]
    x = ey["X_Grid_mid"].values
    y = ey["Y_Grid_mid"].values
    times = ey["time0"].values
    it = int(np.argmin(np.abs(times - P.T_PEAK_AT_FOCUS)))
    field = ey.isel(time0=it).values
    env = np.abs(hilbert(field, axis=0))
    w = waist_scan_1d(env, x, y)
    x_focus, w0 = find_minimum(x, w)
    return dict(x=x, w=w, x_focus=x_focus, w0=w0, t=float(times[it]))


cells = {name: analyse_cell(name) for name in
        ("amp_deck_phase_deck", "amp_file_phase_deck", "amp_file_phase_file")}

x_theory = np.linspace(0, 2 * P.X_SPOT, 400)
w_theory = P.w_of_xi(x_theory - P.X_SPOT)

fig, ax = plt.subplots(figsize=(9, 5.5))
ax.plot(x_theory * 1e6, w_theory * 1e6, "k-", lw=1.6, label="theory")
for (name, c), style in zip(cells.items(), ["C0o-", "C1s--", "C2^:"]):
    ax.plot(c["x"] * 1e6, c["w"] * 1e6, style, lw=1.4, label=name)
ax.axvline(P.X_SPOT * 1e6, color="grey", ls=":")
ax.set(xlabel="x (um)", ylabel="beam radius w (um)",
      title="Test A (2D): injector isolation chain")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIGURES / "waist_vs_x.png", dpi=140)
plt.close(fig)

rows = [",".join(["cell", "x_focus_um", "w0_um", "x_focus_err_um", "w0_rel_err"])]
for name, c in cells.items():
    rows.append(",".join(str(v) for v in (
        name, c["x_focus"] * 1e6, c["w0"] * 1e6,
        abs(c["x_focus"] - P.X_SPOT) * 1e6, abs(c["w0"] - P.W0) / P.W0)))
(RESULTS / "metrics.csv").write_text("\n".join(rows) + "\n")

base = cells["amp_deck_phase_deck"]
amp_only = cells["amp_file_phase_deck"]
both = cells["amp_file_phase_file"]
valid = np.isfinite(base["w"]) & np.isfinite(amp_only["w"]) & np.isfinite(both["w"])
amp_rms = float(np.sqrt(np.nanmean(((amp_only["w"][valid]-base["w"][valid])/base["w"][valid])**2)))
total_rms = float(np.sqrt(np.nanmean(((both["w"][valid]-base["w"][valid])/base["w"][valid])**2)))
phase_rms = float(np.sqrt(np.nanmean(((both["w"][valid]-amp_only["w"][valid])/amp_only["w"][valid])**2)))

summary = f"""Test A (2D) -- injector 2x2 (3-cell) isolation
================================================
Beam: lambda0={P.LAMBDA0*1e6:.3f}um w0={P.W0*1e6:.3f}um NA~{P.NA:.4f}
(1D-transverse slab, HALF Gouy phase -- see physics_params.py).

x_focus / w0 per cell (theory: x_focus={P.X_SPOT*1e6:.4f}um, w0={P.W0*1e6:.4f}um):
"""
for name, c in cells.items():
    summary += f"  {name:22s}: x_focus={c['x_focus']*1e6:8.4f}um  w0={c['w0']*1e6:7.4f}um\n"
summary += f"""
w(x) RMS relative error:
  amplitude-injector-only    = {amp_rms*100:.3f}%
  phase-injector-additional  = {phase_rms*100:.3f}%
  total                      = {total_rms*100:.3f}%
"""
(RESULTS / "summary.txt").write_text(summary)
print(summary)
