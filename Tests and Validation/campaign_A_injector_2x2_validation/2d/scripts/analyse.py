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

    # Bracket the ideal focal time (t_peak@focus = t_centre + x_spot/c)
    # between the two nearest AVAILABLE snapshots and linearly interpolate
    # the FITTED w(x) curve to that exact time, instead of snapping to
    # whichever snapshot happens to be nearest. Interpolating the raw
    # field itself would be unsound here -- the carrier oscillates several
    # times per dt_snapshot, so a naive linear interpolation between
    # distant-in-phase snapshots wouldn't reconstruct the intermediate
    # field correctly. The fitted width w(x), by contrast, varies slowly
    # over the ~tens-of-fs pulse duration and interpolates cleanly.
    it_lo = int(np.clip(np.searchsorted(times, P.T_PEAK_AT_FOCUS) - 1, 0, len(times) - 2))
    it_hi = it_lo + 1
    t_lo, t_hi = float(times[it_lo]), float(times[it_hi])
    frac = (P.T_PEAK_AT_FOCUS - t_lo) / (t_hi - t_lo)

    env_lo = np.abs(hilbert(ey.isel(time0=it_lo).values, axis=0))
    env_hi = np.abs(hilbert(ey.isel(time0=it_hi).values, axis=0))
    w_lo = waist_scan_1d(env_lo, x, y)
    w_hi = waist_scan_1d(env_hi, x, y)
    w_interp = (1.0 - frac) * w_lo + frac * w_hi
    x_focus, w0 = find_minimum(x, w_interp)

    # Nearest-snapshot result too, for comparison (this was the previous
    # method's number, kept so the size of the quantisation error is visible).
    nearest_is_lo = abs(t_lo - P.T_PEAK_AT_FOCUS) <= abs(t_hi - P.T_PEAK_AT_FOCUS)
    t_nearest = t_lo if nearest_is_lo else t_hi
    x_focus_nearest, w0_nearest = find_minimum(x, w_lo if nearest_is_lo else w_hi)

    return dict(x=x, w=w_interp, x_focus=x_focus, w0=w0, t=P.T_PEAK_AT_FOCUS,
               t_lo=t_lo, t_hi=t_hi, frac=frac,
               x_focus_nearest=x_focus_nearest, w0_nearest=w0_nearest,
               t_nearest=t_nearest)


cells = {name: analyse_cell(name) for name in
        ("amp_deck_phase_deck", "amp_file_phase_deck", "amp_file_phase_file")}

x_theory = np.linspace(0, 2 * P.X_SPOT, 400)
w_theory = P.w_of_xi(x_theory - P.X_SPOT)

fig, ax = plt.subplots(figsize=(9, 5.5))
ax.plot(x_theory * 1e6, w_theory * 1e6, "k-", lw=1.2, label="theory")
for (name, c), style in zip(cells.items(), ["C0o-", "C1s--", "C2^:"]):
    ax.plot(c["x"] * 1e6, c["w"] * 1e6, style, lw=1.0, ms=1.5, label=name)
ax.axvline(P.X_SPOT * 1e6, color="grey", ls=":")
ax.set(xlabel="x (um)", ylabel="beam radius w (um)",
      title="Test A (2D): injector isolation chain")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIGURES / "waist_vs_x.png", dpi=140)
plt.close(fig)

rows = [",".join(["cell", "x_focus_um", "w0_um", "x_focus_err_um", "w0_rel_err",
                  "x_focus_nearest_snapshot_um", "w0_nearest_snapshot_um"])]
for name, c in cells.items():
    rows.append(",".join(str(v) for v in (
        name, c["x_focus"] * 1e6, c["w0"] * 1e6,
        abs(c["x_focus"] - P.X_SPOT) * 1e6, abs(c["w0"] - P.W0) / P.W0,
        c["x_focus_nearest"] * 1e6, c["w0_nearest"] * 1e6)))
(RESULTS / "metrics.csv").write_text("\n".join(rows) + "\n")

base = cells["amp_deck_phase_deck"]
amp_only = cells["amp_file_phase_deck"]
both = cells["amp_file_phase_file"]
valid = np.isfinite(base["w"]) & np.isfinite(amp_only["w"]) & np.isfinite(both["w"])
amp_rms = float(np.sqrt(np.nanmean(((amp_only["w"][valid]-base["w"][valid])/base["w"][valid])**2)))
total_rms = float(np.sqrt(np.nanmean(((both["w"][valid]-base["w"][valid])/base["w"][valid])**2)))
phase_rms = float(np.sqrt(np.nanmean(((both["w"][valid]-amp_only["w"][valid])/amp_only["w"][valid])**2)))

# The plot above compares each cell against the IDEALISED THEORY curve, which
# every cell (including the native-deck baseline) sits ~0.3-0.4um off of --
# a finite-resolution/numerical-dispersion effect common to the whole chain,
# not an injector-pipeline error. That's not what this test actually checks
# for (injector fidelity = cells matching EACH OTHER, not matching theory
# exactly, which no finite-grid PIC run does). This second plot makes the
# real pass/fail metric visible directly: relative w(x) difference between
# each file-injected cell and the native-deck baseline, which the RMS
# numbers above show is ~1e-4 to 1e-3 % -- i.e. the curves overlap almost
# exactly once you're looking at the right quantity.
fig, ax = plt.subplots(figsize=(9, 4.5))
x_um = base["x"][valid] * 1e6
ax.plot(x_um, (amp_only["w"][valid] - base["w"][valid]) / base["w"][valid] * 100,
       "C1-", lw=1.0, label="amp_file_phase_deck vs baseline")
ax.plot(x_um, (both["w"][valid] - base["w"][valid]) / base["w"][valid] * 100,
       "C2-", lw=1.0, label="amp_file_phase_file vs baseline")
ax.axhline(0, color="grey", ls=":", lw=1)
ax.set(xlabel="x (um)", ylabel="relative w(x) difference vs baseline (%)",
      title="Test A (2D): injector-pipeline error (the actual pass/fail metric)")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIGURES / "injector_error_vs_x.png", dpi=140)
plt.close(fig)

any_cell = next(iter(cells.values()))
summary = f"""Test A (2D) -- injector 2x2 (3-cell) isolation
================================================
Beam: lambda0={P.LAMBDA0*1e6:.3f}um w0={P.W0*1e6:.3f}um NA~{P.NA:.4f}
(1D-transverse slab, HALF Gouy phase -- see physics_params.py).

Snapshot timing: ideal t_peak@focus = {P.T_PEAK_AT_FOCUS*1e15:.3f} fs, but
output only exists every dt_snapshot -- bracketing snapshots at
{any_cell['t_lo']*1e15:.1f} fs and {any_cell['t_hi']*1e15:.1f} fs
(interpolation fraction {any_cell['frac']:.3f}). x_focus/w0 below use the
fitted w(x) curve LINEARLY INTERPOLATED between those two snapshots to
the exact ideal time (not the raw field, which oscillates too fast per
dt_snapshot for linear interpolation to be valid -- see analyse_cell's
docstring comment). Old method (snap to nearest snapshot,
t={any_cell['t_nearest']*1e15:.1f} fs) numbers included alongside for
comparison, to show how much this quantisation correction moves things.

x_focus / w0 per cell (theory: x_focus={P.X_SPOT*1e6:.4f}um, w0={P.W0*1e6:.4f}um):
"""
for name, c in cells.items():
    summary += (f"  {name:22s}: x_focus={c['x_focus']*1e6:10.6f}um  w0={c['w0']*1e6:9.6f}um"
               f"   (nearest-snapshot was: x_focus={c['x_focus_nearest']*1e6:10.6f}um  "
               f"w0={c['w0_nearest']*1e6:9.6f}um)\n")
summary += f"""
w(x) RMS relative error:
  amplitude-injector-only    = {amp_rms*100:.3e}%
  phase-injector-additional  = {phase_rms*100:.3e}%
  total                      = {total_rms*100:.3e}%
"""
(RESULTS / "summary.txt").write_text(summary)
print(summary)
