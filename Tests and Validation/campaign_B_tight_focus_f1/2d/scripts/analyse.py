"""
Test B (2D) analysis -- LASY (non-paraxial) vs paraxial closed-form, f/1
tight-focusing demo. Unlike Campaign A's injector-isolation chain, there is
no native-deck "ground truth" cell here (both paraxial/ and lasy/ inject
via file only) -- the comparison is against each other and against the
paraxial theory curve, to see whether genuine non-paraxial (angular-
spectrum) propagation differs measurably from the closed-form formula at
NA=0.5 where the paraxial approximation is expected to break down.

Same 1D-transverse waist-scan method as injector_2x2_validation/2d/
analyse.py (Hilbert envelope along x at the snapshot nearest peak-at-focus,
1D Gaussian fit per x-plane, parabolic interpolation to the minimum).

Usage: python analyse.py [base_dir]
"""
import sys
import glob
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

T_PEAK_AT_FOCUS = P.T_CENTRE + P.X_SPOT / P.C_LIGHT


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

    # Bracket the ideal focal time between the two nearest AVAILABLE
    # snapshots and linearly interpolate the FITTED w(x) curve to that
    # exact time, instead of snapping to whichever snapshot happens to be
    # nearest. Interpolating the raw field itself would be unsound -- the
    # carrier oscillates several times per dt_snapshot, so linear
    # interpolation between distant-in-phase snapshots wouldn't
    # reconstruct the intermediate field correctly. The fitted width
    # w(x), by contrast, varies slowly over the pulse duration and
    # interpolates cleanly.
    it_lo = int(np.clip(np.searchsorted(times, T_PEAK_AT_FOCUS) - 1, 0, len(times) - 2))
    it_hi = it_lo + 1
    t_lo, t_hi = float(times[it_lo]), float(times[it_hi])
    frac = (T_PEAK_AT_FOCUS - t_lo) / (t_hi - t_lo)

    env_lo = np.abs(hilbert(ey.isel(time0=it_lo).values, axis=0))
    env_hi = np.abs(hilbert(ey.isel(time0=it_hi).values, axis=0))
    w_lo = waist_scan_1d(env_lo, x, y)
    w_hi = waist_scan_1d(env_hi, x, y)
    w_interp = (1.0 - frac) * w_lo + frac * w_hi
    x_focus, w0 = find_minimum(x, w_interp)

    # Nearest-snapshot result too, for comparison (the previous method's
    # number, kept so the size of the quantisation correction is visible).
    nearest_is_lo = abs(t_lo - T_PEAK_AT_FOCUS) <= abs(t_hi - T_PEAK_AT_FOCUS)
    t_nearest = t_lo if nearest_is_lo else t_hi
    x_focus_nearest, w0_nearest = find_minimum(x, w_lo if nearest_is_lo else w_hi)

    return dict(x=x, w=w_interp, x_focus=x_focus, w0=w0, t=T_PEAK_AT_FOCUS,
               t_lo=t_lo, t_hi=t_hi, frac=frac,
               x_focus_nearest=x_focus_nearest, w0_nearest=w0_nearest,
               t_nearest=t_nearest)


cells = {name: analyse_cell(name) for name in ("paraxial", "lasy")}

x_theory = np.linspace(0, 2 * P.X_SPOT, 400)
w_theory = P.W0 * np.sqrt(1.0 + ((x_theory - P.X_SPOT) / P.X_R) ** 2)

fig, ax = plt.subplots(figsize=(9, 5.5))
ax.plot(x_theory * 1e6, w_theory * 1e6, "k-", lw=1.2, label="paraxial theory")
for (name, c), style in zip(cells.items(), ["C0o-", "C1s--"]):
    ax.plot(c["x"] * 1e6, c["w"] * 1e6, style, lw=1.0, ms=1.5, label=name)
ax.axvline(P.X_SPOT * 1e6, color="grey", ls=":", label=f"theory focus ({P.X_SPOT*1e6:.3f} um)")
ax.set(xlabel="x (um)", ylabel="beam radius w (um)",
      title="Test B (2D): LASY vs paraxial, f/1 (NA=0.5)")
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

paraxial = cells["paraxial"]
lasy = cells["lasy"]
valid = np.isfinite(paraxial["w"]) & np.isfinite(lasy["w"])
paraxial_vs_lasy_rms = float(np.sqrt(np.nanmean(
    ((paraxial["w"][valid] - lasy["w"][valid]) / lasy["w"][valid]) ** 2)))

# The RMS number alone hides WHERE along x the two models diverge -- this
# plot shows the relative difference directly as a function of x, which is
# the actual finding this test is after (unlike Campaign A, divergence from
# theory/each-other here is the EXPECTED result, not a bug signature).
fig, ax = plt.subplots(figsize=(9, 4.5))
x_um = paraxial["x"][valid] * 1e6
ax.plot(x_um, (paraxial["w"][valid] - lasy["w"][valid]) / lasy["w"][valid] * 100,
       "C3-", lw=1.0, label="paraxial vs LASY")
ax.axhline(0, color="grey", ls=":", lw=1)
ax.set(xlabel="x (um)", ylabel="relative w(x) difference, paraxial vs LASY (%)",
      title="Test B (2D): paraxial-approximation breakdown vs x")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIGURES / "paraxial_vs_lasy_diff_vs_x.png", dpi=140)
plt.close(fig)

any_cell = next(iter(cells.values()))
summary = f"""Test B (2D) -- LASY vs paraxial closed-form, f/1 tight-focusing demo
======================================================================
Beam: lambda0={P.LAMBDA0*1e6:.3f}um f/{P.F_NUMBER:.1f} NA={P.NA:.4f}
(1D-transverse slab, HALF Gouy phase -- see physics_params.py).
Paraxial theory: x_focus={P.X_SPOT*1e6:.4f}um  w0={P.W0*1e6:.4f}um

Snapshot timing: ideal t_peak@focus = {T_PEAK_AT_FOCUS*1e15:.3f} fs, but
output only exists every dt_snapshot -- bracketing snapshots at
{any_cell['t_lo']*1e15:.1f} fs and {any_cell['t_hi']*1e15:.1f} fs
(interpolation fraction {any_cell['frac']:.3f}). x_focus/w0 below use the
fitted w(x) curve LINEARLY INTERPOLATED between those two snapshots to
the exact ideal time (not the raw field -- see analyse_cell's comment).
Old method (snap to nearest snapshot, t={any_cell['t_nearest']*1e15:.1f} fs)
numbers included alongside for comparison.

x_focus / w0 per cell:
"""
for name, c in cells.items():
    err_x = abs(c["x_focus"] - P.X_SPOT) * 1e6
    err_w = abs(c["w0"] - P.W0) / P.W0 * 100
    err_x_nearest = abs(c["x_focus_nearest"] - P.X_SPOT) * 1e6
    err_w_nearest = abs(c["w0_nearest"] - P.W0) / P.W0 * 100
    summary += (f"  {name:10s}: x_focus={c['x_focus']*1e6:8.4f}um (vs theory {err_x:+.4f}um)  "
               f"w0={c['w0']*1e6:7.4f}um (vs theory {err_w:+.2f}%)\n"
               f"  {'':10s}  nearest-snapshot was: x_focus={c['x_focus_nearest']*1e6:8.4f}um "
               f"(vs theory {err_x_nearest:+.4f}um)  w0={c['w0_nearest']*1e6:7.4f}um "
               f"(vs theory {err_w_nearest:+.2f}%)\n")
summary += f"""
w(x) RMS relative difference, paraxial vs LASY (over the interior valid
range) = {paraxial_vs_lasy_rms*100:.3f}%
This is the key Test B number: how much the closed-form paraxial formula
diverges from LASY's genuine non-paraxial propagation at NA=0.5 -- a large
value here indicates the paraxial approximation is breaking down at this
f-number, as expected going in.
"""
(RESULTS / "summary.txt").write_text(summary)
print(summary)
