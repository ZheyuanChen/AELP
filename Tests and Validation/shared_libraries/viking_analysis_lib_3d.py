"""
viking_analysis_lib_3d.py
==========================
Shared, head-less analysis helpers for the 3D custom-laser-injection
verification tests (Viking_results/3D). Companion to the 2D package's
``common/viking_analysis_lib.py`` -- kept as a separate module because the
3D tests need genuinely 3D-specific analysis (2D Gaussian beam-waist
fitting per x-plane, transverse power integrals) that the 2D lib has no
use for.

Design notes
------------
* Field data is loaded with ``sdf_xarray.open_mfdataset(..., separate_times=True)``
  (project convention, CLAUDE.md) rather than the ``xr.open_mfdataset`` +
  ``SDFPreprocess`` pattern used for moving-window runs -- these are static
  boxes, so the simpler/more memory-efficient path applies. Pass
  ``data_vars=[...]`` to avoid loading fields you don't need (these are 3D
  arrays and can be large at full Viking resolution).
* No animations (head-less, matches the 2D package). Static PNGs into
  ``results/figures/``.
* Coordinates kept in SI internally.
* Matplotlib uses the ``Agg`` backend.

Required packages: numpy, scipy, matplotlib, xarray, sdf-xarray.
"""

import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sdf_xarray as sdfxr
from scipy.signal import hilbert
from scipy.optimize import curve_fit

EPSILON0 = 8.8541878128e-12
C_LIGHT = 299792458.0
MU0 = 1.0 / (EPSILON0 * C_LIGHT ** 2)


# ----------------------------------------------------------------------
# I/O helpers (mirrors the 2D lib's conventions)
# ----------------------------------------------------------------------
def find_sdf_dir(run_dir):
    if glob.glob(os.path.join(run_dir, "*.sdf")):
        return run_dir
    sub = os.path.join(run_dir, "sdf_files")
    if glob.glob(os.path.join(sub, "*.sdf")):
        return sub
    raise FileNotFoundError(f"No .sdf files found in {run_dir} or {run_dir}/sdf_files")


def load_fields(run_dir, data_vars):
    """
    Load the requested field variables from a run directory as an xarray
    Dataset. Uses the project-standard separate_times=True loader. The time
    coordinate is named 'time0' by this loader (vs 'time' for the
    xr.open_mfdataset+SDFPreprocess path used elsewhere in this project) --
    callers should use ds['time0'] here.
    """
    sdf_dir = find_sdf_dir(run_dir)
    files = sorted(glob.glob(os.path.join(sdf_dir, "*.sdf")))
    if not files:
        raise FileNotFoundError(f"No .sdf files in {sdf_dir}")
    return sdfxr.open_mfdataset(files, separate_times=True, data_vars=data_vars)


def make_results_dirs(base_dir):
    results = os.path.join(base_dir, "results")
    figures = os.path.join(results, "figures")
    os.makedirs(figures, exist_ok=True)
    return results, figures


def write_text(path, text):
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def write_metrics_csv(path, rows, header):
    lines = [",".join(header)]
    for row in rows:
        lines.append(",".join(str(x) for x in row))
    write_text(path, "\n".join(lines) + "\n")


# ----------------------------------------------------------------------
# Envelope extraction
# ----------------------------------------------------------------------
def hilbert_envelope_along_x(field_xyz):
    """
    Amplitude envelope of an oscillating field along axis 0 (propagation
    direction, X_Grid_mid). ``field_xyz`` is a plain numpy array (x,y,z).
    """
    return np.abs(hilbert(np.asarray(field_xyz), axis=0))


# ----------------------------------------------------------------------
# 2D Gaussian beam-waist fitting
# ----------------------------------------------------------------------
def _gauss2d(coords, amp, y0, z0, w):
    Y, Z = coords
    return amp * np.exp(-((Y - y0) ** 2 + (Z - z0) ** 2) / w ** 2)


def fit_transverse_gaussian(envelope_slab, y, z, p0=None):
    """
    Fit envelope_slab(y,z) [the field ENVELOPE, not intensity] to a
    circular Gaussian amp*exp(-((y-y0)^2+(z-z0)^2)/w^2) via least squares.

    Returns dict: amp, y0, z0, w, w_err (fit std error on w), success (bool).
    Falls back to an intensity-weighted second-moment estimate (with
    w_err = nan) if the fit fails to converge -- this keeps a waist scan
    running through low-SNR planes (e.g. near t=0 before the pulse
    arrives) instead of raising.
    """
    Y, Z = np.meshgrid(y, z, indexing="ij")
    data = np.asarray(envelope_slab)
    peak = float(data.max())

    if p0 is None:
        # second-moment initial guess
        I = data ** 2
        tot = I.sum()
        if tot <= 0 or peak <= 0:
            return dict(amp=0.0, y0=np.nan, z0=np.nan, w=np.nan,
                        w_err=np.nan, success=False)
        y0_guess = float((Y * I).sum() / tot)
        z0_guess = float((Z * I).sum() / tot)
        var_guess = float((((Y - y0_guess) ** 2 + (Z - z0_guess) ** 2) * I).sum()
                          / tot)
        w_guess = np.sqrt(max(var_guess, (y[1] - y[0]) ** 2))
        p0 = (peak, y0_guess, z0_guess, w_guess)

    try:
        popt, pcov = curve_fit(
            _gauss2d, (Y.ravel(), Z.ravel()), data.ravel(), p0=p0,
            maxfev=5000)
        perr = np.sqrt(np.diag(pcov))
        return dict(amp=float(popt[0]), y0=float(popt[1]), z0=float(popt[2]),
                    w=float(abs(popt[3])), w_err=float(perr[3]), success=True)
    except Exception:
        return dict(amp=float(p0[0]), y0=float(p0[1]), z0=float(p0[2]),
                    w=float(p0[3]), w_err=np.nan, success=False)


def waist_scan(envelope_xyz, x, y, z, snr_threshold_frac=0.05):
    """
    Fit a circular Gaussian at every x-plane of a (x,y,z) envelope array.
    Planes whose peak is below snr_threshold_frac of the global peak are
    skipped (returned as nan) -- avoids fitting noise where the pulse
    envelope hasn't arrived yet / has already passed.

    Returns a dict of numpy arrays (length len(x)): w, w_err, y0, z0, amp,
    success (bool array).
    """
    global_peak = float(np.abs(envelope_xyz).max())
    n = len(x)
    w = np.full(n, np.nan)
    w_err = np.full(n, np.nan)
    y0 = np.full(n, np.nan)
    z0 = np.full(n, np.nan)
    amp = np.zeros(n)
    success = np.zeros(n, dtype=bool)

    for ix in range(n):
        slab = envelope_xyz[ix]
        if float(slab.max()) < snr_threshold_frac * global_peak:
            continue
        fit = fit_transverse_gaussian(slab, y, z)
        w[ix] = fit["w"]
        w_err[ix] = fit["w_err"]
        y0[ix] = fit["y0"]
        z0[ix] = fit["z0"]
        amp[ix] = fit["amp"]
        success[ix] = fit["success"]

    return dict(w=w, w_err=w_err, y0=y0, z0=z0, amp=amp, success=success)


def find_waist_minimum(x, w):
    """Parabolic interpolation around the minimum of w(x) for a sub-cell
    estimate of x_focus and w0. Returns (x_focus, w0) or (nan, nan) if
    too few valid points."""
    valid = np.isfinite(w)
    if valid.sum() < 3:
        return np.nan, np.nan
    xi = x[valid]
    wi = w[valid]
    i = int(np.argmin(wi))
    if i == 0 or i == len(wi) - 1:
        return float(xi[i]), float(wi[i])
    x0, x1, x2 = xi[i - 1], xi[i], xi[i + 1]
    y0, y1, y2 = wi[i - 1], wi[i], wi[i + 1]
    denom = (x0 - x1) * (x0 - x2) * (x2 - x1)
    if denom == 0:
        return float(x1), float(y1)
    a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom
    b = (x2 ** 2 * (y0 - y1) + x1 ** 2 * (y2 - y0) + x0 ** 2 * (y1 - y2)) / denom
    if a == 0:
        return float(x1), float(y1)
    x_focus = -b / (2 * a)
    # evaluate the parabola at x_focus for the fitted w0
    c = y0 - a * x0 ** 2 - b * x0
    w0 = a * x_focus ** 2 + b * x_focus + c
    return float(x_focus), float(w0)


# ----------------------------------------------------------------------
# Energy / power conservation
# ----------------------------------------------------------------------
def transverse_intensity_power(envelope_xyz, y, z):
    """
    P(x) = integral |envelope|^2 dy dz at each x-plane -- proportional to
    the transverse-plane optical power (up to the constant c*epsilon0/2).
    Returns a 1D array of length len(x). Should be flat vs x for an
    energy-conserving injected wavefront (aside from small numerical
    dispersion / absorbing-boundary losses).
    """
    dy = y[1] - y[0]
    dz = z[1] - z[0]
    return np.sum(np.asarray(envelope_xyz) ** 2, axis=(1, 2)) * dy * dz


def expected_envelope_shape(x, t_snapshot, t_centre, tau, c=C_LIGHT):
    """
    Predicted (unnormalised) envelope^2 shape vs x at a fixed snapshot
    time, for a deck `t_profile = gauss(time, t_centre, tau)` temporal
    envelope injected at x=0: envelope(x) = gauss(t_snapshot - x/c,
    t_centre, tau), i.e. the retarded-time value of the boundary's own
    temporal profile. A genuinely pulsed (non-CW) beam does NOT have flat
    transverse power P(x) at a single snapshot -- different x-planes
    sample different points of the temporal envelope, delayed by the
    light travel time x/c. Dividing measured P(x) by this predicted shape
    (see energy_conservation_ratio) isolates real conservation violations
    from this expected, benign pulse-envelope effect.
    """
    return np.exp(-((t_snapshot - np.asarray(x) / c - t_centre) / tau) ** 2) ** 2


def energy_conservation_ratio(P_measured, x, t_snapshot, t_centre, tau,
                              c=C_LIGHT):
    """
    Ratio of measured transverse power P(x) to the expected retarded-time
    envelope-squared shape, normalised to 1 at the shape's own peak. Flat
    at ~1.0 confirms energy conservation (self-consistency of the injected
    wavefront); systematic drift indicates spurious loss/gain. Expect
    degradation right at the injection boundary (near-field transient) and
    at the far absorbing boundary (imperfect absorption) -- mask those
    when interpreting.
    """
    shape = expected_envelope_shape(x, t_snapshot, t_centre, tau, c)
    predicted = shape / shape.max() * np.nanmax(P_measured)
    return P_measured / predicted


def poynting_x_power(ey_xyz, ez_xyz, by_xyz, bz_xyz, y, z):
    """
    Time-averaged-free instantaneous x-directed Poynting flux, integrated
    over each transverse plane: P(x) = integral (Ey*Bz - Ez*By)/mu0 dy dz.
    A more physically complete energy-transport check than the envelope
    power alone (uses both E and B, so it is not fooled by a
    standing-wave / partially-reflected component near the boundary).
    """
    Sx = (np.asarray(ey_xyz) * np.asarray(bz_xyz)
          - np.asarray(ez_xyz) * np.asarray(by_xyz)) / MU0
    dy = y[1] - y[0]
    dz = z[1] - z[0]
    return np.sum(Sx, axis=(1, 2)) * dy * dz


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------
def save_waist_comparison(x_um, curves, x_focus_theory_um, w0_theory_um,
                          out_path, title="Beam waist vs x"):
    """
    curves: list of (label, x_um, w_um, style) tuples (x_um may differ in
    length per curve, e.g. numerical vs analytical vs theory).
    """
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for label, xc, wc, style in curves:
        ax.plot(xc, wc, style, lw=1.0, ms=1.5, label=label)
    ax.axvline(x_focus_theory_um, color="grey", ls=":", lw=1,
               label=f"theory focus ({x_focus_theory_um:.3f} um)")
    ax.axhline(w0_theory_um, color="grey", ls=":", lw=1)
    ax.set(xlabel="x (um)", ylabel="beam radius w (um)", title=title)
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def save_power_conservation(x_um, curves, out_path,
                            title="Transverse power vs x (energy conservation)"):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for label, xc, pc, style in curves:
        pc_norm = np.asarray(pc) / np.nanmax(np.abs(pc))
        ax.plot(xc, pc_norm, style, lw=1.6, label=label)
    ax.axhline(1.0, color="grey", ls=":", lw=1)
    ax.set(xlabel="x (um)", ylabel="P(x) / max(P)", title=title)
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def save_field_slice(field_yz, y_um, z_um, title, out_path, cmap="RdBu_r"):
    fig, ax = plt.subplots(figsize=(6, 5.2))
    vmax = float(np.abs(field_yz).max())
    vmax = vmax if vmax > 0 else 1.0
    im = ax.imshow(field_yz.T, origin="lower", aspect="auto",
                   extent=[y_um.min(), y_um.max(), z_um.min(), z_um.max()],
                   cmap=cmap, vmin=-vmax, vmax=vmax)
    ax.set(title=title, xlabel="y (um)", ylabel="z (um)")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)