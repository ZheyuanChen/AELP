"""
viking_analysis_lib.py
======================
Shared, head-less analysis helpers for the Viking laser-injection comparison
tests. These replace the interactive Jupyter notebooks used on the laptop:
instead of inline plots and animations, every routine writes PNG figures and
plain-text/CSV results to a structured ``results/`` folder, so the whole
analysis can run unattended on a head-less HPC node.

Design notes
------------
* Field data (E_y) is read with ``sdf_xarray`` (the project convention for
  field/grid data). Particle data is **not** handled here.
* No animations (head-less). Where the notebooks animated over time or z, we
  instead emit a montage of static panels at selected indices.
* Coordinates are kept in SI internally; we only multiply by 1e6 / 1e15 when
  setting plot axes, so we never mutate the dataset (and avoid depending on the
  ``.epoch`` accessor, keeping the scripts portable to a bare Viking env).
* Matplotlib uses the ``Agg`` backend (no display needed).

Required packages: numpy, scipy, matplotlib, xarray, sdf-xarray.
"""

import glob
import os

import matplotlib
matplotlib.use("Agg")  # head-less backend — must be set before pyplot import
import matplotlib.pyplot as plt
import numpy as np
import sdf_xarray as sdfxr
from scipy.signal import hilbert

# ----------------------------------------------------------------------
# Physical constants and a0 normalisation
# ----------------------------------------------------------------------
M_E = 9.1093837015e-31      # electron mass [kg]
C_LIGHT = 299792458.0       # speed of light [m/s]
E_CHARGE = 1.602176634e-19  # elementary charge [C]


def a0_norm(lambda0):
    """Return the E-field [V/m] that corresponds to a0 = 1 for wavelength lambda0."""
    omega0 = 2.0 * np.pi * C_LIGHT / lambda0
    return M_E * C_LIGHT * omega0 / E_CHARGE


# ----------------------------------------------------------------------
# I/O helpers
# ----------------------------------------------------------------------
def find_sdf_dir(run_dir):
    """
    Locate the directory that actually contains the .sdf files for a run.

    EPOCH may be run with the deck in ``run_dir`` (output alongside it) or with
    a dedicated ``sdf_files`` subdirectory. Check both.
    """
    if glob.glob(os.path.join(run_dir, "*.sdf")):
        return run_dir
    sub = os.path.join(run_dir, "sdf_files")
    if glob.glob(os.path.join(sub, "*.sdf")):
        return sub
    raise FileNotFoundError(f"No .sdf files found in {run_dir} or {run_dir}/sdf_files")


def load_ey(run_dir, lambda0, data_var="Electric_Field_Ey"):
    """
    Load E_y from a run directory as an xarray DataArray normalised to a0.

    The returned array is lazy (dask-backed) where possible; call ``.load()`` or
    slice-then-``.load()`` to materialise only what is needed. Coordinates are
    left in SI units.
    """
    sdf_dir = find_sdf_dir(run_dir)
    files = sorted(glob.glob(os.path.join(sdf_dir, "*.sdf")))
    ds = sdfxr.open_mfdataset(files, data_vars=[data_var])
    ey = ds[data_var] / a0_norm(lambda0)
    ey.attrs["units"] = "$a_0$"
    ey.attrs["long_name"] = "$E_y$"
    return ey


# ----------------------------------------------------------------------
# Analysis primitives
# ----------------------------------------------------------------------
def hilbert_envelope_along_x(da):
    """
    Amplitude envelope of an oscillating field along the propagation axis
    (X_Grid_mid), returned as a plain numpy array the same shape as ``da``.
    """
    x_idx = list(da.dims).index("X_Grid_mid")
    return np.abs(hilbert(np.asarray(da.values), axis=x_idx))


def difference_metrics(diff_values, envelope, env_threshold_frac=0.01):
    """
    Compute summary metrics for a difference field.

    Parameters
    ----------
    diff_values : ndarray
        variant - reference (in a0).
    envelope : ndarray
        Reference Hilbert envelope (same shape), used as the relative-error
        denominator; cells below ``env_threshold_frac`` of the peak are masked.

    Returns
    -------
    dict with keys: max_abs, max_rel, mean_rel, peak_env.
    """
    peak_env = float(envelope.max())
    denom = np.where(envelope > env_threshold_frac * peak_env, envelope, np.nan)
    rel = np.abs(diff_values) / denom
    return {
        "max_abs": float(np.abs(diff_values).max()),
        "max_rel": float(np.nanmax(rel)),
        "mean_rel": float(np.nanmean(rel)),
        "peak_env": peak_env,
    }


def orient_xy(da):
    """Transpose a 2D (x,y) slice so X is horizontal when plotted with imshow-like APIs."""
    if "X_Grid_mid" in da.dims and "Y_Grid_mid" in da.dims:
        other = [d for d in da.dims if d not in ("X_Grid_mid", "Y_Grid_mid")]
        return da.transpose(*other, "Y_Grid_mid", "X_Grid_mid")
    return da


# ----------------------------------------------------------------------
# Plotting helpers (all save to file, none display)
# ----------------------------------------------------------------------
def _coord_um(da, name):
    return da.coords[name].values * 1e6


def save_triptych(ref_slice, var_slice, diff_slice, time_fs, out_path,
                  ref_label="Analytical", var_label="Numerical",
                  x_name="X_Grid_mid", y_name="Y_Grid_mid"):
    """
    Save a 3-panel figure: reference, variant, and their difference, for one 2D
    slice (e.g. a z=0 plane or a full 2D run) at a single time.
    """
    fig, ax = plt.subplots(1, 3, figsize=(18, 4.6))
    xu = _coord_um(ref_slice, x_name)
    yu = _coord_um(ref_slice, y_name)
    ext = [xu.min(), xu.max(), yu.min(), yu.max()]

    a = orient_xy(ref_slice).values
    b = orient_xy(var_slice).values
    d = orient_xy(diff_slice).values
    vmax_field = max(np.abs(a).max(), np.abs(b).max())
    vmax_diff = np.abs(d).max() if np.abs(d).max() > 0 else 1e-30

    im0 = ax[0].imshow(a, origin="lower", aspect="auto", extent=ext,
                       cmap="RdBu_r", vmin=-vmax_field, vmax=vmax_field)
    ax[0].set(title=f"{ref_label} [$a_0$]  |  t = {time_fs:.1f} fs",
              xlabel="x (µm)", ylabel="y (µm)"); fig.colorbar(im0, ax=ax[0])
    im1 = ax[1].imshow(b, origin="lower", aspect="auto", extent=ext,
                       cmap="RdBu_r", vmin=-vmax_field, vmax=vmax_field)
    ax[1].set(title=f"{var_label} [$a_0$]  |  t = {time_fs:.1f} fs",
              xlabel="x (µm)", ylabel="y (µm)"); fig.colorbar(im1, ax=ax[1])
    im2 = ax[2].imshow(d, origin="lower", aspect="auto", extent=ext,
                       cmap="RdBu_r", vmin=-vmax_diff, vmax=vmax_diff)
    ax[2].set(title=f"Difference ({var_label} - {ref_label}) [$a_0$]",
              xlabel="x (µm)", ylabel="y (µm)"); fig.colorbar(im2, ax=ax[2])

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def save_relative_diff_panel(envelope_slice, rel_slice, time_fs, out_path,
                             x_name="X_Grid_mid", y_name="Y_Grid_mid",
                             ref_label="Analytical"):
    """Save a 2-panel figure: reference envelope and relative difference."""
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    xu = _coord_um(envelope_slice, x_name)
    yu = _coord_um(envelope_slice, y_name)
    ext = [xu.min(), xu.max(), yu.min(), yu.max()]

    env = orient_xy(envelope_slice).values
    rel = orient_xy(rel_slice).values

    im0 = ax[0].imshow(env, origin="lower", aspect="auto", extent=ext, cmap="inferno")
    ax[0].set(title=f"{ref_label} envelope [$a_0$]  |  t = {time_fs:.1f} fs",
              xlabel="x (µm)", ylabel="y (µm)"); fig.colorbar(im0, ax=ax[0])
    im1 = ax[1].imshow(rel, origin="lower", aspect="auto", extent=ext, cmap="viridis")
    ax[1].set(title=f"Relative difference  |  t = {time_fs:.1f} fs",
              xlabel="x (µm)", ylabel="y (µm)"); fig.colorbar(im1, ax=ax[1])

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def save_lineout(x_um, curves, time_fs, out_path, xlabel="x (µm)",
                 ylabel="$E_y$ [$a_0$]", title=None):
    """
    Save an overlaid 1D line-out plus a residual panel.

    ``curves`` is a list of (label, values, style) tuples; the first curve is
    treated as the reference for the residual panel.
    """
    fig, ax = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                           gridspec_kw={"height_ratios": [2, 1]})
    ref_label, ref_vals, _ = curves[0]
    for label, vals, style in curves:
        ax[0].plot(x_um, vals, style, lw=1.3, label=label)
    ax[0].set(ylabel=ylabel,
              title=title or f"Line-out  |  t = {time_fs:.1f} fs")
    ax[0].legend(loc="upper right"); ax[0].grid(alpha=0.3)

    for label, vals, style in curves[1:]:
        ax[1].plot(x_um, np.asarray(vals) - np.asarray(ref_vals), lw=1.0,
                   label=f"{label} - {ref_label}")
    ax[1].axhline(0, color="grey", lw=0.5, ls=":")
    ax[1].set(xlabel=xlabel, ylabel=r"residual [$a_0$]"); ax[1].grid(alpha=0.3)
    if len(curves) > 2:
        ax[1].legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


# ----------------------------------------------------------------------
# Output-folder helpers
# ----------------------------------------------------------------------
def make_results_dirs(base_dir):
    """Create results/ and results/figures/ under base_dir; return (results, figures)."""
    results = os.path.join(base_dir, "results")
    figures = os.path.join(results, "figures")
    os.makedirs(figures, exist_ok=True)
    return results, figures


def write_text(path, text):
    """Write a UTF-8 text file."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def write_metrics_csv(path, rows, header):
    """Write a simple CSV from a list of row-tuples and a header tuple."""
    lines = [",".join(header)]
    for row in rows:
        lines.append(",".join(str(x) for x in row))
    write_text(path, "\n".join(lines) + "\n")
