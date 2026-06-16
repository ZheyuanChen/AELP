"""
generate_profile.py
-------------------
Generates spatial_temporal_profile.dat for EPOCH2D arbitrary laser injection,
replicating the Gaussian beam defined in input.deck.

The .dat file encodes the normalised field envelope E(t, y) at the x_min boundary.
The phase (wavefront curvature + Gouy phase) is handled separately in the laser block.

Profile convention (must match EPOCH reader):
  Line 1:     n_t  n_y
  Line 2:     y_0  y_1  ...  y_{n_y-1}       [metres]
  Line 3:     t_0  t_1  ...  t_{n_t-1}       [seconds]
  Lines 4+:   E_matrix[i_t, i_y], one row per time step, normalised to peak = 1.0

The amplitude scaling is handled by intensity_w_cm2 in the laser block,
so this file contains only the normalised envelope [0, 1].
"""

import numpy as np

# =============================================================================
# 1. Parameters — must match input.deck constants exactly
# =============================================================================

# Laser
lambda0     = 1.0e-6          # m, laser wavelength
pulse_FWHM  = 2.0e-6          # m, intensity FWHM (spatial, used to define w_0)
pulse_tau   = 15.0e-15        # s, temporal sigma (NOT FWHM)

# w_0: beam waist at focus
# In EPOCH, gauss(y, 0, w_bnd) drives the FIELD amplitude ~ exp(-y^2 / w_bnd^2)
# For field Gaussian, intensity ~ exp(-2y^2/w_0^2), so intensity FWHM relates
# to field 1/e radius w_0 as:  FWHM_intensity = w_0 * sqrt(2*ln2)
# The deck uses w_0 = pulse_FWHM / 1.665, where 1.665 ~ 2*sqrt(ln2)
# This interprets pulse_FWHM as the *field* FWHM → w_0 = field_FWHM / (2*sqrt(ln2))
# Keep consistent with deck:
w_0 = pulse_FWHM / 1.665      # m, beam waist at focus (field 1/e radius)

# Geometry
x_min         = -40.0e-6      # m
pulse_defocus = 10.0e-6       # m, distance from x_min to focus
x_focus       = x_min + pulse_defocus   # m, focus position in sim coords = -30 um
x_spot        = x_focus - x_min        # m, distance from boundary to focus = 10 um

# Derived beam parameters at x_min boundary
x_R   = np.pi * w_0**2 / lambda0                       # Rayleigh range
w_bnd = w_0 * np.sqrt(1.0 + (x_spot / x_R)**2)        # beam radius at boundary
RC    = x_spot * (1.0 + (x_R / x_spot)**2)            # radius of curvature at boundary
gouy  = np.arctan(x_spot / x_R)                       # Gouy phase at boundary

# Simulation domain (y only — this is the injection boundary axis)
y_min = -25.0e-6   # m
y_max =  25.0e-6   # m

# Temporal window: laser is active from t_start=0 to t_end=6*pulse_tau in deck
# Centre the pulse at 3*pulse_tau (matches t_profile = gauss(time, 3*pulse_tau, pulse_tau))
t_centre = 3.0 * pulse_tau
t_start_prof = 0.0
t_end_prof   = 6.0 * pulse_tau

# =============================================================================
# 2. Grid resolution
# =============================================================================
# Spatial: match EPOCH grid spacing at x_min boundary
# EPOCH has ny=1000 cells over 50 um → dy = 0.05 um. Use the same or finer.
n_y = 1000
# Temporal: resolve the pulse well. Nyquist on the optical cycle isn't needed here
# since the carrier is handled by EPOCH; we only need to resolve the envelope.
# Use enough points to represent the Gaussian envelope accurately.
# Rule of thumb: >= 20 points per pulse_tau.
n_t = 600

y_arr = np.linspace(y_min, y_max, n_y)
t_arr = np.linspace(t_start_prof, t_end_prof, n_t)

print(f"Grid: n_t={n_t}, n_y={n_y}")
print(f"dy = {(y_max-y_min)/(n_y-1)*1e6:.4f} um,  dt = {(t_end_prof-t_start_prof)/(n_t-1)*1e15:.4f} fs")
print()
print(f"Beam parameters:")
print(f"  w_0    = {w_0*1e6:.4f} um")
print(f"  x_spot = {x_spot*1e6:.4f} um")
print(f"  x_R    = {x_R*1e6:.4f} um")
print(f"  w_bnd  = {w_bnd*1e6:.4f} um  (beam radius at injection boundary)")
print(f"  RC     = {RC*1e6:.4f} um")
print(f"  gouy   = {gouy:.6f} rad")
print()
print(f"Temporal parameters:")
print(f"  t_centre = {t_centre*1e15:.2f} fs")
print(f"  pulse_tau = {pulse_tau*1e15:.2f} fs (sigma)")
print(f"  window: [{t_start_prof*1e15:.2f}, {t_end_prof*1e15:.2f}] fs")

# =============================================================================
# 3. Build the normalised field envelope E(t, y)
# =============================================================================
# Spatial envelope: field Gaussian at the boundary
# gauss(y, 0, w_bnd) in EPOCH = exp(-y^2 / w_bnd^2)
#
# Temporal envelope: gauss(time, 3*pulse_tau, pulse_tau)
# EPOCH's gauss(x, mean, sigma) = exp(-(x-mean)^2 / sigma^2)
# NOTE: this is NOT the standard statistics convention exp(-x^2/2sigma^2)
# It is exp(-x^2/sigma^2), i.e. the sigma here is the 1/e half-width of the field.

T, Y = np.meshgrid(t_arr, y_arr, indexing='ij')  # shape (n_t, n_y)

env_t = np.exp(-((T - t_centre) / pulse_tau)**2)   # temporal envelope
env_y = np.exp(-(Y / w_bnd)**2)                    # spatial envelope at boundary

E_matrix = env_t * env_y   # normalised field envelope, peak = 1.0

# Sanity check
print(f"\nE_matrix shape: {E_matrix.shape}")
print(f"E_matrix max:   {E_matrix.max():.6f}  (should be 1.0)")
print(f"E_matrix at t_centre, y=0: {E_matrix[n_t//2, n_y//2]:.6f}")

# =============================================================================
# 4. Write to file
# =============================================================================
filename = "spatial_temporal_profile.dat"

with open(filename, "w") as f:
    # Line 1: grid dimensions
    f.write(f"{n_t} {n_y}\n")

    # Line 2: y-axis values [m], space-delimited
    np.savetxt(f, y_arr.reshape(1, -1), fmt="%.8e", delimiter=" ")

    # Line 3: t-axis values [s], space-delimited
    np.savetxt(f, t_arr.reshape(1, -1), fmt="%.8e", delimiter=" ")

    # Lines 4+: E_matrix, one row per time step
    # Row i_t contains E[i_t, 0], E[i_t, 1], ..., E[i_t, n_y-1]
    np.savetxt(f, E_matrix, fmt="%.8e", delimiter=" ")

print(f"\nWrote '{filename}'  ({n_t} x {n_y} = {n_t*n_y} values)")

# =============================================================================
# 5. Quick verification plot (optional — comment out if running headless)
# =============================================================================
try:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    # Panel 1: full 2D envelope
    im = axes[0].pcolormesh(y_arr*1e6, t_arr*1e15, E_matrix, cmap='inferno', shading='auto')
    axes[0].set_xlabel("y (µm)")
    axes[0].set_ylabel("t (fs)")
    axes[0].set_title("E(t, y) — normalised envelope")
    plt.colorbar(im, ax=axes[0])

    # Panel 2: temporal slice at y=0
    axes[1].plot(t_arr*1e15, E_matrix[:, n_y//2])
    axes[1].set_xlabel("t (fs)")
    axes[1].set_ylabel("E (normalised)")
    axes[1].set_title(f"Temporal envelope at y=0")
    axes[1].axvline(t_centre*1e15, color='r', linestyle='--', label=f't_centre={t_centre*1e15:.1f} fs')
    axes[1].legend()

    # Panel 3: spatial slice at t=t_centre
    i_t_peak = np.argmin(np.abs(t_arr - t_centre))
    axes[2].plot(y_arr*1e6, E_matrix[i_t_peak, :])
    axes[2].set_xlabel("y (µm)")
    axes[2].set_ylabel("E (normalised)")
    axes[2].set_title(f"Spatial envelope at t=t_centre\nw_bnd = {w_bnd*1e6:.3f} µm")
    axes[2].axvline( w_bnd*1e6, color='r', linestyle='--', label=f'±w_bnd')
    axes[2].axvline(-w_bnd*1e6, color='r', linestyle='--')
    axes[2].legend()

    plt.tight_layout()
    plt.savefig("profile_check.png", dpi=150)
    print("Saved verification plot to 'profile_check.png'")
    plt.show()

except ImportError:
    print("matplotlib not available — skipping verification plot")