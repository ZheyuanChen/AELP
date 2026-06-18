"""
generate_profile.py
-------------------
Generates temporal_spatial_profile.dat for EPOCH2D arbitrary laser injection,
using a uniform (flat) spatial profile: E(t, y) = gauss(t) * 1.

This is the simplest possible test case for comparing the .dat file injection
against the analytical profile = 1 in EPOCH.

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

pulse_tau = 15.0e-15        # s, temporal sigma (NOT FWHM)

# Simulation domain (y only — this is the injection boundary axis)
y_min = -12.5e-6   # m
y_max =  12.5e-6   # m

# Temporal window: laser is active from t_start=0 to t_end=6*pulse_tau in deck
# Centre the pulse at 3*pulse_tau (matches t_profile = gauss(time, 3*pulse_tau, pulse_tau))
t_centre = 3.0 * pulse_tau
t_start_prof = 0.0
t_end_prof   = 6.0 * pulse_tau

# =============================================================================
# 2. Grid resolution
# =============================================================================
# Spatial: match EPOCH grid spacing at x_min boundary
# EPOCH has ny=500 cells over 25 um → dy = 0.05 um. Use the same or finer.
n_y = 3000
# Temporal: use enough points to resolve the Gaussian envelope accurately.
n_t = 1900

y_arr = np.linspace(y_min, y_max, n_y)
t_arr = np.linspace(t_start_prof, t_end_prof, n_t)

print(f"Grid: n_t={n_t}, n_y={n_y}")
print(f"dy = {(y_max-y_min)/(n_y-1)*1e6:.4f} um,  dt = {(t_end_prof-t_start_prof)/(n_t-1)*1e15:.4f} fs")
print()
print(f"Temporal parameters:")
print(f"  t_centre = {t_centre*1e15:.2f} fs")
print(f"  pulse_tau = {pulse_tau*1e15:.2f} fs (sigma)")
print(f"  window: [{t_start_prof*1e15:.2f}, {t_end_prof*1e15:.2f}] fs")

# =============================================================================
# 3. Build the normalised field envelope E(t, y)
# =============================================================================
# Spatial envelope: uniform (= 1 everywhere)
#
# Temporal envelope: gauss(time, 3*pulse_tau, pulse_tau)
# EPOCH's gauss(x, mean, sigma) = exp(-(x-mean)^2 / sigma^2)
# NOTE: this is NOT the standard statistics convention exp(-x^2/2sigma^2)
# It is exp(-x^2/sigma^2), i.e. the sigma here is the 1/e half-width of the field.

T, Y = np.meshgrid(t_arr, y_arr, indexing='ij')  # shape (n_t, n_y)

env_t = np.exp(-((T - t_centre) / pulse_tau)**2)   # temporal envelope

# Uniform spatial profile: the entire matrix is just the temporal envelope
E_matrix = env_t   # normalised field envelope, peak = 1.0

# Sanity check
print(f"\nE_matrix shape: {E_matrix.shape}")
print(f"E_matrix max:   {E_matrix.max():.6f}  (should be 1.0)")
print(f"E_matrix at t_centre, y=0: {E_matrix[n_t//2, n_y//2]:.6f}")
print(f"E_matrix at t_centre, y=y_min: {E_matrix[n_t//2, 0]:.6f}  (should equal y=0 for uniform)")

# =============================================================================
# 4. Write to file
# =============================================================================
filename = "temporal_spatial_profile.dat"

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

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Panel 1: full 2D envelope
    im = axes[0].pcolormesh(y_arr*1e6, t_arr*1e15, E_matrix, cmap='inferno', shading='auto')
    axes[0].set_xlabel("y (µm)")
    axes[0].set_ylabel("t (fs)")
    axes[0].set_title("E(t, y) — normalised envelope (uniform)")
    plt.colorbar(im, ax=axes[0])

    # Panel 2: temporal slice at y=0
    axes[1].plot(t_arr*1e15, E_matrix[:, n_y//2])
    axes[1].set_xlabel("t (fs)")
    axes[1].set_ylabel("E (normalised)")
    axes[1].set_title("Temporal envelope at y=0")
    axes[1].axvline(t_centre*1e15, color='r', linestyle='--', label=f't_centre={t_centre*1e15:.1f} fs')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig("profile_check.png", dpi=150)
    print("Saved verification plot to 'profile_check.png'")
    plt.show()

except ImportError:
    print("matplotlib not available — skipping verification plot")
