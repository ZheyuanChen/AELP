"""
Test 1 (2D) profile generator
=============================
Writes ``temporal_spatial_profile.dat`` for the numerical run: the normalised
field envelope E(t, y) at the x_min boundary for the Gaussian beam defined in
input.deck. Amplitude only — the phase (wavefront curvature + Gouy) is supplied
analytically in the laser block, identical to the analytical run.

File format (consumed by epoch2d custom_laser):
    line 1   : n_t  n_y
    line 2   : y_0 ... y_{n_y-1}     [metres]
    line 3   : t_0 ... t_{n_t-1}     [seconds]
    lines 4+ : n_t rows, each n_y values  (row j is time t_j), peak normalised to 1
"""

import numpy as np
from pathlib import Path

# --- Parameters: must match input.deck ---
lambda0    = 1.0e-6
pulse_FWHM = 2.0e-6
pulse_tau  = 15.0e-15
w_0        = pulse_FWHM / 1.665

pulse_defocus = 10.0e-6
x_spot        = pulse_defocus
x_R   = np.pi * w_0**2 / lambda0
w_bnd = w_0 * np.sqrt(1.0 + (x_spot / x_R)**2)   # spot size at boundary

# Transverse domain (matches input.deck y extent) and resolution.
# n_y must be >= the deck's ny so the boundary profile is fully resolved.
y_min, y_max = -12.5e-6, 12.5e-6
n_y = 1250

# Temporal window: laser active from 0 to 6*pulse_tau, centred at 3*pulse_tau.
t_start, t_end = 0.0, 6.0 * pulse_tau
t_centre = 3.0 * pulse_tau
n_t = 800

y_arr = np.linspace(y_min, y_max, n_y)
t_arr = np.linspace(t_start, t_end, n_t)

# Separable envelope: spatial Gaussian at the boundary x temporal Gaussian.
# EPOCH's gauss(u, mu, sig) = exp(-((u-mu)/sig)^2).
T, Y = np.meshgrid(t_arr, y_arr, indexing="ij")          # shape (n_t, n_y)
env = np.exp(-((T - t_centre) / pulse_tau)**2) * np.exp(-(Y / w_bnd)**2)
env /= env.max()

out = Path(__file__).parent / "temporal_spatial_profile.dat"
with open(out, "w") as f:
    f.write(f"{n_t} {n_y}\n")
    np.savetxt(f, y_arr.reshape(1, -1), fmt="%.8e", delimiter=" ")
    np.savetxt(f, t_arr.reshape(1, -1), fmt="%.8e", delimiter=" ")
    np.savetxt(f, env, fmt="%.8e", delimiter=" ")

print(f"Wrote {out}")
print(f"  grid: n_t={n_t}, n_y={n_y};  w_bnd={w_bnd*1e6:.3f} um;  peak={env.max():.4f}")
