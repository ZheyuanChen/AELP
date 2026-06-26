"""
Test 2 generator — hand-written ANALYTICAL amplitude (.dat)
==========================================================
Writes ``amplitude_numerical.dat`` (the analytical Gaussian envelope E(t,y),
amplitude only) into the run directories that consume it:
    numerical/                 (reference run)
    numerical_amp_lasy_phase/  (isolate-phase cell)

This MUST share the identical (y, t) grid with the LASY profiles
(generate_lasy_profiles.py), because the isolate-phase cell pairs this
amplitude with the LASY phase on the same grid.
"""

import numpy as np
from pathlib import Path

# --- Parameters: must match the test2 decks (paraxial: pulse_FWHM = 10 um) ---
lambda0    = 1.0e-6
pulse_FWHM = 10.0e-6
pulse_tau  = 15.0e-15
w_0        = pulse_FWHM / 1.665
x_spot     = 10.0e-6
x_R        = np.pi * w_0**2 / lambda0
w_bnd      = w_0 * np.sqrt(1.0 + (x_spot / x_R)**2)

# Shared grid (KEEP IN SYNC with generate_lasy_profiles.py)
Y_MIN, Y_MAX, N_Y = -25.0e-6, 25.0e-6, 2000
T_START, T_END, N_T = 0.0, 6.0 * pulse_tau, 1000
t_centre = 3.0 * pulse_tau

y_arr = np.linspace(Y_MIN, Y_MAX, N_Y)
t_arr = np.linspace(T_START, T_END, N_T)

T, Y = np.meshgrid(t_arr, y_arr, indexing="ij")          # (n_t, n_y)
env = np.exp(-((T - t_centre) / pulse_tau)**2) * np.exp(-(Y / w_bnd)**2)
env /= env.max()

HERE = Path(__file__).parent
DESTS = [HERE / "numerical", HERE / "numerical_amp_lasy_phase"]
for dest in DESTS:
    dest.mkdir(exist_ok=True)
    out = dest / "amplitude_numerical.dat"
    with open(out, "w") as f:
        f.write(f"{N_T} {N_Y}\n")
        np.savetxt(f, y_arr.reshape(1, -1), fmt="%.8e", delimiter=" ")
        np.savetxt(f, t_arr.reshape(1, -1), fmt="%.8e", delimiter=" ")
        np.savetxt(f, env, fmt="%.8e", delimiter=" ")
    print(f"Wrote {out}")

print(f"  grid: n_t={N_T}, n_y={N_Y};  w_bnd={w_bnd*1e6:.3f} um;  peak={env.max():.4f}")
