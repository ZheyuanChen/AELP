"""
Test 2 generator — LASY amplitude AND phase (.dat)
=================================================
Creates the Gaussian beam in LASY, propagates it from focus back to the x_min
injection boundary, and writes two files:
    amplitude_lasy.dat  -> lasy/, lasy_amp_deck_phase/
    phase_lasy.dat      -> lasy/, numerical_amp_lasy_phase/

Field-convention chain (see daily log, 25 June):
    LASY:  E_phys = Re[E_env e^{-i w0 t}] = |E| cos(w0 t - phi),  phi = angle(E_env)
    EPOCH: E_phys = amp * profile * sin(w0 t + phase)
    => profile = |E| / peak,   phase = -phi + pi/2

CARRIER-OFFSET FIX: LASY's propagator imprints an arbitrary global (piston)
phase on the envelope (on-axis ~ -pi here). Left in, it would appear as a
constant phase offset between this run and the deck-phase reference and produce
a huge pointwise error (the "137 %" bug). We therefore reference the phase to
its on-axis value at the pulse-peak time and pin it to the deck's on-axis phase
(-gouy), which removes the piston while preserving the physical wavefront.

KEEP the (y, t) grid in sync with generate_numerical_profile.py.
"""

import numpy as np
from pathlib import Path
from scipy.interpolate import RegularGridInterpolator
from lasy.laser import Laser
from lasy.profiles.gaussian_profile import GaussianProfile

# --- Parameters: must match the test2 decks ---
lambda0    = 1.0e-6
pulse_FWHM = 10.0e-6
pulse_tau  = 15.0e-15
w_0        = pulse_FWHM / 1.665
x_spot     = 10.0e-6
x_R   = np.pi * w_0**2 / lambda0
w_bnd = w_0 * np.sqrt(1.0 + (x_spot / x_R)**2)
gouy  = np.arctan(x_spot / x_R)

# Shared grid (KEEP IN SYNC with generate_numerical_profile.py)
Y_MIN, Y_MAX, N_Y = -25.0e-6, 25.0e-6, 2000
T_START, T_END, N_T = 0.0, 6.0 * pulse_tau, 1000
t_centre = 3.0 * pulse_tau

# --- Build and propagate the LASY beam (rt geometry, mode 0) ---
profile = GaussianProfile(lambda0, (1, 0), 1.0, w_0, pulse_tau, t_centre)
laser = Laser(dim="rt", lo=(0, T_START), hi=(5 * w_bnd, T_END),
              npoints=(N_Y // 2, N_T), profile=profile)
print("Propagating from focus to the injection boundary...")
laser.propagate(-x_spot)

env = laser.grid.get_temporal_field()[0, :, :]   # mode 0, complex (n_r, n_t)
r_axis = np.real(laser.grid.axes[0])
t_axis = np.real(laser.grid.axes[1])

# Interpolators on (r, t): amplitude, and phase unwrapped along both axes
amp_interp = RegularGridInterpolator((r_axis, t_axis), np.abs(env),
                                     bounds_error=False, fill_value=0.0)
phase_unwrapped = np.unwrap(np.unwrap(np.angle(env), axis=0), axis=1)
phase_interp = RegularGridInterpolator((r_axis, t_axis), phase_unwrapped,
                                       bounds_error=False, fill_value=0.0)

# Resample onto the full (t, y) grid using r = |y| (cylindrical symmetry)
y_arr = np.linspace(Y_MIN, Y_MAX, N_Y)
t_arr = np.real(t_axis)
Tg, Yg = np.meshgrid(t_arr, y_arr, indexing="ij")
Rg = np.abs(Yg)

amp = amp_interp((Rg, Tg))
amp /= amp.max()

phase_lasy = phase_interp((Rg, Tg))
it_ref = int(np.argmin(np.abs(t_arr - t_centre)))
iy_ref = int(np.argmin(np.abs(y_arr - 0.0)))
phi_ref = phase_lasy[it_ref, iy_ref]                 # the propagator's piston
phase_epoch = -(phase_lasy - phi_ref) + (-gouy)      # piston removed, pinned to deck

print(f"  w_bnd={w_bnd*1e6:.3f} um;  amp peak={amp.max():.4f}")
print(f"  phi_ref(on-axis,t_c)={phi_ref:.4f} rad;  "
      f"phase range=[{phase_epoch.min():.3f}, {phase_epoch.max():.3f}] rad")


def write_dat(path, matrix):
    with open(path, "w") as f:
        f.write(f"{N_T} {N_Y}\n")
        np.savetxt(f, y_arr.reshape(1, -1), fmt="%.8e", delimiter=" ")
        np.savetxt(f, t_arr.reshape(1, -1), fmt="%.8e", delimiter=" ")
        np.savetxt(f, matrix, fmt="%.8e", delimiter=" ")
    print(f"Wrote {path}")


HERE = Path(__file__).parent
for d in ("lasy", "lasy_amp_deck_phase"):
    (HERE / d).mkdir(exist_ok=True)
    write_dat(HERE / d / "amplitude_lasy.dat", amp)
for d in ("lasy", "numerical_amp_lasy_phase"):
    (HERE / d).mkdir(exist_ok=True)
    write_dat(HERE / d / "phase_lasy.dat", phase_epoch)
