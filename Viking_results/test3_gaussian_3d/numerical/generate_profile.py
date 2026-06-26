"""
Test 3 (3D) profile generator
=============================
Writes ``spatial_profile.dat``: the 2D (y, z) Gaussian spatial envelope at the
x_min boundary, matching the analytical gauss(y,0,w0) * gauss(z,0,w0).

File format (consumed by epoch3d custom_laser_spatial_setup):
    line 1   : n1  n2                 (ny_profile, nz_profile)
    line 2   : y_0 ... y_{n1-1}       [metres]   (first axis)
    line 3   : z_0 ... z_{n2-1}       [metres]   (second axis)
    lines 4+ : one row per z value; each row has n1 values over y
"""

import numpy as np
from pathlib import Path

# --- Parameters: must match input.deck ---
w0 = 3.0e-6
y_min, y_max = -8.0e-6, 8.0e-6
z_min, z_max = -8.0e-6, 8.0e-6

# .dat resolution (>= the deck's ny/nz so the boundary plane is fully resolved)
n1 = 600   # along y
n2 = 600   # along z

y_coords = np.linspace(y_min, y_max, n1)
z_coords = np.linspace(z_min, z_max, n2)

# Separable Gaussian: exp(-(y/w0)^2) * exp(-(z/w0)^2)
profile = np.outer(np.exp(-(y_coords / w0) ** 2),
                   np.exp(-(z_coords / w0) ** 2))   # shape (n1, n2) = (y, z)

out = Path(__file__).parent / "spatial_profile.dat"
with open(out, "w") as f:
    f.write(f"{n1}  {n2}\n")
    f.write("  ".join(f"{y:e}" for y in y_coords) + "\n")
    f.write("  ".join(f"{z:e}" for z in z_coords) + "\n")
    # row j contains profile(:, j) — all n1 y-values at z(j)
    for j in range(n2):
        f.write("  ".join(f"{profile[i, j]:e}" for i in range(n1)) + "\n")

print(f"Wrote {out}")
print(f"  grid: {n1} x {n2};  w0={w0*1e6:.1f} um;  peak={profile.max():.4f}, edge={profile[0,0]:.3e}")
