# This is a simple script to generate a spatial profile for the laser profile wrapper in the EPOCH simulation. The profile is a super-Gaussian shape, which is commonly used to model laser beams with a flat top and steep edges. The generated profile will be saved in a file called "spatial_profile.dat" in the format that EPOCH expects.

import numpy as np

# 1. Define configuration
num_points = 600
y_min = -10.0e-6  # -10 microns in meters
y_max = 10.0e-6   # 10 microns in meters
beam_waist = 3.5e-6  # Profile starts dropping off around 3.5 microns
order = 6  # Higher order = flatter top, steeper edges

# 2. Generate coordinates and the Super-Gaussian shape
y_coords = np.linspace(y_min, y_max, num_points)
# Formula: exp(-2 * (y/w)^order)
profile_values = np.exp(-2 * (y_coords / beam_waist) ** order)

# 3. Force exact zeros at the absolute boundaries to prevent clipping fields
profile_values[0] = 0.0
profile_values[-1] = 0.0

# 4. Write to spatial_profile.dat in the format EPOCH expects
with open("spatial_profile.dat", "w") as f:
    # First line must be the total number of entries
    f.write(f"{num_points}\n")
    
    # Subsequent lines are: coordinate  value
    for y, val in zip(y_coords, profile_values):
        f.write(f"{y:e}  {val:e}\n")

print(f"Successfully generated spatial_profile.dat with {num_points} points!")