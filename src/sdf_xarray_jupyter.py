import sdf_xarray as sdfxr
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from IPython.display import HTML
import glob
import os
from matplotlib.animation import FuncAnimation
import matplotlib.colors as colors

def read_multiple_sdf(file_parent_path, read_range=None, data_vars=None, convert_units=True, **kwargs):
    '''
    Read multiple SDF files into a single xarray Dataset using sdf_xarray's open_mfdataset.
    Return ds, which is an xarray Dataset containing all the data from the specified SDF files.
    Parameters:
    - file_parent_path: The directory containing the SDF files.
    - read_range: Optional tuple (lower, upper) to specify a range of files to read (e.g., (0, 100) to read files 0000.sdf to 0099.sdf). If None, all .sdf files in the directory will be read.
    - data_vars: Optional list of variable names to load. If None, all variables will be loaded.
    - convert_units: If True, automatically rescale time to femtoseconds and spatial coordinates to micrometers for easier interpretation. Set to False to keep original units.
    - **kwargs: Additional keyword arguments to pass to open_mfdataset (e.g., chunks={'time': 10} for Dask parallelisation).
    '''
    
    print("Starting to load files...")
    
    if read_range is None:
        file_pattern = os.path.join(file_parent_path, "*.sdf")
        files_to_load = sorted(glob.glob(file_pattern))
    else:
        lower, upper = read_range
        # Build the list and filter out non-existent files in one go
        files_to_load = [
            os.path.join(file_parent_path, f"{i:04d}.sdf") 
            for i in range(lower, upper)
        ]
        files_to_load = [f for f in files_to_load if os.path.exists(f)]

    if not files_to_load:
        raise FileNotFoundError("No files found matching your criteria. Check your paths.")
    
    print(f"Attempting to load {len(files_to_load)} files...")
    ds = sdfxr.open_mfdataset(files_to_load, data_vars=data_vars, **kwargs)
    
    if convert_units:
        ds = ds.epoch.rescale_coords(1e15, "fs", "time")
        ds = ds.epoch.rescale_coords(1e6, "µm", ["X_Grid_mid", "Y_Grid_mid"])
    
    print("Files loaded successfully!")
    return ds

def simple_animation(ds, variable_name, save_parent_path=None, convert_units=False):
    '''
    Create a simple animation of a 2D variable across time using xarray's built-in animation capabilities.
    Parameters:
    - ds: The xarray Dataset containing the data.
    - variable_name: The name of the variable to animate.
    - save_path: Optional path to save the animation. If None, the animation will be displayed.
    - convert_units: If True, automatically rescale time to femtoseconds and spatial coordinates to micrometers for easier interpretation. Set to False to keep original units.
    '''
    if convert_units:
        # Chain the rescales cleanly
        ds_plot = ds.epoch.rescale_coords(1e15, "fs", "time")
        ds_plot = ds_plot.epoch.rescale_coords(1e6, "µm", ["X_Grid_mid", "Y_Grid_mid"])
        da = ds_plot[variable_name]
    else:
        da = ds[variable_name]
    
    anim = da.epoch.animate()
    
    if save_parent_path:
        save_path = os.path.join(save_parent_path, f"{variable_name}_animation.mp4")
        anim.save(save_path)
        print(f"Animation saved to {save_path}")
    else:
        anim.show() # plt.show() might be better?
        
    return anim # Returning the object is highly recommended if you use Jupyter Notebooks

def plot_laser_abs_frac(ds, save_parent_path=None, save_name="laser_absorption_fraction.png",does_save=True):

    ds["Laser_Absorption_Fraction_in_Simulation"] = (
        (ds["Total_Particle_Energy_in_Simulation"] - ds["Total_Particle_Energy_in_Simulation"][0])
        / ds["Absorption_Total_Laser_Energy_Injected"]
        ) *100

    # We can also manipulate the units and other attributes
    ds["Laser_Absorption_Fraction_in_Simulation"].attrs["units"] = "%"
    ds["Laser_Absorption_Fraction_in_Simulation"].attrs["long_name"] = "Laser Absorption Fraction"

    ds["Laser_Absorption_Fraction_in_Simulation"].epoch.plot()
    plt.title("Laser absorption fraction in simulation")

    if does_save:
        if save_parent_path is not None:
            save_path = os.path.join(save_parent_path, save_name)
        else:
            save_path = os.path.join(os.getcwd(), save_name)
        plt.savefig(save_path)
        print(f"Plot saved to {save_path}")
        
    plt.show()

def convert_energy_to_MeV(ds, variable_name):
    '''
    Convert an energy variable in the dataset from its original units to MeV using pint for unit handling.
    The function assumes that the variable has proper units attached. If not, it will raise an error.
    '''
    energy_da = ds[variable_name].pint.quantify()  # Quantify with pint
    energy_MeV = energy_da.pint.to("MeV")  # Convert to MeV
    energy_MeV = energy_MeV.pint.dequantify(format="~P")  # Remove units for easier plotting
    return energy_MeV

def animate_energy_spectrum(ds, dist_type="en", species="Electron", time_is_fs=True):
    '''
    Animates the energy spectrum (dN/dE).
    
    Parameters:
    - ds: The xarray Dataset.
    - dist_type: The dimensions of the dist_fn (e.g., "en" or "px_en").
    - species: The species name (e.g., "Electron", "pos_lbw").
    - time_is_fs: Set to True if ds['time'] is already in femtoseconds.
    '''
    # 1. Dynamic Variable Construction
    var_name = f"dist_fn_{dist_type}_{species}"
    if var_name not in ds.data_vars:
        raise KeyError(f"Variable '{var_name}' not found in dataset. Check dist_type and species.")
    
    data_fn = ds[var_name]

    # 2. Handling 1D vs 2D (Summing over momentum if necessary)
    # Search for any dimension containing 'px'
    px_dims = [d for d in data_fn.dims if 'px' in d.lower()]
    
    if px_dims:
        print(f"Detected 2D distribution. Summing over {px_dims[0]}...")
        dn_de = data_fn.sum(dim=px_dims[0])
    else:
        print("Detected 1D distribution (Energy only).")
        dn_de = data_fn

    # 3. Identify the Energy Dimension and Convert Units
    # Finds the dimension string that contains 'energy'
    en_dim = [d for d in dn_de.dims if 'energy' in d.lower()][0]
    
    joules_to_mev = 1.0 / 1.602176634e-13
    energy_mev = dn_de[en_dim].values * joules_to_mev
    
    # 4. Handle Time Units
    time_vals = dn_de['time'].values
    time_label = "fs" if time_is_fs else "s"

    # 5. Setup Plotting
    fig, ax = plt.subplots(figsize=(8, 6))
    line, = ax.plot([], [], lw=2, color='blue')
    
    # Set Limits
    ax.set_xlim(energy_mev.min(), energy_mev.max())
    y_max = float(dn_de.max().values)
    if y_max <= 0: y_max = 1e-10 # Prevent axis errors if empty
    ax.set_ylim(0, y_max * 1.1)

    ax.set_xlabel("Energy (MeV)", fontsize=12)
    ax.set_ylabel(r"$dN/dE$ (Arbitrary Units)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.6)
    title = ax.set_title(f"{species} Energy Spectrum", fontsize=14)

    # 6. Animation Logic
    def init():
        line.set_data([], [])
        title.set_text("")
        return line, title

    def update(frame_idx):
        y_data = dn_de.isel(time=frame_idx).values
        line.set_data(energy_mev, y_data)
        title.set_text(f"{species} Energy Spectrum | Time: {time_vals[frame_idx]:.1f} {time_label}")
        return line, title

    anim = FuncAnimation(
        fig, update, frames=len(time_vals), 
        init_func=init, blit=False, interval=100
    )

    plt.close(fig) # Prevents the extra static plot from appearing
    return HTML(anim.to_jshtml())


def animate_energy_direction(ds, dist_type="px_en", species="Electron", time_is_fs=True, split_direction=False):
    '''
    Animates the energy spectrum. If split_direction=True and data is 2D, 
    it plots forward (+px) and backward (-px) spectra separately.
    '''
    # 1. Variable Setup
    var_name = f"dist_fn_{dist_type}_{species}"
    if var_name not in ds.data_vars:
        raise KeyError(f"Variable '{var_name}' not found. Did you output {dist_type}?")
    
    data_fn = ds[var_name]
    px_dims = [d for d in data_fn.dims if 'px' in d.lower()]
    en_dim = [d for d in data_fn.dims if 'energy' in d.lower()][0]
    
    joules_to_mev = 1.0 / 1.602176634e-13
    energy_mev = data_fn[en_dim].values * joules_to_mev
    time_vals = data_fn['time'].values
    time_unit = "fs" if time_is_fs else "s"

    # 2. Logic for Splitting Direction
    if split_direction and px_dims:
        px_dim = px_dims[0]
        # Slice: Forward (px > 0) and Backward (px < 0)
        # sel() with slice(None, 0) gets all values up to 0
        dn_de_fwd = data_fn.sel({px_dim: slice(0, None)}).sum(dim=px_dim)
        dn_de_bwd = data_fn.sel({px_dim: slice(None, 0)}).sum(dim=px_dim)
        mode = "split"
    else:
        # Standard total sum
        dn_de_total = data_fn.sum(dim=px_dims[0]) if px_dims else data_fn
        mode = "total"

    # 3. Setup Plotting
    fig, ax = plt.subplots(figsize=(8, 6))
    
    if mode == "split":
        line_fwd, = ax.plot([], [], lw=2, color='crimson', label='Forward ($+p_x$)')
        line_bwd, = ax.plot([], [], lw=2, color='royalblue', label='Backward ($-p_x$)')
        y_max = max(float(dn_de_fwd.max()), float(dn_de_bwd.max())) # pyright: ignore[reportPossiblyUnboundVariable]
        ax.legend(frameon=False)
    else:
        line_total, = ax.plot([], [], lw=2, color='black', label='Total')
        y_max = float(dn_de_total.max()) # pyright: ignore[reportPossiblyUnboundVariable]

    # Styling
    ax.set_xlim(energy_mev.min(), energy_mev.max())
    ax.set_ylim(0, (y_max if y_max > 0 else 1e-10) * 1.1)
    ax.set_xlabel("Energy (MeV)", fontsize=12)
    ax.set_ylabel(r"$dN/dE$", fontsize=12)
    ax.grid(True, alpha=0.3)
    title = ax.set_title(f"{species} Directional Spectrum", fontsize=14)

    # 4. Animation Logic
    def init():
        if mode == "split":
            line_fwd.set_data([], []) # pyright: ignore[reportPossiblyUnboundVariable]
            line_bwd.set_data([], []) # pyright: ignore[reportPossiblyUnboundVariable]
            return line_fwd, line_bwd, title # pyright: ignore[reportPossiblyUnboundVariable]
        else:
            line_total.set_data([], []) # pyright: ignore[reportPossiblyUnboundVariable]
            return line_total, title # pyright: ignore[reportPossiblyUnboundVariable]

    def update(frame_idx):
        if mode == "split":
            line_fwd.set_data(energy_mev, dn_de_fwd.isel(time=frame_idx).values) # pyright: ignore[reportPossiblyUnboundVariable]
            line_bwd.set_data(energy_mev, dn_de_bwd.isel(time=frame_idx).values) # pyright: ignore[reportPossiblyUnboundVariable]
            title.set_text(f"{species} Energy | Time: {time_vals[frame_idx]:.1f} {time_unit}")
            return line_fwd, line_bwd, title # pyright: ignore[reportPossiblyUnboundVariable]
        else:
            line_total.set_data(energy_mev, dn_de_total.isel(time=frame_idx).values) # pyright: ignore[reportPossiblyUnboundVariable]
            title.set_text(f"{species} Energy | Time: {time_vals[frame_idx]:.1f} {time_unit}")
            return line_total, title # pyright: ignore[reportPossiblyUnboundVariable]

    anim = FuncAnimation(fig, update, frames=len(time_vals), init_func=init, blit=False, interval=100)
    plt.close(fig)
    return HTML(anim.to_jshtml())

def animate_2d_dist(ds, dist_type="px_py", species="Photon", x_label=None, y_label=None, plot_title=None, time_is_fs=True):
    '''
    Animates a 2D distribution function heatmap over time.
    Automatically converts SI momentum to MeV/c for easier interpretation.
    The canonical example is px_py for photons. 
    '''
    # 1. Setup variable and dimensions
    var_name = f"dist_fn_{dist_type}_{species}"
    if var_name not in ds.data_vars:
        raise KeyError(f"Variable '{var_name}' not found.")
    
    data = ds[var_name]
    dims = data.dims[1:]  # Skip 'time'
    
    # 2. Convert SI Momentum to MeV/c
    # p [MeV/c] = (p [kg*m/s] * c) / q_e / 1e6
    c = 299792458.0
    mev_conv = c / (1.602176634e-13)
    
    x_axis = data[dims[1]].values * mev_conv
    y_axis = data[dims[0]].values * mev_conv
    time_vals = data['time'].values
    time_unit = "fs" if time_is_fs else "s"

    # 3. Setup Plotting
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Calculate global max for a stable colorbar
    v_max = float(data.max())
    v_min = v_max * 1e-6 if v_max > 0 else 1e-10

    # Initialize heatmap with LogNorm
    # We use a dummy array for init; origin='lower' is critical for momentum plots
    quad = ax.pcolormesh(x_axis, y_axis, np.zeros((len(y_axis), len(x_axis))), 
                         norm=colors.LogNorm(vmin=v_min, vmax=v_max),
                         shading='auto', cmap='magma')
    
    plt.colorbar(quad, ax=ax, label='Density (arb. units)')
    
    ax.set_xlabel(x_label if x_label else f"{dims[1]} (MeV/c)")
    ax.set_ylabel(y_label if y_label else f"{dims[0]} (MeV/c)")
    title = ax.set_title(plot_title if plot_title else f"{species} {dist_type} Distribution", fontsize=14)

    # 4. Animation Logic
    def update(frame_idx):
        # Transpose might be needed depending on how EPOCH/xarray orders the grid
        frame_data = data.isel(time=frame_idx).values
        quad.set_array(frame_data.ravel())
        title.set_text(f"{species} {dist_type} | Time: {time_vals[frame_idx]:.1f} {time_unit}")
        return quad, title

    anim = FuncAnimation(fig, update, frames=len(time_vals), blit=False, interval=120)
    
    plt.close(fig)
    return HTML(anim.to_jshtml())