import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from matplotlib.colors import LogNorm
import matplotlib
from scipy.spatial.transform import Rotation as R

data_dir = "./outputs"
sequence_dir = "latest_sequence"
end_time = 10  # Gyr
step = 0.01  # Gyr
time = np.arange(0, end_time + step, step)

dataset = sorted(os.listdir(data_dir))
numbins = 200

mid_pause_rotation_interval = np.linspace(0, 2, 150)


def nav_fre_whi(r, rho_0, r_scale):
    rho = rho_0 / ((r / r_scale) * (1 + (r / r_scale)) ** 2)
    return rho


def surface_density(x_coord, y_coord, z_coord, masses, radius, num_bins):
    """
    Gets 3d density profile of a DM halo.
    """

    starting_point = 0.001  # in parsecs
    r = np.geomspace(starting_point, radius, num=num_bins, endpoint=True)

    all_positions = np.vstack((x_coord, y_coord, z_coord)).T

    distances = np.sqrt(np.sum(np.square(all_positions), axis=1))
    mass_per_bin, bin_edges = np.histogram(distances, bins=r, weights=masses)
    count_per_bin, _ = np.histogram(distances, bins=r)

    # mask out empty bins
    mask = count_per_bin > 0
    mass_per_bin = mass_per_bin[mask]
    count_per_bin = count_per_bin[mask]

    # getting bin properties
    right_edges = bin_edges[1:]
    left_edges = bin_edges[:-1]
    bin_ctrs = 0.5 * (left_edges + right_edges)[mask]
    # calculate the mass per thin surface area of a sphere
    d_volume = (4 / 3) * np.pi * (right_edges ** 3 - left_edges ** 3)[mask]
    mass_density = mass_per_bin / d_volume
    avg_star_masses = mass_per_bin / count_per_bin
    err_surf_mass_density = np.sqrt(count_per_bin) * (avg_star_masses / d_volume)

    return bin_ctrs, mass_density, err_surf_mass_density


sn_to_plot = [
    600,
]
with plt.style.context("dark_background"):
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)

    for idx, file in enumerate(dataset):
        print(file)
        if idx not in sn_to_plot:
            continue
        else:
            snapshot_data = pd.read_csv(os.path.join(data_dir, file))
            # snapshot_data = snapshot_data.to_numpy(snapshot_data)
            snapshot_data = np.loadtxt(
                os.path.join(data_dir, file),
                delimiter=",",
                usecols=[0, 1, 2, 3, 4, 5, 6],
            )
            output_number = file.split("_")[1].split(".")[0]
            mass = snapshot_data[:, 0]
            x = snapshot_data[:, 1]
            y = snapshot_data[:, 2]
            z = snapshot_data[:, 3]
            vir_rad = 400
            r, rho, err = surface_density(
                x_coord=x,
                y_coord=y,
                z_coord=z,
                masses=mass,
                radius=vir_rad,
                num_bins=50,
            )
            ax.errorbar(
                r,
                rho,
                yerr=err,
                fmt="o",
                capsize=3,
                capthick=1,
                elinewidth=1,
                ms=3,
                alpha=1,
                c="magenta",
            )
            concentration_param = 4
            rho_not = 1e-2
            r_s = vir_rad / concentration_param
            rho_not = 0.20 * np.sum(mass) / r_s ** 3
            nfw = nav_fre_whi(
                np.geomspace(r.min(), r.max() * 1.5, 100),
                rho_not,
                vir_rad / concentration_param,
            )
            ax.plot(
                np.geomspace(r.min(), r.max() * 1.5, 100),
                nfw,
                label=r"$ R_{{\mathrm{{vir}}}} \:/ \:r_{{\mathrm{{s}}}} = {:.1f} $"
                "\n"
                r"$\rho_0 = {:.3f}$".format(concentration_param, rho_not),
                color="white",
                ls="--",
            )
            ax.text(
                0.05,
                0.05,
                (r"$\mathrm{{t_{{sim}} = {:.2f} \: Gyr}}$").format(time[idx]),
                ha="left",
                va="bottom",
                color="white",
                transform=ax.transAxes,
            )
            ax.legend()
            ax.set(xscale="log", yscale="log", ylabel=r"$\rho$", xlabel="radius")
