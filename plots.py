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

dataset = sorted(os.listdir(data_dir))  # [:3:1]  # [::20]
numbins = 200

mid_pause_rotation_interval = np.linspace(0, 2, 150)
#%%
def draw_frame(x, y, mass, t):
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)
    # ax.scatter(x,y, s=0.2)
    # ax.set(xlim=(-100,100), ylim=(-100,100))

    xy, _, _ = np.histogram2d(
        x,
        y,
        weights=mass,
        bins=numbins,
        normed=False,
        range=[[-100, 100], [-100, 100]],
    )
    surf_density = ax.imshow(
        xy,
        cmap="cubehelix",
        interpolation="gaussian",
        origin="lower",
        extent=[-100, 100, -100, 100],
        # vmin=0,
        # vmax=150,
        norm=LogNorm(),
    )

    ax.set_facecolor(matplotlib.cm.cubehelix(0))

    ax.set(
        xlim=(-80, 80),
        ylim=(-80, 80),
        xticklabels=[],
        yticklabels=[],
    )

    cbar_ax = ax.inset_axes([0.10, 0.10, 0.45, 0.040])
    dens_cbar = fig.colorbar(
        surf_density,
        cax=cbar_ax,
        pad=-1,
        orientation="horizontal",
    )
    cbar_ax.set_title(
        label=r"$\log \: \mathrm{\Sigma}$" "\n",
        fontsize=10,
        pad=-15,
    )
    ax.text(
        0.05,
        0.95,
        (r"$\mathrm{{t_{{sim}} = {:.2f} \: Gyr}}$").format(t),
        ha="left",
        va="top",
        color="white",
        transform=ax.transAxes,
    )
    ax.tick_params(axis="both", labeltop="on", direction="in", which="both")


pause_and_rotate = [175, 375]
with plt.style.context("dark_background"):
    for idx, file in enumerate(dataset):
        print(file)
        if idx not in pause_and_rotate:
            continue
        snapshot_data = pd.read_csv(os.path.join(data_dir, file))
        # snapshot_data = snapshot_data.to_numpy(snapshot_data)
        snapshot_data = np.loadtxt(
            os.path.join(data_dir, file), delimiter=",", usecols=[0, 1, 2, 3, 4, 5, 6]
        )
        output_number = file.split("_")[1].split(".")[0]
        mass = snapshot_data[:, 0]
        x = snapshot_data[:, 1]
        y = snapshot_data[:, 2]
        z = snapshot_data[:, 3]
        # vx = snapshot_data[:, 4]
        # vy = snapshot_data[:, 5]
        # vz = snapshot_data[:, 6]
        if idx in pause_and_rotate:
            # reset the star positions every loop
            print("panning")
            for pan_idx, pause_rot_angle in enumerate(mid_pause_rotation_interval):
                print(pause_rot_angle)
                dm_to_rotate_x = x
                dm_to_rotate_y = y
                dm_to_rotate_z = z
                # new angle each loop
                rotation_angle = pause_rot_angle * np.pi
                # along (x,y,z) axis
                r = R.from_rotvec(rotation_angle * np.array([1, 0, 0]))
                rotation_matrix = r.as_matrix()
                # rotate stars
                rotated_star_positions = np.dot(
                    np.vstack([dm_to_rotate_x, dm_to_rotate_y, dm_to_rotate_z]).T,
                    rotation_matrix.T,
                )
                plot_x, plot_y, plot_z = rotated_star_positions.T
                draw_frame(plot_x, plot_y, mass, time[idx])
                plt.savefig(
                    "./{}/{}_{}.png".format(sequence_dir, output_number, str(pan_idx).zfill(3)),
                    dpi=300,
                    bbox_inches="tight",
                    pad_inches=0.0,
                )
                plt.close()

        else:
            draw_frame(x, y, mass, time[idx])

        plt.savefig(
            "./{}/{}.png".format(sequence_dir, output_number),
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.0,
        )
        # plt.close()
#%%


def nav_fre_whi(r, rho_0, r_scale):
    rho = rho_0 / ((r / r_scale) * (1 + (r / r_scale)) ** 2)
    return rho


def surface_density(x_coord, y_coord, z_coord, masses, radius, num_bins):
    """
    Gets 3d density profile of a DM halo.
    """

    starting_point = 0.01  # in parsecs
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
    d_volume = (4 / 3) * np.pi * (right_edges**3 - left_edges**3)[mask]
    mass_density = mass_per_bin / d_volume
    avg_star_masses = mass_per_bin / count_per_bin
    err_surf_mass_density = np.sqrt(count_per_bin) * (avg_star_masses / d_volume)

    return bin_ctrs, mass_density, err_surf_mass_density


#%%
# #%%
# import glob
# from PIL import Image

# frames_dir = "./images_el/"
# output = "test"

# fp_in = "./{}*.png".format(frames_dir)
# fp_out = "./{}.gif".format(output)
# print(fp_in)
# # https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html#gif
# img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
# print("Processing Frames to a gif...")
# img.save(fp=fp_out, format="GIF", append_images=imgs, save_all=True, duration=4, loop=0)
# print("Image saved to", fp_out)
