import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from matplotlib.colors import LogNorm
import matplotlib
from scipy.spatial.transform import Rotation as R

data_dir = "./outputs"
sequence_dir = "multi_panel"
end_time = 10  # Gyr
step = 0.01  # Gyr
time = np.arange(0, end_time + step, step)

dataset = sorted(os.listdir(data_dir))  # [:800:50]  # [::20]
numbins = 200


def draw_frame(x, y, mass, t, axis_object):

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
    surf_density = axis_object.imshow(
        xy,
        cmap="cubehelix",
        interpolation="gaussian",
        origin="lower",
        extent=[-100, 100, -100, 100],
        # vmin=0,
        # vmax=150,
        norm=LogNorm(),
    )

    axis_object.set_facecolor(matplotlib.cm.cubehelix(0))

    axis_object.set(
        xlim=(-80, 80),
        ylim=(-80, 80),
        xticklabels=[],
        yticklabels=[],
    )

    cbar_ax = axis_object.inset_axes([0.10, 0.13, 0.45, 0.040])
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
    axis_object.text(
        0.05,
        0.95,
        (r"$\mathrm{{t_{{sim}} = {:.2f} \: Gyr}}$").format(t),
        ha="left",
        va="top",
        color="white",
        transform=axis_object.transAxes,
    )
    axis_object.tick_params(axis="both", labeltop="on", direction="in", which="both")


sn_to_plot = [0, 100, 150, 200, 250, 300, 350, 400, 450, 600, 650, 900]
with plt.style.context("dark_background"):
    nrows = 4
    ncols = 3
    fig, ax = plt.subplots(nrows, ncols, figsize=(8, 12), dpi=300)
    ax = ax.ravel()
    plt_idx = 0
    for idx, file in enumerate(dataset):
        # if idx > nrows * ncols:
        #     continue

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
            # vx = snapshot_data[:, 4]
            # vy = snapshot_data[:, 5]
            # vz = snapshot_data[:, 6]

            draw_frame(x, y, mass, time[idx], ax[plt_idx])
            plt_idx += 1
            plt.subplots_adjust(hspace=0, wspace=0)

    plt.savefig(
        "./{}/{}.png".format(sequence_dir, output_number),
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.0,
    )
    # plt.close()
