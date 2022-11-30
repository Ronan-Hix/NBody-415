import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from matplotlib.colors import LogNorm
import matplotlib


<<<<<<< Updated upstream
data_dir = '/Users/elainetaylor/Desktop/ASTR415/Final Project /NBody-415/outputs'

dataset = sorted(os.listdir(data_dir)) #[::20]

for file in dataset:
    snapshot_data = pd.read_csv(os.path.join(data_dir, file))
    # snapshot_data = snapshot_data.to_numpy(snapshot_data)
    snapshot_data = np.loadtxt(os.path.join(data_dir, file), 
                               delimiter=",",
                               usecols=[1,2,3,4,5,6])
    output_number = file.split("_")[1].split(".")[0]
    x = snapshot_data[:,0]
    y = snapshot_data[:,1]
    z = snapshot_data[:,2]
    vx = snapshot_data[:,3]
    vy = snapshot_data[:,4]
    vz = snapshot_data[:,5]
    
    fig, ax = plt.subplots(1, 1, figsize=(4,4), dpi=300) 
    #ax.scatter(x,y, s=0.2)
    #ax.set(xlim=(-100,100), ylim=(-100,100))
    
    xy, _, _ = np.histogram2d(
                    x,
                    y,
                    bins=500,
                    normed=False,
                    range=[[-100, 100], [-100, 100]],
                    
                )
    ax.imshow(xy,
                    cmap="inferno",
                    interpolation="gaussian",
                    origin="lower",
                    extent=[-100, 100, -100, 100],
                    norm= LogNorm()
                )
    ax.set_facecolor(matplotlib.cm.Greys_r(0))
    ax.set(xlim=(-100,100), ylim=(-100,100))
    plt.savefig("./images/{}.png".format(output_number), dpi=300)
    
=======
data_dir = "./outputs"

dataset = sorted(os.listdir(data_dir))  # [:20:10]  # [::20]
numbins = 200
with plt.style.context("dark_background"):
    for file in dataset:
        snapshot_data = pd.read_csv(os.path.join(data_dir, file))
        # snapshot_data = snapshot_data.to_numpy(snapshot_data)
        snapshot_data = np.loadtxt(
            os.path.join(data_dir, file), delimiter=",", usecols=[1, 2, 3, 4, 5, 6]
        )
        output_number = file.split("_")[1].split(".")[0]
        x = snapshot_data[:, 0]
        y = snapshot_data[:, 1]
        z = snapshot_data[:, 2]
        vx = snapshot_data[:, 3]
        vy = snapshot_data[:, 4]
        vz = snapshot_data[:, 5]

        fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=300)
        # ax.scatter(x,y, s=0.2)
        # ax.set(xlim=(-100,100), ylim=(-100,100))

        xy, _, _ = np.histogram2d(
            x,
            y,
            bins=numbins,
            normed=False,
            range=[[-100, 100], [-100, 100]],
        )
        surf_density = ax.imshow(
            xy,
            cmap="viridis",
            # interpolation="gaussian",
            origin="lower",
            extent=[-100, 100, -100, 100],
            # vmin=0,
            # vmax=150,
            norm=LogNorm(),
        )

        ax.set_facecolor(matplotlib.cm.viridis(0))
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
            pad=-3,
        )

        ax.tick_params(axis="both", labeltop="on", direction="in", which="both")
        plt.savefig(
            "./images/{}.png".format(output_number),
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.0,
        )
        plt.close()

#%%
import glob
from PIL import Image

frames_dir = "./images/"
output = "test"

fp_in = "./{}*.png".format(frames_dir)
fp_out = "./{}.gif".format(output)
print(fp_in)
# https://pillow.readthedocs.io/en/stable/handbook/image-file-formats.html#gif
img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
print("Processing Frames to a gif...")
img.save(fp=fp_out, format="GIF", append_images=imgs, save_all=True, duration=4, loop=0)
print("Image saved to", fp_out)
>>>>>>> Stashed changes
