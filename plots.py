import numpy as np 
import matplotlib.pyplot as plt 
import os 
import pandas as pd
from  matplotlib.colors import LogNorm
import matplotlib
import glob 
import PIL
from PIL import Image


data_dir = './outputs'

dataset = sorted(os.listdir(data_dir))[::20]

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

fp_in = "./{}*.png".format("./images/")
fp_out = "./{}_{}.gif".format(x, y)
print(fp_in)
img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
print("Processing Frames to a gif...")
img.save(
     fp=fp_out, format="GIF", append_images=imgs, save_all=True, duration=4, loop=0
 )
print("Image saved to", fp_out)