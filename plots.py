import numpy as np 
import matplotlib.pyplot as plt 
import os 
import pandas as pd

data_dir = './outputs'

dataset = sorted(os.listdir(data_dir))#[::5]

for file in dataset:
    # snapshot_data = pd.read_csv(os.path.join(data_dir, file))
    # snapshot_data = snapshot_data.to_numpy(snapshot_data)
    snapshot_data = np.loadtxt(os.path.join(data_dir, file), 
                               delimiter=",",
                               usecols=[1,2,3,4,5,6])

    x = snapshot_data[:,0]
    y = snapshot_data[:,1]
    z = snapshot_data[:,2]
    vx = snapshot_data[:,3]
    vy = snapshot_data[:,4]
    vz = snapshot_data[:,5]
    
    fig, ax = plt.subplots(1, 1, figsize=(4,4), dpi=300) 
    ax.scatter(x,y, s=0.2)
    ax.set(xlim=(-100,100), ylim=(-100,100))
    
    