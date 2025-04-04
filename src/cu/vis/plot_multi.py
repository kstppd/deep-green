import numpy as np 
import sys
import os
import matplotlib.pyplot as plt
from multiprocessing import Pool
import h5py as h5
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors


available = [
    "primitives/rho", "primitives/vx", "primitives/vy", 
    "primitives/vz", "primitives/pressure", "conserved/mass", 
    "conserved/energy"
]

units = {
    "primitives/rho": r"$\frac{kg}{m^{3}}$",
    "primitives/vx": r"$m \, s^{-1}$",
    "primitives/vy": r"$m \, s^{-1}$",
    "primitives/vz": r"$m \, s^{-1}$",
    "primitives/pressure": r"$kg \, m \, s^{-1}$", 
    "conserved/mass": r"$kg$",
    "conserved/energy": r"$J \, m^{-3}$",
}

limits = {
    "primitives/rho": [0.7, 1.6],
    "primitives/vx": [None, None],
    "primitives/vy": [None, None],
    "primitives/vz": [None, None],
    "primitives/pressure": [None, None], 
    "conserved/mass": [None, None],
    "conserved/energy": [None, None]
}        

def read_var(file, var):
    with h5.File(file, 'r') as f:
        data = f[var][:]
    
    nz, ny, nx = np.shape(data)
    data = data[2:-2, 2:-2, 2:-2]  
    return data

def plotFile(input):
    file, cnt = input
    name = f"{file}_{str(cnt).zfill(7)}.png"
    vars = [
        "primitives/rho",
        "primitives/vx",
        "primitives/vy",
        "primitives/vz",
        "primitives/p",
        "sdf_object",
    ]

    data = [read_var(file, v) for v in vars]
    nz, ny, nx = np.shape(data[0])

    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    axes = axes.flatten()
    inside_object = np.where((data[5]) <=1)

    for i in range(len(data) - 1):  # Skip sdf_object itself
        data[i][inside_object] = np.nan

    for i in range(len(vars)):
        if (vars[i]=="primitives/p"):
            im = axes[i].imshow(data[i][:, ny//2, :].T, cmap='viridis',norm=colors.LogNorm())
        else:
            im = axes[i].imshow(data[i][:, ny//2, :].T, cmap='viridis')
        axes[i].set_title(vars[i])
        divider = make_axes_locatable(axes[i])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
    
    for i in range(len(vars), len(axes)):  
        axes[i].axis('off')  

    plt.tight_layout()
    plt.savefig(name, dpi=200)
    plt.close(fig) 
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: python3 {sys.argv[0]} <file sequence>")
        sys.exit(1)

    files = sys.argv[1:]
    pool = Pool(8)
    pool.map(plotFile, zip(files, range(len(files))))
