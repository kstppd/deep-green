import numpy as np 
import sys ,os
import matplotlib.pyplot as plt
from multiprocessing import Pool
import h5py as h5
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
plt.style.use('dark_background')


available=[ "rho", "vx", "vy", "vz", "p", "mass", "momentum_x", "momentum_y", "momentum_x", "energy", "temperature"]
units={
        "rho": "\\frac{kg}{m^{3}}",
        "vx" : "m s^{-1}",
        "vy" : "m s^{-1}",
        "vz" : "m s^{-1}",
        "p" :  "kg m s^{-1}", 
        "mass" :"kg",
        "momentum_x" :"kg m s^{-1}",
        "momentum_y": "kg m s^{-1}",
        "momentum_x": "kg m s^{-1}",
        "energy" : "J m^{-3}",
        "temperature":"C"
        }

def plotFile(input):
    var,file,cnt=input
    data=h5.File(file)[var][:]
    print(data.shape)
    nz,ny,nx=np.shape(data)
    data=data[2:-2,2:-2,2:-2]
    nx-=4
    ny-=4
    nz-=4
    print(np.shape(data))
    data= np.flip(data,0)
    #plt.imshow(data[:,ny//2,:],cmap='turbo',norm=colors.LogNorm(vmin=0.7, vmax=1.6))
    #plt.imshow(data[:,ny//2,:],cmap='seismic',norm=colors.LogNorm(vmin=0.5, vmax=2))
    plt.imshow(data[:,ny//2,:],cmap='hot')
    plt.colorbar()
    plt.xlabel("x")
    plt.ylabel("z")
    plt.title(f"{var} [${units[var]}$] ")
    plt.tight_layout()
    plt.savefig(var+"_"+str(cnt).zfill(7)+".png")
    plt.close()


if (len(sys.argv)<3):
    print(f"Usage: python3 {sys.argv[0]} <var> <file sequence>")
    sys.exit(1)

var=sys.argv[1]
if ( not var in available):
    print(f"Invalid var {var} requested!")
    print("\tvar options : rho, vx, vy, vz, p, mass, momentum_x, momentum_y, momentum_x, energy, temperature")
    sys.exit(1)

files=sys.argv[2::]
pool=Pool(8)
index=np.arange(0,len(files))
pool.map(plotFile,zip(np.repeat(var,len(index)),files,index))
