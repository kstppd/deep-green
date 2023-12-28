import numpy as np 
import sys ,os
import matplotlib.pyplot as plt
from multiprocessing import Pool
import h5py as h5
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.style.use('dark_background')


available=[ "rho", "vx", "vy", "vz", "p", "mass", "momentum_x", "momentum_y", "momentum_x", "energy", "temperature"]
units={
        "rho": "kg m^{-3}",
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
    nz,ny,nx=np.shape(data)
    data=np.reshape(data,(nx,ny,nz),order='C')
    data=data[2:-2,2:-2,2:-2]
    nx-=4
    ny-=4
    nz-=4
    print(np.shape(data))
    fig, (ax1, ax2,ax3) = plt.subplots(1, 3)
    fig.suptitle(f"Y and Z slices of {var}")
    im1=ax1.imshow(data[:,:,nz//2],cmap='plasma')#,vmin=1.1,vmax=1.23)
    im2=ax2.imshow(data[:,ny//2,:],cmap='plasma')#,vmin=1.1,vmax=1.23)
    im3=ax3.imshow(data[nx//2,:,:],cmap='plasma')#,vmin=1.1,vmax=1.23)
    ax1.invert_yaxis() 
    ax2.invert_yaxis()
    #ax1.set_title(f"{var}, Z = 0 slice")
    #ax2.set_title(f"{var}, Y = 0 slice")
    ax1.set_xlabel("y")
    ax1.set_ylabel("z")
    ax2.set_xlabel("x")
    ax2.set_ylabel("z")
    ax3.set_xlabel("z")
    ax3.set_ylabel("y")
    divider1 = make_axes_locatable(ax1)
    divider2 = make_axes_locatable(ax2)
    divider3 = make_axes_locatable(ax3)
    cax1 = divider1.append_axes('right', size='5%', pad=0.05)
    cax2 = divider2.append_axes('right', size='5%', pad=0.05)
    cax3 = divider3.append_axes('right', size='5%', pad=0.05)
    cb1=fig.colorbar(im1, cax=cax1)
    cb2=fig.colorbar(im2, cax=cax2)
    cb3=fig.colorbar(im3, cax=cax3)
    cb1.ax.set_title(f"${units[var]}$ ")
    cb2.ax.set_title(f"${units[var]}$ ")
    cb3.ax.set_title(f"${units[var]}$ ")
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
