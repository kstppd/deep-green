import os,sys
import matplotlib.pyplot as plt
import h5py as h5
import numpy as np
import scipy.ndimage as ndimage


def get_drag_force(p,sdf,flow=np.array([1, 0, 0]),ds=1.0):
    nx,ny,nz=p.shape
    surface = np.argwhere(np.abs(sdf) <= 1)
    print(len(surface))
    grad_px, grad_py, grad_pz = np.gradient(sdf.astype(float),ds,ds,ds)
    norm_magnitude = np.sqrt(grad_px**2 + grad_py**2 + grad_pz**2)
    normal_x = np.zeros_like(sdf)
    normal_y = np.zeros_like(sdf)
    normal_z = np.zeros_like(sdf)
    
    valid = norm_magnitude > 1e-8 
    normal_x[valid] = grad_px[valid] / norm_magnitude[valid]
    normal_y[valid] = grad_py[valid] / norm_magnitude[valid]
    normal_z[valid] = grad_pz[valid] / norm_magnitude[valid]
    
    
    dA = ds * ds  
    F_D = 0.0
    for i, j, k in surface:
        phere = p[i, j, k]
        n = np.array([normal_x[i, j, k], normal_y[i, j, k], normal_z[i, j, k]])
        force_here = phere * np.dot(n, flow) * dA
        F_D += force_here
    return F_D

def get_drag_coeff(file,upstream_vel,upstream_rho,radious):
    rho=h5.File(file)["primitives/rho"][:]
    vx,vy,vz=h5.File(file)["primitives/vx"][:],h5.File(file)["primitives/vy"][:],h5.File(file)["primitives/vz"][:]
    p=h5.File(file)["primitives/p"][:]
    sdf=h5.File(file)["sdf_object"][:]
    count = np.sum(np.abs(sdf) < 1)
    fd=get_drag_force(p,sdf)
    print("--->",fd)
    S =np.pi * radious**2
    cd=2.0*fd/(upstream_rho*upstream_vel*upstream_vel*S)
    print(f"Cd={cd}")

files = sys.argv[1::]

for file in files:
    get_drag_coeff(file,10.0,1.225,32.0)


