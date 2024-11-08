import numpy as np
import os
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl

Q_directory = "../3D_solver/TBL/data"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)

nx0 = int()
nx  = int()
Nm  = int(0.8*Ny)

Q_files.sort(key=extract_number)

rho_path = os.path.join(Q_directory, "rho.npy")
u_path   = os.path.join(Q_directory, "u.npy")
v_path   = os.path.join(Q_directory, "v.npy")
w_path   = os.path.join(Q_directory, "w.npy")
p_path   = os.path.join(Q_directory, "p.npy")

rho = np.load(rho_path)
u   = np.load(u_path)
v   = np.load(v_path)
w   = np.load(w_path)
p   = np.load(p_path)

rhom = rho[:,:,nx0:nx0+nx]
um   =   u[:,:,nx0:nx0+nx]
vm   =   v[:,:,nx0:nx0+nx]
wm   =   w[:,:,nx0:nx0+nx]
pm   =   p[:,:,nx0:nx0+nx]

del rho, u, v, w, p

Q    = np.zeros((5,Nz,Ny,nx), dtype=np.float32)

rho0 = np.mean(rhom[:,Nm,:])
u0   = np.mean(um[:,Nm,:])
p0   = np.mean(pm[:,Nm,:])

# non dim
Q[0,:,:,:], Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:], Q[4,:,:,:] = rhom, um, vm, wm, pm
yp, tw, _, _, ut, up = non_dim_tbl(Q, x, y, z)

# Cf
Cf  = np.mean(2.e0 * tw / (rho0 * u0**2), axis=0)
# p / p0
pp0 = np.mean(pm[:,0,:], axis=0) / p0

for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_directory, Q_file)
  p   = getScalar(file_path, Nx, Ny, Nz, 'p')
  p2 += p**2

# p rms
p2  /= np.float32(num_files)
prms = np.mean(np.sqrt(p2[:,0,:] - pm[:,0,:]**2), axis=0) / p0

save_path = os.path.join(Q_directory, "wall_prop.d")
delta     = 2.e-3
with open(save_path, "w", encoding="UTF-8") as f:
  print("# (x-x0)/d Cf       p / p0      prms", file=f)
  for i in range(nx0, nx0+nx):
    print(f'{(x[i]-x[Nx1])/delta:.3e}', f'{Cf[i]:.3e}', f'{pp0[i]:.3e}', f'{prms[i]:.3e}', file=f)

