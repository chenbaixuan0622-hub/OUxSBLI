import numpy as np
import os
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl

Q_directory = "../../a100"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)

nx0 = int(0.1*Nx)
nx  = int(0.2*Nx)
Nm  = int(0.8*Ny)

rho_path = os.path.join(Q_directory, "rho.npy")
u_path   = os.path.join(Q_directory, "u.npy")
v_path   = os.path.join(Q_directory, "v.npy")
w_path   = os.path.join(Q_directory, "w.npy")
p_path   = os.path.join(Q_directory, "p.npy")

rhom = np.load(rho_path)
um   = np.load(u_path)
vm   = np.load(v_path)
wm   = np.load(w_path)
pm   = np.load(p_path)

rho0 = np.mean(rhom[:,Nm,nx0:nx0+nx])
u0   = np.mean(  um[:,Nm,nx0:nx0+nx])
p0   = np.mean(  pm[:,Nm,nx0:nx0+nx])

# non dim
Q    = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
Q[0,:,:,:], Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:], Q[4,:,:,:] = rhom, um, vm, wm, pm
yp, tw, _, _, ut, up = non_dim_tbl(Q, x, y, z)
del Q

# Cf
Cf  = np.mean(2.e0 * tw / (rho0 * u0**2), axis=0)
# p / p0
pp0 = np.mean(pm[:,0,:], axis=0) / p0

itr = 0.e0
p2 = np.zeros((Nz,Ny,Nx), dtype=np.float32)
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_directory, Q_file)
  if file_path == os.path.join(Q_directory, "TKE.vtr"):
    continue
  p   = getScalar(file_path, Nx, Ny, Nz, 'p')
  p2 += p**2
  itr += 1.e0

# p rms
p2  /= itr
prms = np.mean(np.sqrt(p2[:,0,:] - pm[:,0,:]**2), axis=0) / p0

save_path = os.path.join(Q_directory, "wall_prop.d")
delta     = 2.e-3
with open(save_path, "w", encoding="UTF-8") as f:
  print("# x/d      Cf       p / p0      prms", file=f)
  for i in range(Nx):
    print(f'{x[i]/delta:.3e}', f'{Cf[i]:.3e}', f'{pp0[i]:.3e}', f'{prms[i]:.3e}', file=f)

