import numpy as np
import os
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl, Sutherland

#Q_directory = "../../../../../mnt/data1/TBL_HRSLAU2_LES"
Q_directory = "../3D_solver/TBL/data"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)

# parameters
Rgas  = 287.03e0
gamma = 1.4e0
delta = 0.5e-3

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

nx1 = int(0.1*Nx)
nx2 = int(0.4*Nx)
Q = np.zeros((5,Nz,Ny,nx2-nx1), dtype=np.float32)

Nm   = int(0.8*Ny)
rhow = np.mean(rho[:,0, nx1:nx2])
rho0 = np.mean(rho[:,Nm,nx1:nx2])
u0   = np.mean(u[:,Nm,nx1:nx2])
Tw   = np.mean(p[:,0, nx1:nx2] / (Rgas * rho[:,0, nx1:nx2]))
T0   = np.mean(p[:,Nm,nx1:nx2] / (Rgas * rho[:,Nm,nx1:nx2]))
muw  = Sutherland(Tw)
mu0  = Sutherland(T0)

# non dim
Q[0,:,:,:], Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:], Q[4,:,:,:] = rho[:,:,nx1:nx2], u[:,:,nx1:nx2], v[:,:,nx1:nx2], w[:,:,nx1:nx2], p[:,:,nx1:nx2]
yp, _, _, tw, ut, up = non_dim_tbl(Q, x[nx1:nx2], y, z)

theta = 0.e0
for j in range(1,Ny):
  theta += np.mean(rho[:,j,nx1:nx2] * u[:,j,nx1:nx2]) / (rho0 * u0) * (1.e0 - np.mean(u[:,j,nx1:nx2]) / u0) * (-y[j-1] + y[j])

Retau    = rhow * ut * delta / muw
Retheta  = rho0 * u0 * theta / mu0
Redelta2 = rho0 * u0 * theta / muw
Cf       = 2.e0 * tw / (rho0 * u0**2)
Mt       = ut / np.sqrt(gamma * Rgas * Tw)

save_path = os.path.join(Q_directory, "Re.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# Retau     Retheta       Redelta2     Cf       Mt", file=f)
  print(f'{Retau:.3e}', f'{Retheta:.3e}', f'{Redelta2:.3e}', f'{Cf:.3e}', f'{Mt:.3e}', file=f)

