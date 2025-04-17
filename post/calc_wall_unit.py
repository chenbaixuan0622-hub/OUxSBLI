import numpy as np
import os
import re
import matplotlib.pyplot as plt
from numba import jit
from mod.mod_read import getGrid, getVector, getScalar, extract_number
from mod.mod_turb_stat import non_dim_tbl, Sutherland
from mod.mod_plot import set_Params


set_Params()

Q_directory = "../3D_solver/TBL/data"
file_path = os.path.join(Q_directory, "Q00500.vtr")
save_path = os.path.join(Q_directory, "yplus.d")
hist_path = os.path.join(Q_directory, "hist.png")
u0    = 506.8e0
gamma = 1.4e0
Rgas  = 287.03e0

Nx, Ny, Nz, x, y, z = getGrid(file_path)

Q = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)

Q[0,:,:,:]                          = getScalar(file_path, Nx, Ny, Nz, 'rho')
Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:]  = getVector(file_path, Nx, Ny, Nz, 'velocity')
Q[4,:,:,:]                          = getScalar(file_path, Nx, Ny, Nz, 'p')
nx1 = Nx // 10
nx2 = Nx // 10 * 4
yp, tw, ut, twm, utm, up = non_dim_tbl(Q[:,:,:,nx1:nx2], x[nx1:nx2], y, z)
Mt = np.mean(ut / np.sqrt(gamma * Q[4,:,0,nx1:nx2] / Q[0,:,0,nx1:nx2]))
save_path = os.path.join(Q_directory, "yp.npy")
np.save(save_path, yp)

print("tw is ", twm)
print("ut is ", utm)
save_path = os.path.join(Q_directory, "yplus.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# yplus  uplus", file=f)
  for j in range(Ny):
    print(f'{yp[j]:.3e}', f'{up[j]:.3e}', file=f)
fig, axs = plt.subplots(1, 2, figsize=(8,4))
axs[0].hist(tw.ravel(), bins=30, density=True, color='blue')
axs[0].set_xlabel(r'$\it{\tau_{w}}$', fontsize=18)
axs[0].set_ylabel(r'$\it{Frequency}$', fontsize=18)

axs[1].hist(ut.ravel(), bins=30, density=True, color='blue')
axs[1].set_xlabel(r'$\it{u_{\tau}}$', fontsize=18)
axs[1].set_ylabel(r'$\it{Frequency}$', fontsize=18)

plt.tight_layout()
plt.savefig(hist_path)
for j in range(Ny):
  if np.mean(Q[1,:,j,0]) >= 0.99e0 * u0:
    tblin = y[j] - (-y[j-1] + y[j]) * \
    (np.mean(Q[1,:,j,0]) - 0.99e0 * u0) \
    / (-np.mean(Q[1,:,j-1,0]) + np.mean(Q[1,:,j,0]) + 1.e-20)
    break
for j in range(Ny):
  if np.mean(Q[1,:,j,nx1:nx2]) >= 0.99e0 * u0:
    tblre = y[j] - (-y[j-1] + y[j]) * \
    (np.mean(Q[1,:,j,nx1:nx2]) - 0.99e0 * u0) \
    / (-np.mean(Q[1,:,j-1,nx1:nx2]) + np.mean(Q[1,:,j,nx1:nx2]) + 1.e-20)
    break
print("tblin is ", tblin*1e3, "[mm]")
print("tblre is ", tblre*1e3, "[mm]")

mu = Sutherland(Q[4,:,0,nx1:nx2] / (Rgas * Q[0,:,0,nx1:nx2]))
save_path = os.path.join(Q_directory, "wall_unit.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# rhow       ut         mu       tau       Mt", file=f)
  print(f'{np.mean(Q[0,:,0,nx1:nx2]):.3e}', f'{utm:.3e}', f'{np.mean(mu):.3e}', f'{twm:.3e}', f'{Mt:.3e}', file=f)

