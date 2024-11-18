import numpy as np
import os
from mod.mod_read import getGrid, getVector, getScalar

Q_directory = "../3D_solver/TBL/data"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)

# parameters
Rgas  = 287.03e0
delta = 2.e-3

rho_path = os.path.join(Q_directory, "rho.npy")
u_path   = os.path.join(Q_directory, "u.npy")
p_path   = os.path.join(Q_directory, "p.npy")

rhom = np.load(rho_path)
um   = np.load(u_path)
pm   = np.load(p_path)

nx1 = int(0.5*Nx)
nx2 = int(0.9*Nx)
ny1 = int(0.8*Ny)

rho0 = np.mean(rhom[:,ny1,nx1:nx2])
u0   = np.mean(um[:,ny1,nx1:nx2])
T0   = np.mean(pm[:,ny1,nx1:nx2] / (Rgas * rhom[:,ny1,nx1:nx2]))

rhou = np.zeros(Ny)
u    = np.zeros(Ny)
T    = np.zeros(Ny)

print(np.mean(um[:,0,:]))

for j in range(Ny):
  rhou[j] = np.mean(rhom[:,j,nx1:nx2] * um[:,j,nx1:nx2]) / (rho0 * u0)
  u[j]    = np.mean(um[:,j,nx1:nx2]) / u0
  T[j]    = np.mean(pm[:,j,nx1:nx2] / (Rgas * rhom[:,j,nx1:nx2])) / T0

save_path = os.path.join(Q_directory, "turb_mean.d")

with open(save_path, "w", encoding="UTF-8") as f:
  print("# y      rhou       u          T", file=f)
  for j in range(Ny):
    print(f'{y[j]/delta:.3e}', f'{rhou[j]:.3e}', f'{u[j]:.3e}', f'{T[j]:.3e}', file=f)

