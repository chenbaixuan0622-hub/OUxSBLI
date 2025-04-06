import numpy as np
import os
from mod.mod_read import getGrid, getScalar

Q_directory = "../../SBLI/SBLI_16delta/2.2ms/9"

file_path = os.path.join(Q_directory, "yp.npy")
yp = np.load(file_path)

file_path           = os.path.join(Q_directory, "ReynoldsStress.vtr")
Nx, Ny, Nz, x, y, z = getGrid(file_path)
ruu = getScalar(file_path, Nx, Ny, Nz, 'ruu')
rvv = getScalar(file_path, Nx, Ny, Nz, 'rvv')
rww = getScalar(file_path, Nx, Ny, Nz, 'rww')
ruv = getScalar(file_path, Nx, Ny, Nz, 'ruv')

Lx1 = 28.e-3
Lx2 = 29.e-3
#nx1 = 0
#nx2 = Nx
for i in range(Nx):
  if x[i] > Lx1:
    nx1 = i
    break
for i in range(Nx):
  if x[i] > Lx2:
    nx2 = i
    break
Lz1 = 4.5e-3
Lz2 = 5.5e-3
for k in range(Nz):
  if z[k] > Lz1:
    nz1 = k
    break
for k in range(Nz):
  if z[k] > Lz2:
    nz2 = k
    break

save_path = os.path.join(Q_directory, "Reynolds_Stress_z4.5_z5.5.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# yp      ruu        rvv       rww        ruv", file=f)
  for j in range(Ny):
    print(f'{yp[j]:.3e}', f'{np.mean(ruu[nz1:nz2,j,nx1:nx2]):.3e}', f'{np.mean(rvv[nz1:nz2,j,nx1:nx2]):.3e}', \
          f'{np.mean(rww[nz1:nz2,j,nx1:nx2]):.3e}', f'{np.mean(ruv[nz1:nz2,j,nx1:nx2]):.3e}', file=f)

