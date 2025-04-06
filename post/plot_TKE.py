import numpy as np
import os
from mod.mod_read import getGrid, getScalar

Q_directory = "../../SBLI/SBLI_4delta/2.2ms"

file_path = os.path.join(Q_directory, "yp.npy")
yp = np.load(file_path)

file_path           = os.path.join(Q_directory, "TKE.vtr")
Nx, Ny, Nz, x, y, z = getGrid(file_path)
P   = getScalar(file_path, Nx, Ny, Nz, 'P')
T   = getScalar(file_path, Nx, Ny, Nz, 'T')
PI  = getScalar(file_path, Nx, Ny, Nz, 'PI')
D   = getScalar(file_path, Nx, Ny, Nz, 'D')
eps = getScalar(file_path, Nx, Ny, Nz, 'eps')

Lx1 = 12.5e-3#32.e-3
Lx2 = 18.75e-3#33.e-3
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

Lz1 = 0.e-3
Lz2 = 8.e-3
nz1 = 0
nz2 = Nz
for k in range(Nz):
  if z[k] > Lz1:
    nz1 = k
    break
for k in range(Nz):
  if z[k] > Lz2:
    nz2 = k
    break
save_path = os.path.join(Q_directory, "TKE_Budget_mean_x12.5_18.75.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# yp      P          T         PI         D      eps", file=f)
  for j in range(1,Ny):
    print(f'{yp[j]:.3e}', f'{np.mean(P[nz1:nz2,j-1,nx1:nx2]):.3e}', f'{np.mean(T[nz1:nz2,j-1,nx1:nx2]):.3e}', f'{np.mean(PI[nz1:nz2,j-1,nx1:nx2]):.3e}', \
          f'{np.mean(D[nz1:nz2,j-1,nx1:nx2]):.3e}', f'{np.mean(eps[nz1:nz2,j-1,nx1:nx2]):.3e}', file=f)

