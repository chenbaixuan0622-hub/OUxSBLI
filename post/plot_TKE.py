import numpy as np
import os
from mod.mod_read import getGrid, getScalar

Q_directory = "../3D_solver/TBL/SBLI_05delta/stat03ms_09ms"

file_path = os.path.join(Q_directory, "yp.npy")
yp = np.load(file_path)

file_path           = os.path.join(Q_directory, "TKE.vtr")
Nx, Ny, Nz, x, y, z = getGrid(file_path)
P   = getScalar(file_path, Nx, Ny, Nz, 'P')
T   = getScalar(file_path, Nx, Ny, Nz, 'T')
PI  = getScalar(file_path, Nx, Ny, Nz, 'PI')
D   = getScalar(file_path, Nx, Ny, Nz, 'D')
eps = getScalar(file_path, Nx, Ny, Nz, 'eps')

delta = 25.e0

nx1 = int(21.e0 / delta * Nx)
nx2 = int(23.e0 / delta * Nx)

save_path = os.path.join(Q_directory, "TKE_Budget_relaxation_region.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# yp      P          T         PI         D      eps", file=f)
  for j in range(1,Ny):
    print(f'{yp[j]:.3e}', f'{np.mean(P[:,j-1,nx1:nx2]):.3e}', f'{np.mean(T[:,j-1,nx1:nx2]):.3e}', f'{np.mean(PI[:,j-1,nx1:nx2]):.3e}', \
          f'{np.mean(D[:,j-1,nx1:nx2]):.3e}', f'{np.mean(eps[:,j-1,nx1:nx2]):.3e}', file=f)

