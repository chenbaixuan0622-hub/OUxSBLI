import numpy as np
import os
from mod.mod_read import getGrid, getScalar

Q_directory = "../../../../../../media/user/HD-EDS-E/hatayama/TBL/SBLI_05delta/stat03ms_09ms_SLAU"

file_path = os.path.join(Q_directory, "yp.npy")
yp = np.load(file_path)

file_path           = os.path.join(Q_directory, "ReynoldsStress.vtr")
Nx, Ny, Nz, x, y, z = getGrid(file_path)
ruu = getScalar(file_path, Nx, Ny, Nz, 'ruu')
rvv = getScalar(file_path, Nx, Ny, Nz, 'rvv')
rww = getScalar(file_path, Nx, Ny, Nz, 'rww')
ruv = getScalar(file_path, Nx, Ny, Nz, 'ruv')

delta = 25.e0

nx1 = int(21.e0 / delta * Nx)
nx2 = int(23.e0 / delta * Nx)

save_path = os.path.join(Q_directory, "ReynoldsStress_relaxation_region.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# yp      ruu        rvv       rww        ruv", file=f)
  for j in range(Ny):
    print(f'{yp[j]:.3e}', f'{np.mean(ruu[:,j,nx1:nx2]):.3e}', f'{np.mean(rvv[:,j,nx1:nx2]):.3e}', \
          f'{np.mean(rww[:,j,nx1:nx2]):.3e}', f'{np.mean(ruv[:,j,nx1:nx2]):.3e}', file=f)

