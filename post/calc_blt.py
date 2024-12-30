import numpy as np
import os
from tqdm import tqdm
from mod.mod_read import getGrid, getVector


Q_dir  = "../../a100"
Q_path = os.path.join(Q_dir, "Qmean.vtr")
Lzs    = [5.75e-3, 6.25e-3]
u0     = 506.8e0

Nx, Ny, Nz, x, y, z = getGrid(Q_path)
U, _, _ = getVector(Q_path, Nx, Ny, Nz, 'velocity')

nz1 = int(Lzs[0] / z[-1] * Nz)
nz2 = int(Lzs[1] / z[-1] * Nz)

blt = np.zeros(Nx)

# calc boundary layer thickness
for i in tqdm(range(Nx)):
  for j in range(Ny):
    u = np.mean(U[nz1:nz2,:,i], axis=0)
    if u[j] >= 0.99e0 * u0:
      blt[i] = y[j] - (-y[j-1] + y[j]) * (u[j] - 0.99e0 * u0) / (-u[j-1] + u[j])
      break


name = "blt_" + "z" + f'{Lzs[0]*1e3:.2f}' + "_" + f'{Lzs[1]*1e3:.2f}'
save_path = os.path.join(Q_dir, name + ".d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# x[mm]   blt[mm]", file=f)
  for i in range(Nx):
    print(f'{x[i]*1e3:.3e}', f'{blt[i]*1e3:.3e}', file=f)

