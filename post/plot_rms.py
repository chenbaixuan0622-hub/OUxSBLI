import numpy as np
import os
from mod.mod_read import getGrid, getScalar, getVector

Q_directory = "../../SBLI/SBLI_16delta/2.2ms/5"

file_path = os.path.join(Q_directory, "yp.npy")
yp = np.load(file_path)

gamma = 1.4e0
rhow  = 1.886e-1
tau   = 1.075e2
Mt    = 6.978e-2

file_path           = os.path.join(Q_directory, "Qrms.vtr")
Nx, Ny, Nz, x, y, z = getGrid(file_path)
rho = getScalar(file_path, Nx, Ny, Nz, 'rho')
p   = getScalar(file_path, Nx, Ny, Nz, 'p')
print(np.min(rho), np.max(rho))
print(np.min(p), np.max(p))
rhorms = rho / (gamma * rhow * Mt**2)
prms   = p / tau

nx1 = 0
nx2 = Nx
'''
Lx1 = 18.75e-3
Lx2 = 25.e-3
for i in range(Nx):
  if x[i] > Lx1:
    nx1 = i
    break
for i in range(Nx):
  if x[i] > Lx2:
    nx2 = i
    break
'''
nz1 = 0
nz2 = Nz
'''
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
'''
save_path = os.path.join(Q_directory, "rms_mean.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# yp      rhorms      prms", file=f)
  for j in range(1,Ny):
    print(f'{yp[j]:.3e}', f'{np.mean(rhorms[nz1:nz2,j-1,nx1:nx2]):.3e}', \
          f'{np.mean(prms[nz1:nz2,j-1,nx1:nx2]):.3e}', file=f)

