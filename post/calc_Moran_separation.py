import numpy as np
import os
from libpysal.weights import lat2W
from esda.moran import Moran
import matplotlib.pyplot as plt
from tqdm import tqdm
from mod.mod_read import getGrid, getScalar, getVector, extract_number


# parameter
u_dir = "../../../../../mnt/data2/SBLI_16delta_2mm/u"
delta = 0.5e-3
Lx1   = 14.5e0 * delta
Lx2   = 15.5e0 * delta
# wall unit
rhow  = 0.1895e0
ut    = 27.41e0
mu    = 1.72e-5
nu    = mu / rhow

files = [f for f in os.listdir(u_dir) if f.endswith(".vtr")]
files.sort(key=extract_number)
Nx, Ny, Nz, x, y, z = getGrid(os.path.join(u_dir, files[0]))
nx1 = int(Nx * (-x[0] + Lx1) / (-x[0] + x[-1]))
nx2 = int(Nx * (-x[0] + Lx2) / (-x[0] + x[-1]))
for j in range(Ny):
  yp = ut * y[j] / nu
  if yp > 10.e0:
    ny = j
    break

dz = -z[0] + z[1]
X, Z = np.meshgrid(x[nx1:nx2], z)
Global_moran = 0.e0
for itr, file in tqdm(enumerate(files)):
  filepath = os.path.join(u_dir, file)
  U = getScalar(filepath, Nx, Ny, Nz, 'u')
  u = U[:,ny,nx1:nx2]
  w = lat2W(Nz, nx2-nx1)
  global_moran = Moran(u.flatten(), w)
  Global_moran = (itr * Global_moran + global_moran.I) / (itr + 1.e0)

print(Global_moran)

