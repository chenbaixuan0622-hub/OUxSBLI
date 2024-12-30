import numpy as np
from mod.mod_Gaussian import Gaussian
from mod.mod_read import getGrid

Q_path = "../../Q00020.vtr"
Nx, Ny, Nz, x, y, z = getGrid(Q_path)

delta = 25.e0

sigma = 3.e0

nx1 = int(10.e0 / delta * Nx)
nx2 = int(22.e0 / delta * Nx)
ny1 = 0
ny2 = 256
nz1 = 0
nz2 = Nz

Gaussian(nx1,nx2,ny1,ny2,nz1,nz2,Q_path,sigma)


