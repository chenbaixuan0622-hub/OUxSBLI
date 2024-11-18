import numpy as np
import matplotlib.pyplot as plt
from mod.mod_read import getGrid, getVector, getScalar


Q_path = "Q00020.vtr"
Nx, Ny, Nz, X, Y, Z = getGrid(Q_path)
u, v, w = getVector(file_path, Nx, Ny, Nz, 'velocity')

dz = -Z[0] + Z[1]
lz = Lz / 99
u = filter_low_freq(U, dz, lz)



