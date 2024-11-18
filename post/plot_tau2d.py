import os
import re
import numpy as np
import matplotlib.pyplot as plt
import mod.mod_plot as myplt
from mod.mod_read import getGrid, getScalar, getVector
from mod.mod_turb_stat import tau_2d

Q_directory = "../../../../../media/user/HD-EDS-E/hatayama/TBL/TBL20240819_KEEP6thVisc4th"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

Q_path              = os.path.join(Q_directory, Q_files[50])
Nx, Ny, Nz, x, y, z = getGrid(Q_path)

Q = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
Q[0,:,:,:]                         = getScalar(Q_path, Nx, Ny, Nz, 'rho'     )
Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:] = getVector(Q_path, Nx, Ny, Nz, 'velocity')
Q[4,:,:,:]                         = getScalar(Q_path, Nx, Ny, Nz, 'p'       )
# tw[nz,nx], ut[nz,nx]
tw, ut = tau_2d(Q, x, y, z)

X, Z = np.meshgrid(x, z)
plt.contourf(1.e3 * X, 1.e3 * Z, tw, cmap="turbo")
plt.colorbar()
plt.savefig("tw2d.png")
plt.close()

tw1d = np.mean(tw, axis=0)
plt.plot(1.e3 * x, tw1d)
plt.savefig("tw1d.png")
plt.close()

