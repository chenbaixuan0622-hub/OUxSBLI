import numpy as np
import os
import re
from numba import jit
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl

Q_directory = "../../../../../media/user/HD-EDS-E/hatayama/TBL20240813"
target_path = os.path.join(Q_directory, "Q01000.vtr")
save_path   = os.path.join(Q_directory, "yplus.d")

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

Q = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)

for Q_file in Q_files:
  file_path = os.path.join(Q_directory, Q_file)
  if file_path == target_path:
    Q[0,:,:,:]                          = getScalar(file_path, Nx, Ny, Nz, 'rho')
    Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:]  = getVector(file_path, Nx, Ny, Nz, 'velocity')
    Q[4,:,:,:]                          = getScalar(file_path, Nx, Ny, Nz, 'p')
    yp, tw, Qp, Qvd = non_dim_tbl(Q, x, y, z)
    print("tw is", np.mean(tw))
    with open(save_path, "w", encoding="UTF-8") as f:
      print("# yplus  uplus", file=f)
      for j in range(Ny):
        print(f'{yp[j]:.3e}', f'{np.mean(Qvd[0,:,j,:]):.3e}', file=f)

