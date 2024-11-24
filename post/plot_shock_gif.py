import numpy as np
import os
import re
from tqdm import tqdm
from mod.mod_FFT import calc_FFT_save_animation
from mod.mod_shock import reflected_shock
from mod.mod_read import getGrid, getVector


def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_dir = ["../../a100/yp590"]
Lx1   = 33.e-3
Lx2   = 35.e-3
name  = "shock"


Q_files = [f for f in os.listdir(Q_dir[0]) if f.endswith(".vtr")]
first_path = os.path.join(Q_dir[0], Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)
Nk    = len(Q_dir)
Nt    = len(Q_files)

data  = np.zeros((Nk,Nt,Nz))

for k in range(Nk):
  Q_files = [f for f in os.listdir(Q_dir[k]) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)
  nx1 = int(Lx1 / x[-1] * Nx)
  nx2 = int(Lx2 / x[-1] * Nx)
  dz = -z[0] + z[1]

  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir[k], Q_file)
    if file_path == os.path.join(Q_dir[k], "TKE.vtr") or \
       file_path == os.path.join(Q_dir[k], "ReynoldsStress.vtr") or \
       file_path == os.path.join(Q_dir[k], "Qmean.vtr"):
      continue
    u, v, w = getVector(file_path, Nx, Ny, Nz, 'velocity')
    ls, _ = reflected_shock(x[nx1:nx2], y, z, 1, u[:,:,nx1:nx2], v[:,:,nx1:nx2], w[:,:,nx1:nx2])
    ls = ls - np.mean(ls)
    data[k,itr,:] = ls*1e3
    itr += 1
  
calc_FFT_save_animation(Nk, Nt, Nz, z, data, Q_dir[0], name)

