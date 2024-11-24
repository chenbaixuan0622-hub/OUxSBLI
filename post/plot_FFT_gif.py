import numpy as np
import os
import re
from tqdm import tqdm
from mod.mod_FFT import calc_FFT_save_animation
from mod.mod_read import getGrid, getVector


def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_dir = ["../../a100/yp5", "../../a100/yp30", "../../a100/yp80", "../../a100/yp590"]
Lx    = 33.e-3
name  = "u_x=" + str(Lx)


Q_files = [f for f in os.listdir(Q_dir[0]) if f.endswith(".vtr")]
first_path = os.path.join(Q_dir[0], Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)
Nk    = len(Q_dir)
Nt    = len(Q_files)

data  = np.zeros((Nk,Nt,Nz))

for k in range(Nk):
  Q_files = [f for f in os.listdir(Q_dir[k]) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)
  nx = int(Lx / x[-1] * Nx)
  dz = -z[0] + z[1]

  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir[k], Q_file)
    if file_path == os.path.join(Q_dir[k], "TKE.vtr") or \
       file_path == os.path.join(Q_dir[k], "ReynoldsStress.vtr") or \
       file_path == os.path.join(Q_dir[k], "Qmean.vtr"):
      continue
    u, _, _ = getVector(file_path, Nx, Ny, Nz, 'velocity')
    data[k,itr,:] = u[:,1,nx]
    itr += 1

calc_FFT_save_animation(Nk, Nt, Nz, z, data, Q_dir[0], name)

