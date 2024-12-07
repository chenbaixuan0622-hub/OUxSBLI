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

Q_dir = ["../../a100/yp5", "../../a100/yp30", "../../a100/yp80", "../../a100/yp269", "../../a100/yp590"]
Lx    = 33.e-3
Lz    = [4.5e-3, 6.e-3]
endT  = 0.1e-3
name  = "u_x=" + str(Lx)


Q_files = [f for f in os.listdir(Q_dir[0]) if f.endswith(".vtr")]
first_path = os.path.join(Q_dir[0], Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)
Nk    = len(Q_dir)
Nkz   = len(Lz)
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


Nt = itr
t  = np.linspace(0.e0, endT, Nt)

for kz in range(Nkz):
  nz   = int(Lz[kz]  / z[-1] * Nz)
  name = "u_" + "x=" + str(Lx) + "_z=" + str(Lz[kz])
  save_path = os.path.join(Q_dir[0], name + ".d")
  with open(save_path, "w", encoding="UTF-8") as f:
    print("# t[ms]   u[m/s]", file=f)
    for i in range(Nt):
      print(f'{t[i]*1e3:.3e}', f'{data[0,i,nz]:.3e}', f'{data[1,i,nz]:.3e}', \
            f'{data[2,i,nz]:.3e}', f'{data[3,i,nz]:.3e}', f'{data[4,i,nz]:.3e}', file=f)

