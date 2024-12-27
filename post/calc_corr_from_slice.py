import numpy as np
import os
import re
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, extract_number


#Q_dir = ["../../a100/yp5", "../../a100/yp30", "../../a100/yp80", "../../a100/yp269", "../../a100/yp590"]
Q_dir = ["../3D_solver/TBL/data"]
Lxs   = [0.e-3, 20.e-3]
Lzs   = [0.e-3, 4.e-3]
Ly    = 0.56e0 * 2.e-3

Q_files = [f for f in os.listdir(Q_dir[0]) if f.endswith(".vtr")]
first_path = os.path.join(Q_dir[0], Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)
Nk = len(Q_dir)

nx1 = int(Lxs[0] / x[-1] * Nx)
nx2 = int(Lxs[1] / x[-1] * Nx)
nz1 = int(Lzs[0] / z[-1] * Nz)
nz2 = int(Lzs[1] / z[-1] * Nz)

for j in range(Ny):
  if y[j] > Ly:
    ny = j
    break

um  = np.zeros(nx2-nx1)
vm  = np.zeros(nx2-nx1)
wm  = np.zeros(nx2-nx1)

Ruu = np.zeros((Nk,nx2-nx1))
Rvv = np.zeros((Nk,nx2-nx1))
Rww = np.zeros((Nk,nx2-nx1))

for k in range(Nk):
  Q_files = [f for f in os.listdir(Q_dir[k]) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)

  # calc average
  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir[k], Q_file)
    U, V, W = getVector(file_path, Nx, Ny, Nz, 'velocity')
    um += np.mean(U[nz1:nz2,ny,nx1:nx2], axis=0)
    vm += np.mean(V[nz1:nz2,ny,nx1:nx2], axis=0)
    wm += np.mean(W[nz1:nz2,ny,nx1:nx2], axis=0)
    itr += 1
  um /= float(itr)
  vm /= float(itr)
  wm /= float(itr)

  # calc corr
  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir[k], Q_file)
    U, V, W = getVector(file_path, Nx, Ny, Nz, 'velocity')
    u = np.mean(U[nz1:nz2,ny,nx1:nx2], axis=0)
    v = np.mean(V[nz1:nz2,ny,nx1:nx2], axis=0)
    w = np.mean(W[nz1:nz2,ny,nx1:nx2], axis=0)
    for i in range(nx2-nx1):
      Ruu[k,i] += (u[0] - um[0]) * (u[i] - um[i])
      Rvv[k,i] += (v[0] - vm[0]) * (v[i] - vm[i])
      Rww[k,i] += (w[0] - wm[0]) * (w[i] - wm[i])
    itr += 1
  Ruu[k,:] /= float(itr)
  Rvv[k,:] /= float(itr)
  Rww[k,:] /= float(itr)
  Ruu[k,:] = Ruu[k,:] / Ruu[k,0]
  Rvv[k,:] = Rvv[k,:] / Rvv[k,0]
  Rww[k,:] = Rww[k,:] / Rww[k,0]

  name = "corr_" + "x" + f'{Lxs[0]*1e3:.2f}' + "_" + f'{Lxs[1]*1e3:.2f}' \
                + "_z" + f'{Lzs[0]*1e3:.2f}' + "_" + f'{Lzs[1]*1e3:.2f}'
  save_path = os.path.join(Q_dir[k], name + ".d")
  with open(save_path, "w", encoding="UTF-8") as f:
    print("# x[mm]   Ruu       Rvv       Rww", file=f)
    for i in range(nx2-nx1):
      print(f'{x[i]*1e3:.3e}', f'{Ruu[k,i]:.3e}', \
            f'{Rvv[k,i]:.3e}', f'{Rww[k,i]:.3e}', file=f)

