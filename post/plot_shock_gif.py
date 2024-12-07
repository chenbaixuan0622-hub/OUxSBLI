import numpy as np
import os
import re
import matplotlib.pyplot as plt
from tqdm import tqdm
from mod.mod_FFT import calc_FFT_save_animation
from mod.mod_shock import reflected_shock
from mod.mod_read import getGrid, getScalar, getVector


def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_dir = ["../../a100/yp590"]
Lx1   = 33.e-3
Lx2   = 35.e-3
Lz    = [4.5e-3, 6.e-3]
endT  = 0.1e-3
name  = "shock"


Q_files = [f for f in os.listdir(Q_dir[0]) if f.endswith(".vtr")]
first_path = os.path.join(Q_dir[0], Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)
Nk    = len(Q_dir)
Nkz   = len(Lz)
Nt    = len(Q_files)

data  = np.zeros((Nk,Nt,Nz))

# shock strength
rho01 = np.zeros(Nkz)
u01   = np.zeros(Nkz)
p01   = np.zeros(Nkz)

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
    Rho     = getScalar(file_path, Nx, Ny, Nz, 'rho')
    u, v, w = getVector(file_path, Nx, Ny, Nz, 'velocity')
    P       = getScalar(file_path, Nx, Ny, Nz, 'p')
    ls, _, ix, _ = reflected_shock(x[nx1:nx2], y, z, 1, u[:,:,nx1:nx2], v[:,:,nx1:nx2], w[:,:,nx1:nx2])
    data[k,itr,:] = ls*1e3
    # calc shock strength
    rho = Rho[:,:,nx1:nx2]
    p   =   P[:,:,nx1:nx2]
    for kz in range(Nkz):
      nz = int(Lz[kz]  / z[-1] * Nz)
      rho01[kz] += rho[nz,0,ix[nz]] / rho[nz,0,ix[nz]-10]
      p01[kz]   +=   p[nz,0,ix[nz]] /   p[nz,0,ix[nz]-10]
    itr += 1

calc_FFT_save_animation(Nk, Nt, Nz, z, data, Q_dir[0], name)

rho01 /= float(itr)
p01   /= float(itr)

Nt = itr
t  = np.linspace(0.e0, endT, Nt)

for kz in range(Nkz):
  # print shock strength
  name = "shock_strength_" + "z=" + str(Lz[kz])
  save_path = os.path.join(Q_dir[0], name + ".d")
  with open(save_path, "w", encoding="UTF-8") as f:
    print("# rho01   p01", file=f)
    print(f'{rho01[kz]:.3e}', f'{p01[kz]:.3e}', file=f)

  # print time series data
  nz   = int(Lz[kz]  / z[-1] * Nz)
  name = "shock_" + "z=" + str(Lz[kz])
  save_path = os.path.join(Q_dir[0], name + ".d")
  with open(save_path, "w", encoding="UTF-8") as f:
    print("# t[ms]   x[mm]", file=f)
    for i in range(Nt):
      print(f'{t[i]*1e3:.3e}', f'{data[0,i,nz]:.3e}', file=f)

mean = np.mean(data[0,:,:], axis=0)
name = "shock_" + "mean"
save_path = os.path.join(Q_dir[0], name + ".d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# z[mm]   x[mm]", file=f)
  for k in range(Nz):
    print(f'{z[k]*1e3:.3e}', f'{mean[k]:.3e}', file=f)

