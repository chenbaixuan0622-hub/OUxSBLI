import numpy as np
from numpy.fft import fft
import os
import re
from pyinform.dist import Dist
from pyinform.utils import bin_series
from pyinform import transfer_entropy
from mod.mod_read import getGrid

data_directory = "./data"

# shock region
nsx1 = 
nsx2 = 
nsy1 = 
nsy2 = 

# tbl region
ntx1 = 
ntx2 = 
nty1 = 
nty2 =

Q_files = [f for f in os.listdir(data_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

file_path = os.path.join(data_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(file_path)

# TE
TE_lon12 = np.zeros((nx2-nx1), dtype=np.float32)

for i in range(nx1+1,nx2):
  t = 0
  for Q_file in Q_files:
    file_path = os.path.join(data_directory, Q_file)
    rho       = getScalar(file_path, Nx, Ny, Nz, 'rho')
    p         = getScalar(file_path, Nx, Ny, Nz, 'p')
    u, v, w   = getVector(file_path, Nx, Ny, Nz, 'velocity')

    for k in range(Nz):
      # point 1
      u1[t] += Qp[1,k,ny1,nx1] - Qm[1,k,ny1,nx1]
      v1[t] += Qp[2,k,ny1,nx1] - Qm[2,k,ny1,nx1]
      # point 2
      u2[t] += Qp[1,k,ny1,i]   - Qm[1,k,ny1,i]
      v2[t] += Qp[2,k,ny1,i]   - Qm[2,k,ny1,i]

    u1[t] = u1[t] / float(Nz)
    v1[t] = v1[t] / float(Nz)
    u2[t] = u2[t] / float(Nz)
    v2[t] = v2[t] / float(Nz)

    t += 1
  
  # calc TE
  xs = bin_series(u1[:], b=10)
  ys = bin_series(u2[:], b=10)
  TE_lon12[i-nx1] = transfer_entropy(xs[0], ys[0], k=10)
  TE_lon21[i-nx1] = transfer_entropy(ys[0], xs[0], k=10)
  xs = bin_series(v1[:], b=10)
  ys = bin_series(v2[:], b=10)
  TE_lat12[i-nx1] = transfer_entropy(xs[0], ys[0], k=10)
  TE_lat21[i-nx1] = transfer_entropy(ys[0], xs[0], k=10)


# print output
info_directory = os.path.join(data_directory, "info")
os.makedirs(info_directory, exist_ok = True)

# print TE
file_path = os.path.join(info_directory, "TE_lon_lat.d")
with open(file_path, "w", encoding="UTF-8") as fo:
  for i in range(nx2-nx1):
    print(f'{-xp[nx1]+xp[nx1+i+1]:.3e}', f'{TE_lon12[i]:.3e}', f'{TE_lon21[i]:.3e}', f'{TE_lat12[i]:.3e}', f'{TE_lat21[i]:.3e}', file=fo)

