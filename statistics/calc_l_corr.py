import numpy as np
from numpy.fft import fft
import os
import re
from pyinform.dist import Dist
from pyinform.utils import bin_series
from pyinform.shannon import relative_entropy
from pyinform.shannon import mutual_info
from pyinform import transfer_entropy
from mod.mod_read import gridInfo

data_directory = "./np_data"
Qp_directory   = os.path.join(data_directory, "Qp")
file_name      = "Qm.npy"
file_path      = os.path.join(Qp_directory, file_name)
Qm             = np.load(file_path)

x, y, z, Nx, Ny, Nz = gridInfo(data_directory)

# get y+
file_path = os.path.join(data_directory, "yplus.npy")
yp = np.load(file_path)
file_path = os.path.join(data_directory, "xplus.npy")
xp = np.load(file_path)

nx1 = int(0.5 * Nx)
nx2 = nx1 + 30

yp1 = 40
for j in range(len(yp)):
  if yp[j] <= yp1:
    ny1 = j
  else:
    break

Qp_files = [f for f in os.listdir(Qp_directory) if f.endswith(".npy")]
num_files = len(Qp_files)

def extract_number(filename):
  match = re.search(r'Q(\d+)\.npy$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Qp_files.sort(key=extract_number)

# fluctuating velocity
u1 = np.zeros((num_files), dtype=np.float32)
v1 = np.zeros((num_files), dtype=np.float32)
u2 = np.zeros((num_files), dtype=np.float32)
v2 = np.zeros((num_files), dtype=np.float32)

# longitudinal / lateral correlation
f = np.ones((nx2-nx1), dtype=np.float32)
g = np.ones((nx2-nx1), dtype=np.float32)
# TE
TE_lon12 = np.zeros((nx2-nx1), dtype=np.float32)
TE_lon21 = np.zeros((nx2-nx1), dtype=np.float32)
TE_lat12 = np.zeros((nx2-nx1), dtype=np.float32)
TE_lat21 = np.zeros((nx2-nx1), dtype=np.float32)

for i in range(nx1+1,nx2):
  t = 0
  u1[:] = 0.e0
  v1[:] = 0.e0
  u2[:] = 0.e0
  v2[:] = 0.e0
  for Qp_file in Qp_files:
    file_path = os.path.join(Qp_directory, Qp_file)
    Qp = np.load(file_path)

    for k in range(Nz):
      # point 1
      u1[t] += Qp[1,k,ny1,nx1] - Qm[1,k,ny1,nx1]
      v1[t] += Qp[2,k,ny1,nx1] - Qm[2,k,ny1,nx1]
      # point 2
      u2[t] += Qp[1,k,ny1,i]   - Qm[1,k,ny1,i]
      v2[t] += Qp[2,k,ny1,i]   - Qm[2,k,ny1,i]
      # longitudinal corr
      f[i-nx1] += u1[t] * u2[t] / np.sqrt(u1[t]**2 * u2[t]**2 + 1.e-20)
      # lateral corr
      g[i-nx1] += v1[t] * v2[t] / np.sqrt(v1[t]**2 * v2[t]**2 + 1.e-20)

    u1[t] = u1[t] / float(Nz)
    v1[t] = v1[t] / float(Nz)
    u2[t] = u2[t] / float(Nz)
    v2[t] = v2[t] / float(Nz)

    t += 1
  
  # calc longitudinal and lateral corr
  f[i-nx1] = f[i-nx1] / float(num_files * Nz)
  g[i-nx1] = g[i-nx1] / float(num_files * Nz)

  print("i=", i, "nx1=", nx1, "f=", f[i-nx1], "g=", g[i-nx1])

  # calc TE
  xs = bin_series(u1[:], b=10)
  ys = bin_series(u2[:], b=10)
  TE_lon12[i-nx1] = transfer_entropy(xs[0], ys[0], k=2)
  TE_lon21[i-nx1] = transfer_entropy(ys[0], xs[0], k=2)
  xs = bin_series(v1[:], b=10)
  ys = bin_series(v2[:], b=10)
  TE_lat12[i-nx1] = transfer_entropy(xs[0], ys[0], k=2)
  TE_lat21[i-nx1] = transfer_entropy(ys[0], xs[0], k=2)


# print output
info_directory = os.path.join(data_directory, "info")
os.makedirs(info_directory, exist_ok = True)

# print longitudianl / lateral correlation
file_path = os.path.join(info_directory, "vel_corr.d")
with open(file_path, "w", encoding="UTF-8") as fo:
  for i in range(nx2-nx1):
    print(f'{-xp[nx1]+xp[nx1+i+1]:.3e}', f'{f[i]:.3e}', f'{g[i]:.3e}', file=fo)

# print TE
file_path = os.path.join(info_directory, "TE_lon_lat.d")
with open(file_path, "w", encoding="UTF-8") as fo:
  for i in range(nx2-nx1):
    print(f'{-xp[nx1]+xp[nx1+i+1]:.3e}', f'{TE_lon12[i]:.3e}', f'{TE_lon21[i]:.3e}', f'{TE_lat12[i]:.3e}', f'{TE_lat21[i]:.3e}', file=fo)

