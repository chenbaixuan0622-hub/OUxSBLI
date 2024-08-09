import numpy as np
import os
import re
from pyinform.dist import Dist
from pyinform.utils import bin_series
from pyinform import transfer_entropy
from tqdm import tqdm
from mod.mod_read import getGrid, getMeanVector

data_directory = "../../../../../mnt/data1/TBL/KEEP20240802"
Q_directory    = data_directory

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

nx1 = int(0.5 * Nx)
nx2 = nx1 + 60

ny1 = 50

um_path = os.path.join(data_directory, "um.npy")
vm_path = os.path.join(data_directory, "vm.npy")
if os.path.isfile(um_path) & os.path.isfile(vm_path):
  um = np.load(um_path)
  vm = np.load(vm_path)
else:
  um, vm, _ = getMeanVector(Q_directory, Q_files, Nx, Ny, Nz, "velocity")
  np.save(um_path, um)
  np.save(vm_path, vm)

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

for i in tqdm(range(nx1+1,nx2)):
  t = 0
  u1[:] = 0.e0
  v1[:] = 0.e0
  u2[:] = 0.e0
  v2[:] = 0.e0
  for Q_file in Q_files:
    file_path = os.path.join(Q_directory, Q_file)
    u, v, _ = getScalar(file_path, Nx, Ny, Nz, 'velocity')

    for k in range(Nz):
      # point 1
      u1[t] += u[k,ny1,nx1] - um[k,ny1,nx1]
      v1[t] += v[k,ny1,nx1] - vm[k,ny1,nx1]
      # point 2
      u2[t] += u[k,ny1,i]   - um[k,ny1,i]
      v2[t] += v[k,ny1,i]   - vm[k,ny1,i]
      # longitudinal corr
      f[i-nx1] += u1[t] * u2[t] / np.sqrt(u1[t]**2 * u2[t]**2 + 1.e-20)
      # lateral corr
      g[i-nx1] += v1[t] * v2[t] / np.sqrt(v1[t]**2 * v2[t]**2 + 1.e-20)

    u1[t] = u1[t] / float(Nz)
    v1[t] = v1[t] / float(Nz)
    u2[t] = u2[t] / float(Nz)
    v2[t] = v2[t] / float(Nz)

    del u, v
    t += 1
  
  # calc longitudinal and lateral corr
  f[i-nx1] = f[i-nx1] / float(num_files * Nz)
  g[i-nx1] = g[i-nx1] / float(num_files * Nz)

  print("i=", i, "nx1=", nx1, "f=", f[i-nx1], "g=", g[i-nx1])

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

