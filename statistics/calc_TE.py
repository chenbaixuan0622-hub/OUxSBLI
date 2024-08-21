import numpy as np
import os
import re
from pyinform.dist import Dist
from pyinform.utils import bin_series
from pyinform import transfer_entropy
import statsmodels.api as sm
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getMeanVector

Q_directory = "../../../../../mnt/data1/TBL20240819"

info_directory = os.path.join(Q_directory, "info")
os.makedirs(info_directory, exist_ok = True)

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

# set parameters
nx1 = int(0.5 * Nx)
nx2 = nx1 + 60
ny1 = 48  # log layer
ny0 = 8   # visc sub layer
ut  = 22.e0


um_path = os.path.join(Q_directory, "um.npy")
vm_path = os.path.join(Q_directory, "vm.npy")
wm_path = os.path.join(Q_directory, "wm.npy")
if os.path.isfile(um_path) & os.path.isfile(vm_path):
  um = np.load(um_path)
  vm = np.load(vm_path)
  wm = np.load(wm_path)
else:
  um, vm, wm = getMeanVector(Q_directory, Q_files, Nx, Ny, Nz, "velocity")
  np.save(um_path, um)
  np.save(vm_path, vm)
  np.save(wm_path, wm)


# TE
TE_lon12 = np.zeros((nx2-nx1), dtype=np.float32)
TE_lon21 = np.zeros((nx2-nx1), dtype=np.float32)
TE_lat12 = np.zeros((nx2-nx1), dtype=np.float32)
TE_lat21 = np.zeros((nx2-nx1), dtype=np.float32)

TE_vis12 = np.zeros((nx2-nx1), dtype=np.float32)
TE_vis21 = np.zeros((nx2-nx1), dtype=np.float32)
TE_log12 = np.zeros((nx2-nx1), dtype=np.float32)
TE_log21 = np.zeros((nx2-nx1), dtype=np.float32)

uf_path = os.path.join(Q_directory, "uf.npy")
vf_path = os.path.join(Q_directory, "vf.npy")
tke_vis_path = os.path.join(Q_directory, "tke_vis.npy")
tke_log_path = os.path.join(Q_directory, "tke_log.npy")


if os.path.isfile(uf_path) & os.path.isfile(vf_path) & os.path.isfile(tke_vis_path) & os.path.isfile(tke_log_path):
  uf = np.load(uf_path)
  vf = np.load(vf_path)
  tke_vis = np.load(tke_vis_path)
  tke_log = np.load(tke_log_path)
else:
  # fluctuating velocity
  uf = np.zeros((nx2-nx1, num_files), dtype=np.float32)
  vf = np.zeros((nx2-nx1, num_files), dtype=np.float32)
  # turbulent kinetic energy
  tke_vis = np.zeros((nx2-nx1, num_files), dtype=np.float32)
  tke_log = np.zeros((nx2-nx1, num_files), dtype=np.float32)
  # longitudinal / lateral correlation
  f = np.ones((nx2-nx1), dtype=np.float32)
  g = np.ones((nx2-nx1), dtype=np.float32)
  t = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_directory, Q_file)
    u, v, w   = getVector(file_path, Nx, Ny, Nz, 'velocity')
    Uf = u - um
    Vf = v - vm
    Wf = w - wm
    for i in range(nx1,nx2):
      uf[i-nx1,t] += np.mean(Uf[:,ny1,i])
      vf[i-nx1,t] += np.mean(Vf[:,ny1,i])
      f[i-nx1] += uf[0,t] * uf[i-nx1,t] / np.sqrt(uf[0,t]**2 * uf[i-nx1,t]**2 + 1.e-20)
      g[i-nx1] += vf[0,t] * vf[i-nx1,t] / np.sqrt(vf[0,t]**2 * vf[i-nx1,t]**2 + 1.e-20)
      # K+
      tke_vis[i-nx1,t] = 0.5e0 * np.mean(Uf[:,ny0,i]**2 + Vf[:,ny0,i]**2 + Wf[:,ny0,i]**2) / (ut**2)
      tke_log[i-nx1,t] = 0.5e0 * np.mean(Uf[:,ny1,i]**2 + Vf[:,ny1,i]**2 + Wf[:,ny1,i]**2) / (ut**2)
    t += 1
  np.save(uf_path, uf)
  np.save(vf_path, vf)
  np.save(tke_vis_path, tke_vis)
  np.save(tke_log_path, tke_log)

  f = f / float(num_files)
  g = g / float(num_files)

  file_path = os.path.join(info_directory, "lon_lat_corr.d")
  with open(file_path, "w", encoding="UTF-8") as fo:
    for i in range(nx2-nx1):
      print(f'{-x[nx1]+x[nx1+i]:.3e}', f'{f[i]:.3e}', f'{g[i]:.3e}', file=fo)

  del u, v, f, g

del um, vm


# calc lon lat TE
for i in range(nx2-nx1):
  xs = bin_series(uf[0,:], b=10)
  ys = bin_series(uf[i,:], b=10)
  TE_lon12[i] = transfer_entropy(xs[0], ys[0], k=5)
  TE_lon21[i] = transfer_entropy(ys[0], xs[0], k=5)
  xs = bin_series(vf[0,:], b=10)
  ys = bin_series(vf[i,:], b=10)
  TE_lat12[i] = transfer_entropy(xs[0], ys[0], k=5)
  TE_lat21[i] = transfer_entropy(ys[0], xs[0], k=5)

file_path = os.path.join(info_directory, "TE_lon_lat.d")
with open(file_path, "w", encoding="UTF-8") as fo:
  for i in range(nx2-nx1):
    print(f'{-x[nx1]+x[nx1+i]:.3e}', f'{TE_lon12[i]:.3e}', f'{TE_lon21[i]:.3e}', f'{TE_lat12[i]:.3e}', f'{TE_lat21[i]:.3e}', file=fo)

del uf, vf

# calc tke TE
for i in range(nx2-nx1):
  xs = bin_series(tke_vis[0,:], b=10)
  ys = bin_series(tke_vis[i,:], b=10)
  TE_vis12[i] = transfer_entropy(xs[0], ys[0], k=5)
  TE_vis21[i] = transfer_entropy(ys[0], xs[0], k=5)
  xs = bin_series(tke_log[0,:], b=10)
  ys = bin_series(tke_log[i,:], b=10)
  TE_log12[i] = transfer_entropy(xs[0], ys[0], k=5)
  TE_log21[i] = transfer_entropy(ys[0], xs[0], k=5)

file_path = os.path.join(info_directory, "TE_tke.d")
with open(file_path, "w", encoding="UTF-8") as fo:
  for i in range(nx2-nx1):
    print(f'{-x[nx1]+x[nx1+i]:.3e}', f'{TE_vis12[i]:.3e}', f'{TE_vis21[i]:.3e}', f'{TE_log12[i]:.3e}', f'{TE_log21[i]:.3e}', file=fo)


# calc auto correlation
u_acf = sm.tsa.stattools.acf(uf, nlags=num_files-1)
v_acf = sm.tsa.stattools.acf(vf, nlags=num_files-1)
time  = np.linspace(0.e0, 4.e-3, num_files)

file_path = os.path.join(info_directory, "uv_acf.d")
with open(file_path, "w", encoding="UTF-8") as fo:
  for i in range(num_files):
    print(f'{time[i]:.3e}', f'{u_acf[i]:.3e}', f'{v_acf[i]:.3e}', file=fo)

