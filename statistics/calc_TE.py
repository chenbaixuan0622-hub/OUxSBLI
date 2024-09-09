import numpy as np
import os
import re
from pyinform.dist import Dist
from pyinform.utils import bin_series
from pyinform import transfer_entropy, conditional_entropy
import statsmodels.api as sm
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar, getMeanVector

Q_directory = "../../../../../media/user/HD-EDS-E/hatayama/TBL/TBL20240826_KEEP4thVisc2nd"

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
ny1 = 49  # log layer
ny0 = 8   # visc sub layer
ut  = 22.e0


um_path   = os.path.join(Q_directory, "um.npy"  )
vm_path   = os.path.join(Q_directory, "vm.npy"  )
wm_path   = os.path.join(Q_directory, "wm.npy"  )
rhom_path = os.path.join(Q_directory, "rhom.npy")
pm_path   = os.path.join(Q_directory, "pm.npy"  )
if os.path.isfile(um_path) & os.path.isfile(vm_path):
  um   = np.load(um_path)
  vm   = np.load(vm_path)
  wm   = np.load(wm_path)
  rhom = np.load(rhom_path)
  pm   = np.load(pm_path)

uf_path   = os.path.join(Q_directory, "uf.npy"  )
vf_path   = os.path.join(Q_directory, "vf.npy"  )
wf_path   = os.path.join(Q_directory, "wf.npy"  )
rhof_path = os.path.join(Q_directory, "rhof.npy")
pf_path   = os.path.join(Q_directory, "pf.npy"  )
tke_vis_path = os.path.join(Q_directory, "tke_vis.npy")
tke_log_path = os.path.join(Q_directory, "tke_log.npy")


if os.path.isfile(uf_path) & os.path.isfile(vf_path) & os.path.isfile(tke_vis_path) & os.path.isfile(tke_log_path) &\
   os.path.isfile(rhof_path) & os.path.isfile(pf_path):
  uf   = np.load(uf_path  )
  vf   = np.load(vf_path  )
  wf   = np.load(wf_path  )
  rhof = np.load(rhof_path)
  pf   = np.load(pf_path  )
  tke_vis = np.load(tke_vis_path)
  tke_log = np.load(tke_log_path)
else:
  # fluctuating quantities
  uf   = np.zeros((nx2-nx1, num_files), dtype=np.float32)
  vf   = np.zeros((nx2-nx1, num_files), dtype=np.float32)
  wf   = np.zeros((nx2-nx1, num_files), dtype=np.float32)
  rhof = np.zeros((nx2-nx1, num_files), dtype=np.float32)
  pf   = np.zeros((nx2-nx1, num_files), dtype=np.float32)
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
    rho       = getScalar(file_path, Nx, Ny, Nz, 'rho'     )
    p         = getScalar(file_path, Nx, Ny, Nz, 'p'       )
    Uf   = u   - um
    Vf   = v   - vm
    Wf   = w   - wm
    Rhof = rho - rhom
    Pf   = p   - pm
    for i in range(nx1,nx2):
      uf[i-nx1,t]   += np.mean(Uf[:,ny1,i])
      vf[i-nx1,t]   += np.mean(Vf[:,ny1,i])
      wf[i-nx1,t]   += np.mean(Wf[:,ny1,i])
      rhof[i-nx1,t] += np.mean(Rhof[:,ny1,i])
      pf[i-nx1,t]   += np.mean(Pf[:,ny1,i])
      f[i-nx1] += uf[0,t] * uf[i-nx1,t] / np.sqrt(uf[0,t]**2 * uf[i-nx1,t]**2 + 1.e-20)
      g[i-nx1] += vf[0,t] * vf[i-nx1,t] / np.sqrt(vf[0,t]**2 * vf[i-nx1,t]**2 + 1.e-20)
      # K+
      tke_vis[i-nx1,t] = 0.5e0 * np.mean(Uf[:,ny0,i]**2 + Vf[:,ny0,i]**2 + Wf[:,ny0,i]**2) / (ut**2)
      tke_log[i-nx1,t] = 0.5e0 * np.mean(Uf[:,ny1,i]**2 + Vf[:,ny1,i]**2 + Wf[:,ny1,i]**2) / (ut**2)
    t += 1
  np.save(uf_path,   uf)
  np.save(vf_path,   vf)
  np.save(wf_path,   wf)
  np.save(rhof_path, rhof)
  np.save(pf_path,   pf)
  np.save(tke_vis_path, tke_vis)
  np.save(tke_log_path, tke_log)

  f = f / float(num_files)
  g = g / float(num_files)

  file_path = os.path.join(info_directory, "lon_lat_corr.d")
  with open(file_path, "w", encoding="UTF-8") as fo:
    for i in range(nx2-nx1):
      print(f'{-x[nx1]+x[nx1+i]:.3e}', f'{f[i]:.3e}', f'{g[i]:.3e}', file=fo)

  del u, v, w, rho, p, f, g

del um, vm, wm, rhom, pm, Uf, Vf, Wf, Rhof, Pf


def NTE(x, y, power, history_len):
  # x and y are time-series data
  # TEx_y = H(Yt|Yt-1:t-L) - H(Yt|Yt-1:t-L, Xt-1:t-L)
  xs = bin_series(x, b=10)
  ys = bin_series(y, b=10)
  TEx_y = transfer_entropy(xs[0], ys[0], k=history_len)
  TEy_x = transfer_entropy(ys[0], xs[0], k=history_len)
  # H(Yt|Yt-1)
  Yt   = bin_series(y[1:],  b=10)
  Yt_1 = bin_series(y[:-1], b=10)
  Xt   = bin_series(x[1:],  b=10)
  Xt_1 = bin_series(x[:-1], b=10)
  Hy = conditional_entropy(Yt_1[0], Yt[0]) # H(Yt|Yt-1)
  Hx = conditional_entropy(Xt_1[0], Xt[0]) # H(Xt|Xt-1)
  # E[TExd_y] = H(Yt|Yt-1:t-L) - H(Yt|Yt-1:t-L, Xdt-1:t-L)
  x_shuffled = np.random.permutation(x)
  y_shuffled = np.random.permutation(y)
  xd = bin_series(x_shuffled, b=10)
  yd = bin_series(y_shuffled, b=10)
  Ex_y = transfer_entropy(xd[0], ys[0], k=history_len)
  Ey_x = transfer_entropy(yd[0], xs[0], k=history_len)
  # calc normalized TE
  NTEx_y = (TEx_y - Ex_y) / Hy
  NTEy_x = (TEy_x - Ey_x) / Hx
  return NTEx_y, NTEy_x

# TE
TE_lon12 = np.zeros((nx2-nx1), dtype=np.float32)
TE_lon21 = np.zeros((nx2-nx1), dtype=np.float32)
TE_lat12 = np.zeros((nx2-nx1), dtype=np.float32)
TE_lat21 = np.zeros((nx2-nx1), dtype=np.float32)

# calc lon lat TE
for i in range(nx2-nx1):
  TE_lon12[i], TE_lon21[i] = NTE(uf[0,:], uf[i,:], 0.5e0, 5)
  TE_lat12[i], TE_lat21[i] = NTE(vf[0,:], vf[i,:], 0.5e0, 5)


file_path = os.path.join(info_directory, "TE_lon_lat.d")
with open(file_path, "w", encoding="UTF-8") as fo:
  print('# x     TE_lon  TE_lon  TE_lat  TE_lat', file=fo)
  for i in range(nx2-nx1):
    print(f'{-x[nx1]+x[nx1+i]:.3e}', f'{TE_lon12[i]:.3e}', f'{TE_lon21[i]:.3e}', f'{TE_lat12[i]:.3e}', f'{TE_lat21[i]:.3e}', file=fo)

del TE_lon12, TE_lon21, TE_lat12, TE_lat21

# TE
TE_vis12 = np.zeros((nx2-nx1), dtype=np.float32)
TE_vis21 = np.zeros((nx2-nx1), dtype=np.float32)
TE_log12 = np.zeros((nx2-nx1), dtype=np.float32)
TE_log21 = np.zeros((nx2-nx1), dtype=np.float32)

# calc tke TE
for i in range(nx2-nx1):
  TE_vis12[i], TE_vis21[i] = NTE(tke_vis[0,:], tke_vis[i,:], 0.1e0, 5)
  TE_log12[i], TE_log21[i] = NTE(tke_log[0,:], tke_log[i,:], 0.1e0, 5)

file_path = os.path.join(info_directory, "TE_tke.d")
with open(file_path, "w", encoding="UTF-8") as fo:
  print('# x     TE_vis  TE_vis  TE_lon  TE_lon', file=fo)
  for i in range(nx2-nx1):
    print(f'{-x[nx1]+x[nx1+i]:.3e}', f'{TE_vis12[i]:.3e}', f'{TE_vis21[i]:.3e}', f'{TE_log12[i]:.3e}', f'{TE_log21[i]:.3e}', file=fo)

del TE_vis12, TE_vis21, TE_log12, TE_log21

'''
# calc auto correlation
u_acf = sm.tsa.stattools.acf(uf[0,:], nlags=num_files-1)
v_acf = sm.tsa.stattools.acf(vf[0,:], nlags=num_files-1)
time  = np.linspace(0.e0, 4.e-3, num_files)

file_path = os.path.join(info_directory, "uv_acf.d")
with open(file_path, "w", encoding="UTF-8") as fo:
  print('# time  u_acf  v_acf', file=fo)
  for i in range(num_files):
    print(f'{time[i]:.3e}', f'{u_acf[i]:.3e}', f'{v_acf[i]:.3e}', file=fo)
'''

