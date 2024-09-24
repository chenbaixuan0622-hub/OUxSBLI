import numpy as np
import os
import re
from pyinform.dist import Dist
from pyinform.utils import bin_series
from pyinform import transfer_entropy, conditional_entropy
import statsmodels.api as sm
from tqdm import tqdm
from mod.mod_read import getGrid, getScalar, getVector, getMeanVector
from mod.mod_turb_stat import non_dim_tbl

Q_directory = "../../../../../media/user/HD-EDS-E/hatayama/TBL/TBL20240819_KEEP6thVisc4th"

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

yp_path = os.path.join(Q_directory, "yp.npy")
yp      = np.load(yp_path)

for j in range(Ny):
  if yp[j] > 80.e0:
    ny = j
    break

for j in range(Ny):
  if yp[j] > 5.e0:
    ny1 = j
    break

print("yp = 80 at ", ny, "yp = 5 at ", ny1)

um_path = os.path.join(Q_directory, "um.npy")
vm_path = os.path.join(Q_directory, "vm.npy")
wm_path = os.path.join(Q_directory, "wm.npy")
um = np.load(um_path)
vm = np.load(vm_path)
wm = np.load(wm_path)

Q   = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
kp  = np.zeros((num_files,ny-1), dtype=np.float32)

# TE
TE_12 = np.zeros(ny-1, dtype=np.float32)
TE_21 = np.zeros(ny-1, dtype=np.float32)

t = 0
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_directory, Q_file)
  Q[0,:,:,:]                         = getScalar(file_path, Nx, Ny, Nz, 'rho')
  Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:] = getVector(file_path, Nx, Ny, Nz, 'velocity')
  Q[4,:,:,:]                         = getScalar(file_path, Nx, Ny, Nz, 'p')
  kp[t,:] = 0.5e0 * ((Q[1,int(Nz/2),1:ny,int(Nx/2)] - um[int(Nz/2),1:ny,int(Nx/2)])**2 \
                   + (Q[2,int(Nz/2),1:ny,int(Nx/2)] - vm[int(Nz/2),1:ny,int(Nx/2)])**2 \
                   + (Q[3,int(Nz/2),1:ny,int(Nx/2)] - wm[int(Nz/2),1:ny,int(Nx/2)])**2)
  _, _, ut, _, _ = non_dim_tbl(Q, x, y, z)
  t += 1

kp = kp / np.mean(ut)**2

del Q, um, vm, wm, ut

def NTE(x, y, history_len):
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


print(kp[:,0])

print(kp[:,-1])

print(type(kp))
print(np.shape(kp))

# calc lon lat TE
for j in range(ny-1):
  TE_12[j], TE_21[j] = NTE(kp[:,j], kp[:,ny1], 5)

file_path = os.path.join(info_directory, "TE_tke_y.d")
with open(file_path, "w", encoding="UTF-8") as fo:
  print('# yp    TE_12   TE_21', file=fo)
  for j in range(ny-1):
    print(f'{yp[j+1]:.3e}', f'{TE_12[j]:.3e}', f'{TE_21[j]:.3e}', file=fo)

