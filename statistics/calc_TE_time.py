import numpy as np
import os
import re
from pyinform.dist import Dist
from pyinform.utils import bin_series
from pyinform import transfer_entropy, conditional_entropy
import statsmodels.api as sm
from tqdm import tqdm
from mod.mod_read import getGrid, getScalar

Q_directory = "../3D_solver/KHI/data"

info_directory = os.path.join(Q_directory, "info")
os.makedirs(info_directory, exist_ok = True)

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)

t = np.linspace(0.e0, 5.e0, num_files)

nx1 = int(0.25*Nx)
nx2 = int(0.75*Nx)
ny1 = int(0.7*Ny)
ny2 = int(0.8*Ny)
nz  = int(0.5*Nz)

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

rho = np.zeros((4,num_files), dtype=np.float32)

itr = 0
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_directory, Q_file)
  r = getScalar(file_path, Nx, Ny, Nz, 'rho')
  rho[0,itr] = r[nz,ny1,nx1]
  rho[1,itr] = r[nz,ny2,nx1]
  rho[2,itr] = r[nz,ny1,nx2]
  rho[3,itr] = r[nz,ny2,nx2]
  itr += 1

del r

def NTE(x, y, history_len):
  # x and y are time-series data
  # TEx_y = H(Yt|Yt-1:t-L) - H(Yt|Yt-1:t-L, Xt-1:t-L)
  xs = bin_series(x, b=5)
  ys = bin_series(y, b=5)
  TEx_y = transfer_entropy(xs[0], ys[0], k=history_len)
  TEy_x = transfer_entropy(ys[0], xs[0], k=history_len)
  # H(Yt|Yt-1)
  Yt   = bin_series(y[1:],  b=5)
  Yt_1 = bin_series(y[:-1], b=5)
  Xt   = bin_series(x[1:],  b=5)
  Xt_1 = bin_series(x[:-1], b=5)
  Hy = conditional_entropy(Yt_1[0], Yt[0]) # H(Yt|Yt-1)
  Hx = conditional_entropy(Xt_1[0], Xt[0]) # H(Xt|Xt-1)
  # E[TExd_y] = H(Yt|Yt-1:t-L) - H(Yt|Yt-1:t-L, Xdt-1:t-L)
  x_shuffled = np.random.permutation(x)
  y_shuffled = np.random.permutation(y)
  xd = bin_series(x_shuffled, b=5)
  yd = bin_series(y_shuffled, b=5)
  Ex_y = transfer_entropy(xd[0], ys[0], k=history_len)
  Ey_x = transfer_entropy(yd[0], xs[0], k=history_len)
  # calc normalized TE
  NTEx_y = TEx_y#(TEx_y - Ex_y)# / Hy
  NTEy_x = TEy_x#(TEy_x - Ey_x)# / Hx
  return NTEx_y, NTEy_x

# TE
batch = 20
TE1 = np.zeros(num_files-batch, dtype=np.float32)
TE2 = np.zeros(num_files-batch, dtype=np.float32)
TE3 = np.zeros(num_files-batch, dtype=np.float32)
TE4 = np.zeros(num_files-batch, dtype=np.float32)

# calc lon lat TE
count = 0
for itr in range(0, num_files-batch):
  TE1[count], TE2[count] = NTE(rho[0,itr:itr+batch], rho[1,itr:itr+batch], 10)
  TE3[count], TE4[count] = NTE(rho[2,itr:itr+batch], rho[3,itr:itr+batch], 10)
  print(itr, rho[0,itr], rho[1,itr], TE1[count], TE2[count], TE3[count], TE4[count])
  count += 1

file_path = os.path.join(info_directory, "TE_rho.d")
with open(file_path, "w", encoding="UTF-8") as fo:
  print('# t     rho', file=fo)
  for itr in range(len(TE1)):
    print(f'{0.5e0*(t[itr]+t[itr+batch]):.3e}', f'{TE1[itr]:.3e}', f'{TE2[itr]:.3e}', f'{TE3[itr]:.3e}', f'{TE4[itr]:.3e}', file=fo)

