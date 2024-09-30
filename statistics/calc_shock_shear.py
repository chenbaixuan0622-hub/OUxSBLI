import numpy as np
import os
import re
from tqdm import tqdm
from pyinform.utils import bin_series
from pyinform import transfer_entropy, conditional_entropy, block_entropy
from mod.mod_read import getGrid, getVector, getScalar

Q_directory = "../3D_solver/SSI/data/"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)

nx = 20
ny = 90
nz = int(Nz / 2)

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

u = np.zeros((Nz,Ny,Nx), dtype=np.float32)
v = np.zeros((Nz,Ny,Nx), dtype=np.float32)
w = np.zeros((Nz,Ny,Nx), dtype=np.float32)
p = np.zeros((Nz,Ny,Nx), dtype=np.float32)

ut = np.zeros((nx,num_files), dtype=np.float32)
vt = np.zeros((nx,num_files), dtype=np.float32)
wt = np.zeros((nx,num_files), dtype=np.float32)
pt = np.zeros((nx,num_files),    dtype=np.float32)
sl = np.zeros((num_files),    dtype=np.float32)

um = np.zeros((Nz,nx), dtype=np.float32)
vm = np.zeros((Nz,nx), dtype=np.float32)
wm = np.zeros((Nz,nx), dtype=np.float32)
pm = np.zeros((Nz,nx), dtype=np.float32)
u2 = np.zeros((Nz,nx), dtype=np.float32)
v2 = np.zeros((Nz,nx), dtype=np.float32)
w2 = np.zeros((Nz,nx), dtype=np.float32)
p2 = np.zeros((Nz,nx), dtype=np.float32)

Hu = np.zeros(nx, dtype=np.float32)
Hv = np.zeros(nx, dtype=np.float32)
Hw = np.zeros(nx, dtype=np.float32)
Hp = np.zeros(nx, dtype=np.float32)
TEu1 = np.zeros(nx, dtype=np.float32)
TEu2 = np.zeros(nx, dtype=np.float32)
TEv1 = np.zeros(nx, dtype=np.float32)
TEv2 = np.zeros(nx, dtype=np.float32)
TEw1 = np.zeros(nx, dtype=np.float32)
TEw2 = np.zeros(nx, dtype=np.float32)

Ruu = np.zeros(nx, dtype=np.float32)
Rvv = np.zeros(nx, dtype=np.float32)
Rww = np.zeros(nx, dtype=np.float32)
Rup = np.zeros(nx, dtype=np.float32)
Rvp = np.zeros(nx, dtype=np.float32)
Rwp = np.zeros(nx, dtype=np.float32)

itr = 0
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_directory, Q_file)
  u, v, w = getVector(file_path, Nx, Ny, Nz, 'velocity')
  p       = getScalar(file_path, Nx, Ny, Nz, 'p')
  dp = -p[nz,ny,1:] + p[nz,ny,:-1]
  ns = np.argmax(dp[200:-10]) + 200
  for i in range(ns-10,ns+10):
    ut[i-ns+10,itr] = np.sum(u[nz-1:nz+1,ny-1:ny+1,i-1:i+1])
    vt[i-ns+10,itr] = np.sum(v[nz-1:nz+1,ny-1:ny+1,i-1:i+1])
    wt[i-ns+10,itr] = np.sum(w[nz-1:nz+1,ny-1:ny+1,i-1:i+1])
    pt[i-ns+10,itr] = np.sum(p[nz-1:nz+1,ny-1:ny+1,i-1:i+1])
    # correlation
    Ruu[i-ns+10] += np.mean(u[:,:,ns-10] * u[:,:,i], axis=(0,1))
    Rvv[i-ns+10] += np.mean(v[:,:,ns-10] * v[:,:,i], axis=(0,1))
    Rww[i-ns+10] += np.mean(w[:,:,ns-10] * w[:,:,i], axis=(0,1))
    Rup[i-ns+10] += np.mean(u[:,:,i] * p[:,:,ns], axis=(0,1))
    Rvp[i-ns+10] += np.mean(v[:,:,i] * p[:,:,ns], axis=(0,1))
    Rwp[i-ns+10] += np.mean(w[:,:,i] * p[:,:,ns], axis=(0,1))
  # rms
  um += u[:,ny,ns-10:ns+10]
  vm += v[:,ny,ns-10:ns+10]
  wm += w[:,ny,ns-10:ns+10]
  pm += p[:,ny,ns-10:ns+10]
  u2 += u[:,ny,ns-10:ns+10]**2
  v2 += v[:,ny,ns-10:ns+10]**2
  w2 += w[:,ny,ns-10:ns+10]**2
  p2 += p[:,ny,ns-10:ns+10]**2
  itr += 1

um /= np.float32(num_files)
vm /= np.float32(num_files)
wm /= np.float32(num_files)
pm /= np.float32(num_files)
u2 /= np.float32(num_files)
v2 /= np.float32(num_files)
w2 /= np.float32(num_files)
p2 /= np.float32(num_files)

Ruu /= (np.float32(num_files))
Rvv /= (np.float32(num_files))
Rww /= (np.float32(num_files))
Rup /= (np.float32(num_files))
Rvp /= (np.float32(num_files))
Rwp /= (np.float32(num_files))
Ruu /= (Ruu[0])
Rvv /= (Rvv[0])
Rww /= (Rww[0])
Rup /= np.max(np.abs(Rup))
Rvp /= np.max(np.abs(Rvp))
Rwp /= np.max(np.abs(Rwp))

urms = np.mean(np.sqrt(u2 - um**2), axis=0)
vrms = np.mean(np.sqrt(v2 - vm**2), axis=0)
wrms = np.mean(np.sqrt(w2 - wm**2), axis=0)
prms = np.mean(np.sqrt(p2 - pm**2), axis=0)

del u, v, w, p, dp, um, vm, wm, pm, u2, v2, w2, p2

def NTE(x, y, bin, history_len):
  # x and y are time-series data
  # TEx_y = H(Yt|Yt-1:t-L) - H(Yt|Yt-1:t-L, Xt-1:t-L)
  xs = bin_series(x, b=bin)
  ys = bin_series(y, b=bin)
  TEx_y = transfer_entropy(xs[0], ys[0], k=history_len)
  TEy_x = transfer_entropy(ys[0], xs[0], k=history_len)
  # H(Yt|Yt-1)
  Yt   = bin_series(y[1:],  b=bin)
  Yt_1 = bin_series(y[:-1], b=bin)
  Xt   = bin_series(x[1:],  b=bin)
  Xt_1 = bin_series(x[:-1], b=bin)
  Hy = conditional_entropy(Yt_1[0], Yt[0]) # H(Yt|Yt-1)
  Hx = conditional_entropy(Xt_1[0], Xt[0]) # H(Xt|Xt-1)
  # E[TExd_y] = H(Yt|Yt-1:t-L) - H(Yt|Yt-1:t-L, Xdt-1:t-L)
  Ex_y = np.zeros(100, dtype=np.float32)
  Ey_x = np.zeros(100, dtype=np.float32)
  for i in range(100):
    x_shuffled = np.random.permutation(x)
    y_shuffled = np.random.permutation(y)
    xd = bin_series(x_shuffled, b=bin)
    yd = bin_series(y_shuffled, b=bin)
    Ex_y[i] = transfer_entropy(xd[0], ys[0], k=history_len)
    Ey_x[i] = transfer_entropy(yd[0], xs[0], k=history_len)
  # calc normalized TE
  NTEx_y = (TEx_y - np.mean(Ex_y)) / Hy
  NTEy_x = (TEy_x - np.mean(Ey_x)) / Hx
  return np.max(np.array([NTEx_y, 0.e0])), np.max(np.array([NTEy_x, 0.e0]))

def PSD(uf, t):
  fft  = np.fft.fft(uf)
  freq = np.fft.fftfreq(len(t), d=-t[0]+t[1])
  Amp  = abs(fft / (0.5e0 * float(len(t))))
  ps   = Amp**2
  psd  = ps / float(len(t))
  return freq, psd

bin = 5
history_len = 10
dx = -x[0] + x[1]
'''
for i in tqdm(range(nx)):
  xs = bin_series(ut[i,:], b=bin)
  Hu[i] = block_entropy(xs[0], k=history_len)
  xs = bin_series(vt[i,:], b=bin)
  Hv[i] = block_entropy(xs[0], k=history_len)
  xs = bin_series(wt[i,:], b=bin)
  Hw[i] = block_entropy(xs[0], k=history_len)
  xs = bin_series(pt[i,:], b=bin)
  Hp[i] = block_entropy(xs[0], k=history_len)
  TEu1[i], TEu2[i] = NTE(ut[i,:], pt[9,:], bin, history_len)
  TEv1[i], TEv2[i] = NTE(vt[i,:], pt[9,:], bin, history_len)
  TEw1[i], TEw2[i] = NTE(wt[i,:], pt[9,:], bin, history_len)

save_path = os.path.join(Q_directory, "shock_shear_rms.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# x        urms      vrms      wrms      prms      Hu        Hv        Hw        Hp", file=f)
  for i in range(nx):
    print(f'{dx*np.float32(i-9):.3e}', f'{urms[i]:.3e}', f'{vrms[i]:.3e}', f'{wrms[i]:.3e}', f'{prms[i]:.3e}', \
    f'{Hu[i]:.3e}', f'{Hv[i]:.3e}', f'{Hw[i]:.3e}', f'{Hp[i]:.3e}', file=f)
'''
save_path = os.path.join(Q_directory, "correlation.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# x        Ruu       Rvv       Rww       Rup       Rvp       Rwp", file=f)
  for i in range(nx):
    print(f'{dx*np.float32(i-9):.3e}', f'{Ruu[i]:.3e}', f'{Rvv[i]:.3e}', f'{Rww[i]:.3e}', \
    f'{Rup[i]:.3e}', f'{Rvp[i]:.3e}', f'{Rwp[i]:.3e}', file=f)
'''
save_path = os.path.join(Q_directory, "TE_shock_shear_interaction.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# x        TEu1      TEu2      TEv1      TEv2      TEw1      TEw2", file=f)
  for i in range(nx):
    print(f'{dx*np.float32(i-9):.3e}', f'{TEu1[i]:.3e}', f'{TEu2[i]:.3e}', f'{TEv1[i]:.3e}', f'{TEv2[i]:.3e}', f'{TEw1[i]:.3e}', f'{TEw2[i]:.3e}', file=f)
'''
