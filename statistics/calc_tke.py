import numpy as np
import os
import re
from scipy.fftpack import fft, ifft, fftfreq
from scipy import signal
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl, Sutherland

Q_directory = "../../../../../mnt/data1/TBL_5thHRSLAU2_20240922"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)

# parameters
endT = 2.087e-3
nyp  = 48

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

def PSD(uf, t):
  fft  = np.fft.fft(uf)
  freq = np.fft.fftfreq(len(t), d=-t[0]+t[1])
  Amp  = abs(fft / (0.5e0 * float(len(t))))
  ps   = Amp**2
  psd  = ps / float(len(t))
  return freq, psd

mean_path  = os.path.join(Q_directory, "Qmean.vtr")
rhom       = getScalar(mean_path, Nx, Ny, Nz, 'rho')
um, vm, wm = getVector(mean_path, Nx, Ny, Nz, 'velocity')
pm         = getScalar(mean_path, Nx, Ny, Nz, 'p')

Q  = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)

uf = np.zeros((100, num_files), dtype=np.float32)
pf = np.zeros((100, num_files), dtype=np.float32)

yp     = np.zeros(Ny, dtype=np.float32)
rhorms = np.zeros(Ny, dtype=np.float32)
urms   = np.zeros(Ny, dtype=np.float32)
vrms   = np.zeros(Ny, dtype=np.float32)
uvrms  = np.zeros(Ny, dtype=np.float32)
prms   = np.zeros(Ny, dtype=np.float32)
Mrms   = np.zeros(Ny, dtype=np.float32)
tke    = np.zeros(Ny, dtype=np.float32)
ruu    = np.zeros(Ny, dtype=np.float32)
rvv    = np.zeros(Ny, dtype=np.float32)
ruv    = np.zeros(Ny, dtype=np.float32)
mdudy  = np.ones(Ny, dtype=np.float32)
up     = np.zeros(Ny, dtype=np.float32)

itr = 0
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_directory, Q_file)
  Q[0,:,:,:]                         = getScalar(file_path, Nx, Ny, Nz, 'rho')
  Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:] = getVector(file_path, Nx, Ny, Nz, 'velocity')
  Q[4,:,:,:]                         = getScalar(file_path, Nx, Ny, Nz, 'p')
  yplus, _, _, tw, ut, uplus = non_dim_tbl(Q, x, y, z)
  yp += yplus
  up += uplus
  # calc RMS
  Mt = ut / np.sqrt(1.4e0 * np.mean(Q[4,:,0,:] / Q[0,:,0,:]))
  rhorms += np.mean(np.mean(np.sqrt((Q[0,:,:,:] - rhom)**2) / (1.4e0 * np.mean(Q[0,:,0,:]) * Mt**2), axis=0), axis=-1)
  urms   += np.mean(np.mean(np.sqrt((Q[1,:,:,:] - um)**2) / ut, axis=0), axis=-1)
  vrms   += np.mean(np.mean(np.sqrt((Q[2,:,:,:] - vm)**2) / ut, axis=0), axis=-1)
  uvrms  += np.mean(np.mean((Q[1,:,:,:] - um) * (Q[2,:,:,:]- vm) / ut**2, axis=0), axis=-1)
  prms   += np.mean(np.mean(np.sqrt((Q[4,:,:,:] - pm)**2) / tw, axis=0), axis=-1)
  # calc TKE
  tke += np.mean(np.mean(0.5e0 * ((Q[1,:,:,:] - um)**2 + (Q[2,:,:,:] - vm)**2 + (Q[3,:,:,:] - wm)**2) / (ut**2), axis=0), axis=-1)
  # calc Reynolds Stress
  ruu += np.mean(np.mean(np.sqrt(Q[0,:,:,:] * (Q[1,:,:,:] - um)**2 / tw), axis=0), axis=-1)
  rvv += np.mean(np.mean(np.sqrt(Q[0,:,:,:] * (Q[2,:,:,:] - vm)**2 / tw), axis=0), axis=-1)
  ruv += np.mean(np.mean(Q[0,:,:,:] * (Q[1,:,:,:] - um) * (Q[2,:,:,:] - vm) / tw, axis=0), axis=-1)
  # shear stress balance
  mu    = Sutherland(np.mean(np.mean(Q[4,:,:,:] / (287.03e0 * Q[0,:,:,:]), axis=0), axis=-1))
  dudy  = np.mean(np.mean(-Q[1,:,:-1,:] + Q[1,:,1:,:], axis=0), axis=-1)
  mdudy[1:] += 0.5e0 * (mu[:-1] + mu[1:]) * dudy / (-y[:-1] + y[1:]) / tw
  # PSD
  uf[:,itr] = (Q[1,0:20:2,nyp,0:100:10] - um[0:20:2,nyp,0:100:10]).flatten()
  pf[:,itr] = (Q[4,0:20:2,nyp,0:100:10] - pm[0:20:2,nyp,0:100:10]).flatten()
  itr += 1

yp     = yp     / np.float32(num_files)
up     = up     / np.float32(num_files)
rhorms = rhorms / np.float32(num_files)
urms   = urms   / np.float32(num_files)
vrms   = vrms   / np.float32(num_files)
uvrms  = uvrms  / np.float32(num_files)
prms   = prms   / np.float32(num_files)
tke    = tke    / np.float32(num_files)
ruu    = ruu    / np.float32(num_files)
rvv    = rvv    / np.float32(num_files)
ruv    = ruv    / np.float32(num_files)
mdudy  = mdudy  / np.float32(num_files)

save_path = os.path.join(Q_directory, "turb_stat.d")
delta     = 2.e-3
with open(save_path, "w", encoding="UTF-8") as f:
  print("# y      yplus      uplus      urms      vrms       uvrms     tke     ruu      rvv     ruv    mdudy    rhorms     prms", file=f)
  for j in range(Ny):
    print(f'{y[j]/delta:.3e}', f'{yp[j]:.3e}', f'{up[j]:.3e}', f'{urms[j]:.3e}', f'{vrms[j]:.3e}', \
    f'{uvrms[j]:.3e}', f'{tke[j]:.3e}', f'{ruu[j]:.3e}', f'{rvv[j]:.3e}', f'{ruv[j]:.3e}', f'{mdudy[j]:.3e}', \
    f'{rhorms[j]:.3e}', f'{prms[j]:.3e}', file=f)

del Q, rhom, um, vm, wm, pm

t = np.linspace(0.e0, endT, num_files)

freq   = np.zeros(num_files, dtype=np.float32)
u_psd  = np.zeros(num_files, dtype=np.float32)
p_psd  = np.zeros(num_files, dtype=np.float32)

for i in range(len(uf[:,0])):
  freq, u_psd += PSD(uf[i,:], t)
  freq, p_psd += PSD(pf[i,:], t)

freq  /= np.float32(len(uf[:,0]))
u_psd /= np.float32(len(uf[:,0]))
p_psd /= np.float32(len(uf[:,0]))

save_path = os.path.join(Q_directory, "PSD.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# k       u_PSD       p_psd", file=f)
  for i in range(1,int(0.5*num_files)):
    print(f'{freq[i]:.3e}', f'{u_psd[i]:.3e}', f'{p_psd[i]:.3e}', file=f)

