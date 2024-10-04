import numpy as np
import os
import re
from scipy.fftpack import fft, ifft, fftfreq
from scipy import signal
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl, Sutherland

Q_directory = "../3D_solver/TBL/data"
#Q_directory = "../../../../../mnt/data1/TBL_5thHRSLAU2_20240922"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)

# parameters
Rgas  = 287.03e0
gamma = 1.4e0
endT  = 2.087e-3
nyp   = 48

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

Q = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)

rho = np.zeros((Nz,Ny,Nx), dtype=np.float32)
u   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
v   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
w   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
p   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
uf  = np.zeros((Nz,Ny,Nx), dtype=np.float32)
vf  = np.zeros((Nz,Ny,Nx), dtype=np.float32)
wf  = np.zeros((Nz,Ny,Nx), dtype=np.float32)
# Farvre average
rhoF = np.zeros((Nz,Ny,Nx), dtype=np.float32)
uF   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
vF   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
wF   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
pF   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
TF   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
TtF  = np.zeros((Nz,Ny,Nx), dtype=np.float32)

uf_PSD = np.zeros((100, num_files), dtype=np.float32)
pf_PSD = np.zeros((100, num_files), dtype=np.float32)

rhorms = np.zeros(Ny, dtype=np.float32)
urms   = np.zeros(Ny, dtype=np.float32)
vrms   = np.zeros(Ny, dtype=np.float32)
uvrms  = np.zeros(Ny, dtype=np.float32)
prms   = np.zeros(Ny, dtype=np.float32)
Mrms   = np.zeros(Ny, dtype=np.float32)
tke    = np.zeros(Ny, dtype=np.float32)
ruu    = np.zeros(Ny, dtype=np.float32)
rvv    = np.zeros(Ny, dtype=np.float32)
rww    = np.zeros(Ny, dtype=np.float32)
ruv    = np.zeros(Ny, dtype=np.float32)
mdudy  = np.zeros(Ny, dtype=np.float32)

dx = -x[:-2] + x[2:]
dy = -y[:-2] + y[2:]
dz = -z[:-2] + z[2:]

# non dim
rhow = rhom[:,0,:]
Q[0,:,:,:], Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:], Q[4,:,:,:] = rhom, um, vm, wm, pm
yp, _, _, tw, ut, up = non_dim_tbl(Q, x, y, z)
Mt = ut / np.sqrt(1.4e0 * np.mean(pm[:,0,:] / rhow))
print("ut, tw ",ut, tw)
# shear stress balance
mu = Sutherland(np.mean(np.mean(pm / (Rgas * rhom), axis=0), axis=-1))
du = np.mean(np.mean(-um[:,:-1,:] + um[:,1:,:], axis=0), axis=-1)
mdudy[0] = 1.e0
mdudy[1:] = 0.5e0 * (mu[:-1] + mu[1:]) * du / (-y[:-1] + y[1:]) / tw

itr = 0
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_directory, Q_file)
  if file_path == os.path.join(Q_directory, "Qmean.vtr"):
    continue
  rho     = getScalar(file_path, Nx, Ny, Nz, 'rho')
  u, v, w = getVector(file_path, Nx, Ny, Nz, 'velocity')
  p       = getScalar(file_path, Nx, Ny, Nz, 'p')
  uf = u - um
  vf = v - vm
  wf = w - wm
  # calc RMS
  rhorms += np.mean(np.mean(np.sqrt((rho - rhom)**2) / (gamma * np.mean(rhow) * Mt**2), axis=0), axis=-1)
  urms   += np.mean(np.mean(np.sqrt(uf**2) / ut, axis=0), axis=-1)
  vrms   += np.mean(np.mean(np.sqrt(vf**2) / ut, axis=0), axis=-1)
  uvrms  += np.mean(np.mean(uf * vf / ut**2, axis=0), axis=-1)
  prms   += np.mean(np.mean(np.sqrt((p - pm)**2) / tw, axis=0), axis=-1)
  # calc TKE
  tke += np.mean(np.mean(0.5e0 * (uf**2 + vf**2 + wf**2) / (ut**2), axis=0), axis=-1)
  # calc Reynolds Stress
  ruu += np.mean(np.mean(rhom * uf**2 / (rhow * ut**2), axis=0), axis=-1)
  rvv += np.mean(np.mean(rhom * vf**2 / (rhow * ut**2), axis=0), axis=-1)
  rww += np.mean(np.mean(rhom * wf**2 / (rhow * ut**2), axis=0), axis=-1)
  ruv += np.mean(np.mean(rhom * uf * vf / tw, axis=0), axis=-1)
  # Favre average
  rhoF += rho
  uF   += u
  vF   += v
  wF   += w
  pF   += p
  TF   += (p / (Rgas * rho))
  TtF  += (p / (Rgas * rho)) + 0.5e0 * (gamma - 1.e0) * (u**2 + v**2 + w**2) / gamma
  # PSD
  uf_PSD[:,itr] = (u[0:20:2,nyp,0:100:10] - um[0:20:2,nyp,0:100:10]).flatten()
  pf_PSD[:,itr] = (p[0:20:2,nyp,0:100:10] - pm[0:20:2,nyp,0:100:10]).flatten()
  itr += 1

rhorms /= np.float32(itr)
urms   /= np.float32(itr)
vrms   /= np.float32(itr)
uvrms  /= np.float32(itr)
prms   /= np.float32(itr)
tke    /= np.float32(itr)
ruu    /= np.float32(itr)
rvv    /= np.float32(itr)
rww    /= np.float32(itr)
ruv    /= np.float32(itr)
rhoF   /= np.float32(itr)
uF     /= np.float32(itr)
vF     /= np.float32(itr)
wF     /= np.float32(itr)
pF     /= np.float32(itr)
TF     /= np.float32(itr)
TtF    /= np.float32(itr)

# save Favre average
rhoF_path = os.path.join(Q_directory, "rhoF")
uF_path   = os.path.join(Q_directory, "uF")
vF_path   = os.path.join(Q_directory, "vF")
wF_path   = os.path.join(Q_directory, "wF")
pF_path   = os.path.join(Q_directory, "pF")
TF_path   = os.path.join(Q_directory, "TF")
TtF_path  = os.path.join(Q_directory, "TtF")
np.save(rhoF_path, rhoF)
np.save(uF_path,   uF)
np.save(vF_path,   vF)
np.save(wF_path,   wF)
np.save(pF_path,   pF)
np.save(TF_path,   TF)
np.save(TtF_path,  TtF)

save_path = os.path.join(Q_directory, "turb_stat.d")

delta = 2.e-3
for j in range(Ny):
  if np.mean(um[:,j,:]) >= 0.99e0 * np.mean(um[:,-1,:]):
    delta = y[j]
    break

print("delta ", delta)
with open(save_path, "w", encoding="UTF-8") as f:
  print("# y      yplus      uplus      urms      vrms       uvrms     tke      ruu       rvv      rww      ruv     mdudy    rhorms     prms", file=f)
  for j in range(Ny):
    print(f'{y[j]/delta:.3e}', f'{yp[j]:.3e}', f'{up[j]:.3e}', f'{urms[j]:.3e}', f'{vrms[j]:.3e}', \
    f'{uvrms[j]:.3e}', f'{tke[j]:.3e}', f'{ruu[j]:.3e}', f'{rvv[j]:.3e}', f'{rww[j]:.3e}', f'{ruv[j]:.3e}', f'{mdudy[j]:.3e}', \
    f'{rhorms[j]:.3e}', f'{prms[j]:.3e}', file=f)

del Q, rhom, um, vm, wm, pm, uf, vf

t = np.linspace(0.e0, endT, num_files)

freq   = np.zeros((100,num_files), dtype=np.float32)
u_psd  = np.zeros((100,num_files), dtype=np.float32)
p_psd  = np.zeros((100,num_files), dtype=np.float32)

for i in range(len(uf_PSD[:,0])):
  freq[i,:], u_psd[i,:] = PSD(uf_PSD[i,:], t)
  freq[i,:], p_psd[i,:] = PSD(pf_PSD[i,:], t)

freq  = np.mean(freq,  axis=0)
u_psd = np.mean(u_psd, axis=0)
p_psd = np.mean(p_psd, axis=0)

save_path = os.path.join(Q_directory, "PSD.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# k       u_PSD       p_psd", file=f)
  for i in range(1,int(0.5*num_files)):
    print(f'{freq[i]:.3e}', f'{u_psd[i]:.3e}', f'{p_psd[i]:.3e}', file=f)

