import numpy as np
import os
import re
from scipy.fftpack import fft, ifft, fftfreq
from scipy import signal
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl, Sutherland

Q_directory = "../3D_solver/TBL/data"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)

# parameters
Rgas  = 287.03e0
gamma = 1.4e0
endT  = 2.087e-3

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

rho_path = os.path.join(Q_directory, "rho.npy")
u_path   = os.path.join(Q_directory, "u.npy")
v_path   = os.path.join(Q_directory, "v.npy")
w_path   = os.path.join(Q_directory, "w.npy")
p_path   = os.path.join(Q_directory, "p.npy")

rhom = np.load(rho_path)
um   = np.load(u_path)
vm   = np.load(v_path)
wm   = np.load(w_path)
pm   = np.load(p_path)

Q = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)

rho = np.zeros((Nz,Ny,Nx), dtype=np.float32)
u   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
v   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
w   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
p   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
uf  = np.zeros((Nz,Ny,Nx), dtype=np.float32)
vf  = np.zeros((Nz,Ny,Nx), dtype=np.float32)
wf  = np.zeros((Nz,Ny,Nx), dtype=np.float32)

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

# spanwise correlaion
Rrr = np.zeros(Nz, dtype=np.float32)
Ruu = np.zeros(Nz, dtype=np.float32)
Rvv = np.zeros(Nz, dtype=np.float32)
Rww = np.zeros(Nz, dtype=np.float32)
Rpp = np.zeros(Nz, dtype=np.float32)

dx = -x[:-2] + x[2:]
dy = -y[:-2] + y[2:]
dz = -z[:-2] + z[2:]

nx1 = int(0.5*Nx)
nx2 = int(0.9*Nx)

# non dim
rhow = rhom[:,0,nx1:nx2]
Q[0,:,:,:], Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:], Q[4,:,:,:] = rhom, um, vm, wm, pm
yp, _, _, tw, ut, up = non_dim_tbl(Q, x, y, z)
Mt = ut / np.sqrt(1.4e0 * np.mean(pm[:,0,nx1:nx2] / rhow))
print("ut, tw ",ut, tw)
# shear stress balance
mu = Sutherland(np.mean(pm[:,:,nx1:nx2] / (Rgas * rhom[:,:,nx1:nx2]), axis=(0,2)))
du = np.mean(np.mean(-um[:,:-1,nx1:nx2] + um[:,1:,nx1:nx2], axis=0), axis=-1)
mdudy[0] = 1.e0
mdudy[1:] = 0.5e0 * (mu[:-1] + mu[1:]) * du / (-y[:-1] + y[1:]) / tw

for j in range(Ny):
  if yp[j] > 10.e0:
    nyp = j
    break

itr = 0
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_directory, Q_file)
  if file_path == os.path.join(Q_directory, "Qmean.vtr"):
    continue
  rho     = getScalar(file_path, Nx, Ny, Nz, 'rho')
  u, v, w = getVector(file_path, Nx, Ny, Nz, 'velocity')
  p       = getScalar(file_path, Nx, Ny, Nz, 'p')
  # fluctuation
  rhof = rho - rhom
  uf   = u - um
  vf   = v - vm
  wf   = w - wm
  pf   = p - pm
  # spanwise correlation
  for k in range(Nz):
    Rrr[k] += np.mean(rhof[0,nyp,nx1:nx2] * rhof[k,nyp,nx1:nx2])
    Ruu[k] += np.mean(  uf[0,nyp,nx1:nx2] *   uf[k,nyp,nx1:nx2])
    Rvv[k] += np.mean(  vf[0,nyp,nx1:nx2] *   vf[k,nyp,nx1:nx2])
    Rww[k] += np.mean(  wf[0,nyp,nx1:nx2] *   wf[k,nyp,nx1:nx2])
    Rpp[k] += np.mean(  pf[0,nyp,nx1:nx2] *   pf[k,nyp,nx1:nx2])
  # calc RMS
  rhorms += np.mean(np.sqrt(rhof[:,:,nx1:nx2]**2) / (gamma * np.mean(rhow) * Mt**2), axis=(0,2))
  urms   += np.mean(np.sqrt(uf[:,:,nx1:nx2]**2) / ut, axis=(0,2))
  vrms   += np.mean(np.sqrt(vf[:,:,nx1:nx2]**2) / ut, axis=(0,2))
  uvrms  += np.mean(uf[:,:,nx1:nx2] * vf[:,:,nx1:nx2] / ut**2, axis=(0,2))
  prms   += np.mean(np.sqrt(pf[:,:,nx1:nx2]**2) / tw, axis=(0,2))
  # calc TKE
  tke += np.mean(0.5e0 * (uf[:,:,nx1:nx2]**2 + vf[:,:,nx1:nx2]**2 + wf[:,:,nx1:nx2]**2) / (ut**2), axis=(0,2))
  # calc Reynolds Stress
  ruu += np.mean(rhom[:,:,nx1:nx2] * uf[:,:,nx1:nx2]**2 / (rhow * ut**2), axis=(0,2))
  rvv += np.mean(rhom[:,:,nx1:nx2] * vf[:,:,nx1:nx2]**2 / (rhow * ut**2), axis=(0,2))
  rww += np.mean(rhom[:,:,nx1:nx2] * wf[:,:,nx1:nx2]**2 / (rhow * ut**2), axis=(0,2))
  ruv += np.mean(rhom[:,:,nx1:nx2] * uf[:,:,nx1:nx2] * vf[:,:,nx1:nx2] / tw, axis=(0,2))
  # PSD
  #uf_PSD[:,itr] = (u[0:20:2,nyp,0:100:10] - um[0:20:2,nyp,0:100:10]).flatten()
  #pf_PSD[:,itr] = (p[0:20:2,nyp,0:100:10] - pm[0:20:2,nyp,0:100:10]).flatten()
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
Rrr    /= np.float32(itr)
Ruu    /= np.float32(itr)
Rvv    /= np.float32(itr)
Rww    /= np.float32(itr)
Rpp    /= np.float32(itr)

save_path = os.path.join(Q_directory, "turb_stat.d")

delta = 2.e-3
for j in range(Ny):
  if np.mean(um[:,j,nx1:nx2]) >= 0.99e0 * np.mean(um[:,-1,nx1:nx2]):
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

save_path = os.path.join(Q_directory, "span_corr.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# z       Rrr       Ruu       Rvv       Rww       Rpp", file=f)
  for k in range(Nz):
    print(f'{z[k]/delta:.3e}', f'{Rrr[k]/Rrr[0]:.3e}', f'{Ruu[k]/Ruu[0]:.3e}', f'{Rvv[k]/Rvv[0]:.3e}', f'{Rww[k]/Rww[0]:.3e}', f'{Rpp[k]/Rpp[0]:.3e}', file=f)

'''
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
'''

