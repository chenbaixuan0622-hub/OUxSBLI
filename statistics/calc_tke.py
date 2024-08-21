import numpy as np
import os
import re
from scipy.fftpack import fft, ifft, fftfreq
from scipy import signal
from tqdm import tqdm
from mod.mod_read import getGrid, getVector

Q_directory = "../../../../../mnt/data1/TBL20240819"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)


# set parameters
Nvis = 8
Nlog = 48
dt   = 0.8e-3 / 2000.e0
ut   = 22.e0


def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')


Q_files.sort(key=extract_number)


def PSD(data, dt):
  N    = len(data)
  data = signal.detrend(data)
  domg = 2.e0 * np.pi / float(N * dt)
  yfft = fft(data) / float(N)
  pw   = abs(yfft)**2
  psd  = pw / domg
  # oneside PSD
  psd_oneside = psd[0:int(0.5 * N)+1] * 2.e0
  psd_oneside[int(0.5 * N)] *= 0.5e0
  psd_oneside[0] *= 0.5e0

  xfft   = np.linspace(0.e0, domg * 0.5e0 * float(N), int(0.5 * N)+1)
  xfft_T = 2.e0 * np.pi / xfft
  return xfft_T, psd_oneside


um_path = os.path.join(Q_directory, "um.npy")
vm_path = os.path.join(Q_directory, "vm.npy")
wm_path = os.path.join(Q_directory, "wm.npy")
yp_path = os.path.join(Q_directory, "yp.npy")
yp      = np.load(yp_path)
if os.path.isfile(um_path) & os.path.isfile(vm_path) & os.path.isfile(wm_path):
  um = np.load(um_path)
  vm = np.load(vm_path)
  wm = np.load(wm_path)
else:
  um, vm, wm = getMeanVector(Q_directory, Q_files, Nx, Ny, Nz, "velocity")
  np.save(um_path, um)
  np.save(vm_path, vm)
  np.save(wm_path, wm)

uf  = np.zeros((Nz,Ny,Nx), dtype=np.float32)
vf  = np.zeros((Nz,Ny,Nx), dtype=np.float32)
wf  = np.zeros((Nz,Ny,Nx), dtype=np.float32)
tke = np.zeros((Nz,Ny,Nx), dtype=np.float32)

tke_vis = np.zeros(num_files, dtype=np.float32)
tke_log = np.zeros(num_files, dtype=np.float32)

itr = 0
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_directory, Q_file)
  u, v, w = getVector(file_path, Nx, Ny, Nz, 'velocity')
  uf  = u - um
  vf  = v - vm
  wf  = w - wm
  tke = 0.5e0 * (uf**2 + vf**2 + wf**2) / (ut**2)
  tke_vis[itr] = np.mean(tke[:,Nvis,:])
  tke_log[itr] = np.mean(tke[:,Nlog,:])
  itr += 1


save_path = os.path.join(Q_directory, "tke_yp.d")
delta     = 2.e-3
with open(save_path, "w", encoding="UTF-8") as f:
  print("# y      yplus      tke", file=f)
  for j in range(Ny):
    print(f'{y[j]/delta:.3e}', f'{yp[j]:.3e}', f'{np.mean(tke[:,j,:]):.3e}', file=f)

x_vis, psd_vis = PSD(tke_vis, dt)
x_log, psd_log = PSD(tke_log, dt)

save_path = os.path.join(Q_directory, "tke_psd.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# x_vis  psd_vis x_log  psd_log", file=f)
  for i in range(len(x_vis)):
    print(f'{x_vis[i]:.3e}', f'{psd_vis[i]:.3e}', f'{x_log[i]:.3e}', f'{psd_log[i]:.3e}', file=f)

