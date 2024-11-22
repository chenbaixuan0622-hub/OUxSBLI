import numpy as np
import os
import re
import matplotlib.pyplot as plt
from numba import njit
from tqdm import tqdm
from mod.mod_read import getGrid, getVector
from mod.mod_info import KMeans, causal_map


Q_dir = "../../SBLI/SBLI_4delta/stat03ms_04ms_SLAU"
Q_files    = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
first_path = os.path.join(Q_dir, Q_files[0])
Nx, Ny, Nz, X, Y, Z = getGrid(first_path)

delta = 25.e0
blt   = 2.e-3
endT  = 1.e0

nx1 = int(21.e0 / delta * Nx)
nx2 = int(23.e0 / delta * Nx)

yp1 = 30.e0
history_len = 10
bin = 5

####################################################################

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)
dz = -Z[0] + Z[1]

'''
yp_path = os.path.join(Q_dir, 'yp.npy')
yp = np.load(yp_path)
for j in range(Ny):
  if yp[j] >= yp1:
    ny1 = j
    break
'''
ny1 = 50

@njit(cache=True, fastmath=True, nogil=True)
def corr_z(Nx, Nz, uf, vf, wf, Ruu, Rvv, Rww):
    # z corr
    for k in range(Nz):
      for i in range(Nx):
        Ruu[k,i] += uf[0,i] * uf[k,i]
        Rvv[k,i] += vf[0,i] * vf[k,i]
        Rww[k,i] += wf[0,i] * wf[k,i]

Ruu = np.zeros((Nz,nx2-nx1), dtype=np.float32)
Rvv = np.zeros((Nz,nx2-nx1), dtype=np.float32)
Rww = np.zeros((Nz,nx2-nx1), dtype=np.float32)

uamp = np.zeros((Nz,nx2-nx1), dtype=np.float32)
vamp = np.zeros((Nz,nx2-nx1), dtype=np.float32)
wamp = np.zeros((Nz,nx2-nx1), dtype=np.float32)

itr = 0.e0
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_dir, Q_file)
  if file_path == os.path.join(Q_dir, "TKE.vtr") or \
     file_path == os.path.join(Q_dir, "ReynoldsStress.vtr") or \
     file_path == os.path.join(Q_dir, "Qmean.vtr"):
    continue
  U, V, W = getVector(file_path, Nx, Ny, Nz, 'velocity')

  u = U[:,ny1,nx1:nx2] - np.mean(U[:,ny1,nx1:nx2])
  v = V[:,ny1,nx1:nx2] - np.mean(V[:,ny1,nx1:nx2])
  w = W[:,ny1,nx1:nx2] - np.mean(W[:,ny1,nx1:nx2])

  corr_z(nx2-nx1, Nz, u, v, w, Ruu, Rvv, Rww)

  uamp += abs(np.fft.fft2(u) / (float(Nz) / 2.e0))
  vamp += abs(np.fft.fft2(v) / (float(Nz) / 2.e0))
  wamp += abs(np.fft.fft2(w) / (float(Nz) / 2.e0))
  itr  += 1.e0

Ruu = np.mean(Ruu, axis=-1) / itr
Rvv = np.mean(Rvv, axis=-1) / itr
Rww = np.mean(Rww, axis=-1) / itr
Ruu /= Ruu[0]
Rvv /= Rvv[0]
Rww /= Rww[0]

uamp /= itr
vamp /= itr
wamp /= itr

# calc integral length scale
uils = np.trapz(Ruu[:int(Nz/2)], Z[:int(Nz/2)])
vils = np.trapz(Rvv[:int(Nz/2)], Z[:int(Nz/2)])
wils = np.trapz(Rww[:int(Nz/2)], Z[:int(Nz/2)])

# calc strongest Fourier component
freq = np.fft.fftfreq(Nz, d=dz)
Lz   = 1.e0 / freq[1:int(Nz/2)]
ku = np.argmax(uamp[1:int(Nz/2),0])
kv = np.argmax(vamp[1:int(Nz/2),0])
kw = np.argmax(wamp[1:int(Nz/2),0])

save_path = os.path.join(Q_dir, "corr_spanwise.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# z       Ruu       Rvv        Rww", file=f)
  for k in range(int(Nz/2)):
    print(f'{Z[k]/blt:.3e}', f'{Ruu[k]:.3e}', f'{Rvv[k]:.3e}', f'{Rww[k]:.3e}', file=f)

save_path = os.path.join(Q_dir, "FFT.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# z       uamp      vamp       wamp", file=f)
  for k in range(int(Nz/2)-1):
    print(f'{Lz[k]/blt:.3e}', f'{uamp[k+1,0]:.3e}', f'{vamp[k+1,0]:.3e}', f'{wamp[k+1,0]:.3e}', file=f)

save_path = os.path.join(Q_dir, "spanwise_scale.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("u: integral length scale = ", f'{uils / blt}', file=f)
  print("v: integral length scale = ", f'{vils / blt}', file=f)
  print("w: integral length scale = ", f'{wils / blt}', file=f)
  print("u: the strongest scale = ",   f'{Lz[ku] / blt}', file=f)
  print("v: the strongest scale = ",   f'{Lz[kv] / blt}', file=f)
  print("w: the strongest scale = ",   f'{Lz[kw] / blt}', file=f)

print("u: integral length scale = ", uils / blt)
print("v: integral length scale = ", vils / blt)
print("w: integral length scale = ", wils / blt)
print("u: the strongest scale = ",   Lz[ku] / blt)
print("v: the strongest scale = ",   Lz[kv] / blt)
print("w: the strongest scale = ",   Lz[kw] / blt)

for k in range(1,int(Nz/2)):
  if 1.e0 / freq[k] < uils:
    iu = k
    break

for k in range(1,int(Nz/2)):
  if 1.e0 / freq[k] < vils:
    iv = k
    break

for k in range(1,int(Nz/2)):
  if 1.e0 / freq[k] < wils:
    iw = k
    break

itr = 0
Nt = len(Q_files)
ufreq1 = np.zeros(Nt, dtype=np.float32)
ufreq2 = np.zeros(Nt, dtype=np.float32)
vfreq1 = np.zeros(Nt, dtype=np.float32)
vfreq2 = np.zeros(Nt, dtype=np.float32)
wfreq1 = np.zeros(Nt, dtype=np.float32)
wfreq2 = np.zeros(Nt, dtype=np.float32)
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_dir, Q_file)
  if file_path == os.path.join(Q_dir, "TKE.vtr") or \
     file_path == os.path.join(Q_dir, "ReynoldsStress.vtr") or \
     file_path == os.path.join(Q_dir, "Qmean.vtr"):
    continue
  U, V, W = getVector(file_path, Nx, Ny, Nz, 'velocity')

  u = U[:,ny1,nx1:nx2] - np.mean(U[:,ny1,nx1:nx2])
  v = V[:,ny1,nx1:nx2] - np.mean(V[:,ny1,nx1:nx2])
  w = W[:,ny1,nx1:nx2] - np.mean(W[:,ny1,nx1:nx2])

  ufft = np.fft.fft2(u)
  vfft = np.fft.fft2(v)
  wfft = np.fft.fft2(w)
  ufreq1[itr] = abs(ufft[iu,0])
  ufreq2[itr] = abs(ufft[ku+1,0])
  vfreq1[itr] = abs(vfft[iv,0])
  vfreq2[itr] = abs(vfft[kv+1,0])
  wfreq1[itr] = abs(wfft[iw,0])
  wfreq2[itr] = abs(wfft[kw+1,0])
  itr += 1

del U, V, W, u, v, w, ufft, vfft, wfft, uamp, vamp, wamp

t = np.linspace(0.e0, endT, itr)

save_path = os.path.join(Q_dir, "FFT_signal.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# t       ufreq1    ufreq2    vfreq1    vfreq2    wfreq1    wfreq2", file=f)
  for i in range(itr):
    print(f'{t[i]:.3e}', f'{ufreq1[i]:.3e}', f'{ufreq2[i]:.3e}', f'{vfreq1[i]:.3e}', \
          f'{vfreq2[i]:.3e}', f'{wfreq1[i]:.3e}', f'{wfreq2[i]:.3e}', file=f)

# binned
ufreq1 = KMeans(ufreq1[:itr], bin)
ufreq2 = KMeans(ufreq2[:itr], bin)
vfreq1 = KMeans(vfreq1[:itr], bin)
vfreq2 = KMeans(vfreq2[:itr], bin)
wfreq1 = KMeans(wfreq1[:itr], bin)
wfreq2 = KMeans(wfreq2[:itr], bin)

save_path = os.path.join(Q_dir, "FFT_signal_binned.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# t       ufreq1    ufreq2    vfreq1    vfreq2    wfreq1    wfreq2", file=f)
  for i in range(itr):
    print(f'{t[i]:.3e}', f'{ufreq1[i]:.3e}', f'{ufreq2[i]:.3e}', f'{vfreq1[i]:.3e}', \
          f'{vfreq2[i]:.3e}', f'{wfreq1[i]:.3e}', f'{wfreq2[i]:.3e}', file=f)

V = np.stack([ufreq1, ufreq2, vfreq1, vfreq2, wfreq1, wfreq2], axis=0)

map = causal_map(V, history_len, bin)
save_path = os.path.join(Q_dir, "map.npy")
np.save(save_path, map)
print("causal map")
print(map)

