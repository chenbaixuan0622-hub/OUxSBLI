import numpy as np
import os
import re
import matplotlib.pyplot as plt
from tqdm import tqdm
from mod.mod_read import getGrid, getVector
from mod.mod_shock import reflected_shock


Q_dir = "../../a100/yp590"
Q_files    = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
first_path = os.path.join(Q_dir, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)


Lx1  = 33.e-3
Lx2  = 35.e-3
lz   = 6.0e-3
endT = 0.1e-3


def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

# calc strongest Fourier component
dz   = -z[0] + z[1]
freq = np.fft.fftfreq(Nz, d=dz)

for k in range(1,int(Nz/2)):
  if 1.e0 / freq[k] < lz:
    wn = k
    break

print("target wavelength is ", 1.e0 / freq[wn-1], 1.e0 / freq[wn])
w1 = (1.e0 / freq[wn] - lz) / (-1.e0 / freq[wn-1] + 1.e0 / freq[wn])
w2 = 1.e0 - w1
print(w1, w2)
print("achieved wavelength is ", 1.e0 / (w1 * freq[wn-1] + w2 * freq[wn]))

Nt  = len(Q_files)
nx1 = int(Lx1 / x[-1] * Nx)
nx2 = int(Lx2 / x[-1] * Nx)
amp = np.zeros(Nt, dtype=np.float32)
itr = 0
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_dir, Q_file)
  if file_path == os.path.join(Q_dir, "TKE.vtr") or \
     file_path == os.path.join(Q_dir, "ReynoldsStress.vtr") or \
     file_path == os.path.join(Q_dir, "Qmean.vtr"):
    continue
  u, v, w  = getVector(file_path, Nx, Ny, Nz, 'velocity')
  ls, _    = reflected_shock(x[nx1:nx2], y, z, 1, u[:,:,nx1:nx2], v[:,:,nx1:nx2], w[:,:,nx1:nx2])
  ls       = (ls - np.mean(ls)) * 1e3
  xfft     = np.fft.fft(ls)
  amp[itr] = w1 * abs(xfft[wn-1] / (Nz/2)) + w2 * abs(xfft[wn] / (Nz/2))
  itr += 1

t = np.linspace(0, endT, Nt)
name = "shock_" + "_k=" + str(lz)
save_path = os.path.join(Q_dir, name + ".d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# t       amp", file=f)
  for i in range(itr):
    print(f'{t[i]:.3e}', f'{amp[i]:.3e}', file=f)

