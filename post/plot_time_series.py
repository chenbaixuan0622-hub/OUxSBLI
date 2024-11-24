import numpy as np
import os
import re
import matplotlib.pyplot as plt
from tqdm import tqdm
from mod.mod_read import getGrid, getVector


Q_dir = "../../a100"
Q_files    = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
first_path = os.path.join(Q_dir, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)


lx   = 30.e-3
lz   = 1.5e-3
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

Nt  = len(Q_files)
nx1 = int(lx / x[-1] * Nx)
amp = np.zeros(Nt, dtype=np.float32)
itr = 0
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_dir, Q_file)
  if file_path == os.path.join(Q_dir, "TKE.vtr") or \
     file_path == os.path.join(Q_dir, "ReynoldsStress.vtr") or \
     file_path == os.path.join(Q_dir, "Qmean.vtr"):
    continue
  U, _, _  = getVector(file_path, Nx, Ny, Nz, 'velocity')
  u        = U[:,1,nx1] - np.mean(U[:,1,nx1])
  ufft     = np.fft.fft(u)
  amp[itr] = abs(ufft[wn] / (Nz/2))
  itr += 1

name = "amp_" + "x=" + str(lx) + "_k=" + str(lz)
save_path = os.path.join(Q_dir, name + ".d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# t       amp", file=f)
  for i in range(itr):
    print(f'{t[i]:.3e}', f'{amp[i]:.3e}', file=f)

