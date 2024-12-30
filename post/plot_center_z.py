import numpy as np
import os
import re
import matplotlib.pyplot as plt
from tqdm import tqdm
from mod.mod_read import getGrid, getVector


def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_dir = ["../../a100/yp5", "../../a100/yp30", "../../a100/yp80", "../../a100/yp269", "../../a100/yp590"]
Lx    = 30.e-3
Lz    = 4.e-3
endT  = 0.1e-3
name  = "u_x=" + str(Lx) + "_z=" + str(Lz)


Q_files = [f for f in os.listdir(Q_dir[0]) if f.endswith(".vtr")]
first_path = os.path.join(Q_dir[0], Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)
Nk    = len(Q_dir)
Nt    = len(Q_files)

data  = np.zeros((Nk,Nt))

for k in range(Nk):
  Q_files = [f for f in os.listdir(Q_dir[k]) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)
  nx = int(Lx / x[-1] * Nx)
  nz = int(Lz / z[-1] * Nz)
  dz = -z[0] + z[1]

  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir[k], Q_file)
    if file_path == os.path.join(Q_dir[k], "TKE.vtr") or \
       file_path == os.path.join(Q_dir[k], "ReynoldsStress.vtr") or \
       file_path == os.path.join(Q_dir[k], "Qmean.vtr"):
      continue
    u, _, _ = getVector(file_path, Nx, Ny, Nz, 'velocity')
    data[k,itr] = u[nz,1,nx]
    itr += 1


# store data for animation and calc FFT
Nt   = itr
t    = np.linspace(0.e0, endT, Nt)
dt   = -t[0] + t[1]
freq = np.fft.fftfreq(Nt, d=dt)
amp  = np.zeros((Nk,Nt))

for k in range(Nk):
  amp[k,:] = abs(np.fft.fft(data[k,:Nt]) / (Nt / 2.e0))

# for plot
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['font.size'] = 12

# plot timeseries data
plt.xlabel("t [ms]", fontsize='14', style='italic', color='black')
plt.ylabel("u [m/s]",  fontsize='14', style='italic', color='black')
colors = ['black', 'blue', 'red', 'green', 'purple', 'orange', 'pink']
for k in range(Nk):
  plt.plot(t, data[k,:Nt], color=colors[k])
save_path = os.path.join(dir, name+".png")
plt.savefig(save_path)
plt.show()
plt.close()

# plot FFT
plt.xlabel("Wavelength [ms]", fontsize='14', style='italic', color='black')
plt.ylabel("Amplitude",  fontsize='14', style='italic', color='black')
colors = ['black', 'blue', 'red', 'green', 'purple', 'orange', 'pink']
for k in range(Nk):
  plt.plot(endT / freq[1:Nt//2], amp[k,1:Nt//2], color=colors[k])
#save_path = os.path.join(dir, name+".png")
#plt.savefig(save_path)
plt.show()
plt.close()


