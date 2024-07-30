import numpy as np
from numpy.fft import fft
import os
import re
from pyinform.dist import Dist
from pyinform.utils import bin_series
from pyinform.shannon import relative_entropy
from pyinform.shannon import mutual_info
from pyinform import transfer_entropy
from mod.mod_read import gridInfo

data_directory = "./np_data"
Qp_directory   = os.path.join(data_directory, "Qp")
file_name      = "Qm.npy"
file_path      = os.path.join(Qp_directory, file_name)
Qm             = np.load(file_path)

x, y, z, Nx, Ny, Nz = gridInfo(data_directory)

# get y+
file_path = os.path.join(data_directory, "yplus.npy")
yp = np.load(file_path)

yp1 = 40
yp2 = 80
for j in range(len(yp)):
  if yp[j] <= yp1:
    ny1 = j
  elif yp[j] <= yp2:
    ny2 = j
  else:
    break

Qp_files = [f for f in os.listdir(Qp_directory) if f.endswith(".npy")]
num_files = len(Qp_files)

def extract_number(filename):
  match = re.search(r'Q(\d+)\.npy$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Qp_files.sort(key=extract_number)

# fluctuating velocity
u = np.zeros((num_files,3), dtype=np.float32)
# velocity correlation
uu = np.zeros((3,3), dtype=np.float32)

t  = 0
us = 0.e0
vs = 0.e0
ws = 0.e0
for Qp_file in Qp_files:
  file_path = os.path.join(Qp_directory, Qp_file)
  Qp = np.load(file_path)
  for j in range(ny1,ny2):
    dy = -yp[j] + yp[j+1]
    ud = np.mean(Qp[1,:,j,:] - Qm[1,:,j,:])
    vd = np.mean(Qp[2,:,j,:] - Qm[2,:,j,:])
    wd = np.mean(Qp[3,:,j,:] - Qm[3,:,j,:])
    us += ud * dy
    vs += vd * dy
    ws += wd * dy
    # velocity correlation tensor
    uu[0,0] += (ud * ud) * dy
    uu[0,1] += (ud * vd) * dy
    uu[0,2] += (ud * wd) * dy
    uu[1,0] += (vd * ud) * dy
    uu[1,1] += (vd * vd) * dy
    uu[1,2] += (vd * wd) * dy
    uu[2,0] += (wd * ud) * dy
    uu[2,1] += (wd * vd) * dy
    uu[2,2] += (wd * wd) * dy
  u[t,0] = us / (-yp1 + yp2)
  u[t,1] = vs / (-yp1 + yp2)
  u[t,2] = ws / (-yp1 + yp2)
  t += 1

del Qp, Qm

# print fluctuating velocity
info_directory = os.path.join(data_directory, "info")
os.makedirs(info_directory, exist_ok = True)
file_path = os.path.join(info_directory, "uf.npy")
np.save(file_path, u)
file_path = os.path.join(info_directory, "uf.d")
with open(file_path, "w", encoding="UTF-8") as f:
  for t in range(len(u[:,0])):
    print(f'{u[t,0]:.3e}', f'{u[t,1]:.3e}', f'{u[t,2]:.3e}', file=f)

uu = uu / (float(num_files) * (-yp1 + yp2))

print('velocity correlation\n', uu)

file_path = os.path.join(info_directory, "corr.d")
with open(file_path, "w", encoding="UTF-8") as f:
  print(f'{uu[0,0]:.3e}', f'{uu[0,1]:.3e}', f'{uu[0,2]:.3e}', file=f)
  print(f'{uu[1,0]:.3e}', f'{uu[1,1]:.3e}', f'{uu[1,2]:.3e}', file=f)
  print(f'{uu[2,0]:.3e}', f'{uu[2,1]:.3e}', f'{uu[2,2]:.3e}', file=f)
file_path = os.path.join(info_directory, "corr_masked.d")
with open(file_path, "w", encoding="UTF-8") as f:
  print('nan', f'{uu[0,1]:.3e}', f'{uu[0,2]:.3e}', file=f)
  print(f'{uu[1,0]:.3e}', 'nan', f'{uu[1,2]:.3e}', file=f)
  print(f'{uu[2,0]:.3e}', f'{uu[2,1]:.3e}', 'nan', file=f)

KLD = np.zeros((3,3), dtype=np.float32)

for i in range(3):
  for j in range(3):
    obsx, x = np.histogram(u[:,i], bins=10) 
    obsy, y = np.histogram(u[:,j], bins=10)
    p_x  = Dist(obsx)
    p_y  = Dist(obsy)
    KLD[i,j] = relative_entropy(p_x, p_y)

print('KLD\n', KLD)

file_path = os.path.join(info_directory, "KLD.d")
with open(file_path, "w", encoding="UTF-8") as f:
  print(f'{KLD[0,0]:.3e}', f'{KLD[0,1]:.3e}', f'{KLD[0,2]:.3e}', file=f)
  print(f'{KLD[1,0]:.3e}', f'{KLD[1,1]:.3e}', f'{KLD[1,2]:.3e}', file=f)
  print(f'{KLD[2,0]:.3e}', f'{KLD[2,1]:.3e}', f'{KLD[2,2]:.3e}', file=f)

TE = np.zeros((3,3), dtype=np.float32)

for i in range(3):
  for j in range(3):
    xs = bin_series(u[:,i], b=10)
    ys = bin_series(u[:,j], b=10)
    TE[i,j] = transfer_entropy(xs[0], ys[0], k=2)

print('TE\n', TE)

file_path = os.path.join(info_directory, "TE.d")
with open(file_path, "w", encoding="UTF-8") as f:
  print(f'{TE[0,0]:.3e}', f'{TE[0,1]:.3e}', f'{TE[0,2]:.3e}', file=f)
  print(f'{TE[1,0]:.3e}', f'{TE[1,1]:.3e}', f'{TE[1,2]:.3e}', file=f)
  print(f'{TE[2,0]:.3e}', f'{TE[2,1]:.3e}', f'{TE[2,2]:.3e}', file=f)

