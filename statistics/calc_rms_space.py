import numpy as np
import scipy as sp
import os
import re
import matplotlib.pyplot as plt
# read module
from mod.mod_read import getGrid, getVector

Q_directory = "../3D_solver/KHI/HRSLAU2"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)

u    = np.zeros((Nz,Ny,Nx), dtype=np.float32)
urms = np.zeros(num_files, dtype=np.float32)
t    = np.linspace(0.e0, 5.e0, num_files)

itr = 0
for Q_file in Q_files:
  file_path = os.path.join(Q_directory, Q_file)
  u, _, _  = getVector(file_path, Nx, Ny, Nz, 'velocity')
  urms[itr] = np.sqrt(np.mean(u**2) - np.mean(u)**2)
  itr += 1

save_path = os.path.join(Q_directory, "urms_space.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# t       urms", file=f)
  for i in range(num_files):
    print(f'{t[i]:.3e}', f'{urms[i]/urms[0]:.3e}', file=f)

plt.plot(t,urms/urms[0])
plt.show()

