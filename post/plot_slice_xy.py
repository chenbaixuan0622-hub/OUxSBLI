import numpy as np
import os
import re
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_plot import print_slice


Q_directory = "../../SBLI/SBLI_4delta/stat03ms_04ms_SLAU"
Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

zs = [4.5e-3, 6.e-3]

delta = 2.e-3

def main():
  first_path          = os.path.join(Q_directory, Q_files[0])
  Nx, Ny, Nz, x, y, z = getGrid(first_path)

  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_directory, Q_file)
    if file_path == os.path.join(Q_directory, "TKE.vtr") \
    or file_path == os.path.join(Q_directory, "ReynoldsStress.vtr") \
    or file_path == os.path.join(Q_directory, "Qmean.vtr"):
      continue
    Rho     = getScalar(file_path, Nx, Ny, Nz, 'rho')
    U, V, W = getVector(file_path, Nx, Ny, Nz, 'velocity')
    P       = getScalar(file_path, Nx, Ny, Nz, 'p')

    for i in range(len(zs)):
      for k in range(Nz):
        if z[k] >= zs[i]:
          nz = k
          break

      rho = Rho[nz-1:nz+2,:,:]
      u   =   U[nz-1:nz+2,:,:]
      v   =   V[nz-1:nz+2,:,:]
      w   =   W[nz-1:nz+2,:,:]
      p   =   P[nz-1:nz+2,:,:]

      basename = os.path.splitext(os.path.basename(file_path))[0]
      file_name = 'Q_' + 'z' + f'{z[nz]:.3f}' + '_' + str(itr).zfill(5)
      print_slice(x, y, z[nz-1:nz+2], rho, u, v, w, p, Q_directory, file_name)
    itr += 1

main()

