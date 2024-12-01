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

xs = [34.e-3]

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

    for k in range(len(xs)):
      for i in range(Nx):
        if x[i] >= xs[k]:
          nx = i
          break

      rho = Rho[:,:,nx-1:nx+2]
      u   =   U[:,:,nx-1:nx+2]
      v   =   V[:,:,nx-1:nx+2]
      w   =   W[:,:,nx-1:nx+2]
      p   =   P[:,:,nx-1:nx+2]

      basename = os.path.splitext(os.path.basename(file_path))[0]
      file_name = 'Q_' + 'x' + f'{x[nx]:.4f}' + '_' + str(itr).zfill(5)
      print_slice(x[nx-1:nx+2], y, z, rho, u, v, w, p, Q_directory, file_name)
    itr += 1

main()

