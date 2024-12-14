import numpy as np
import os
import re
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_plot import print_slice


#Q_directory = "../../SBLI/SBLI_4delta/stat03ms_04ms_SLAU"
Q_directory = "../../../../../../media/user/HD-EDS-E/hatayama/TBL/SBLI_05delta/stat03ms_09ms_SLAU"
Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]

yps = [5.e0, 30.e0, 80.e0, 269.e0, 590.e0]

delta = 2.e-3

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

def main():
  yp_path = os.path.join(Q_directory, "yp.npy")
  yp = np.load(yp_path)

  itr = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_directory, Q_file)
    if file_path == os.path.join(Q_directory, "TKE.vtr") \
    or file_path == os.path.join(Q_directory, "ReynoldsStress.vtr") \
    or file_path == os.path.join(Q_directory, "Qmean.vtr"):
      continue
    Nx, Ny, Nz, x, y, z = getGrid(file_path)
    Rho     = getScalar(file_path, Nx, Ny, Nz, 'rho')
    U, V, W = getVector(file_path, Nx, Ny, Nz, 'velocity')
    P       = getScalar(file_path, Nx, Ny, Nz, 'p')

    for i in range(len(yps)):
      for j in range(Ny):
        if yp[j] >= yps[i]:
          ny = j
          break

      eta = y[j] / delta

      rho = Rho[:,ny-1:ny+1,:]
      u   =   U[:,ny-1:ny+1,:]
      v   =   V[:,ny-1:ny+1,:]
      w   =   W[:,ny-1:ny+1,:]
      p   =   P[:,ny-1:ny+1,:]

      basename = os.path.splitext(os.path.basename(file_path))[0]
      file_name = 'Q_' + 'yp' + f'{yp[ny]:.4f}' + '_d' + f'{eta:.4f}' + '_' + str(itr).zfill(5)
      print_slice(x, y[ny-1:ny+1], z, rho, u, v, w, p, Q_directory, file_name)
    itr += 1

main()

