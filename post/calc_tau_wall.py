import numpy as np
import os
import re
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_tau_wall import tau_wall, print_tau


#Q_directory = "../../SBLI/SBLI_4delta/stat03ms_04ms_SLAU"
Q_directory = "../../../../../../media/user/HD-EDS-E/hatayama/TBL/SBLI_1delta/stat03ms_05ms_SLAU"
Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

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
    rho     = getScalar(file_path, Nx, Ny, Nz, 'rho')
    u, _, _ = getVector(file_path, Nx, Ny, Nz, 'velocity')
    p       = getScalar(file_path, Nx, Ny, Nz, 'p')

    tau = tau_wall(Nx, Nz, y, rho, u, p)
    basename = os.path.splitext(os.path.basename(file_path))[0]
    file_name = 'tau' + str(itr).zfill(5)
    print_tau(x, y[0], z, tau, Q_directory, file_name)
    itr += 1

main()

