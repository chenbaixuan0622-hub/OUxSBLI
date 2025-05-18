import numpy as np
import os
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar
from mod.SRA import SRA, print_SRA


Q_directory = "../../../../../mnt/data1/TBL_Re1400"
Q_files = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
Rgas  = 287.03e0
gamma = 1.4e0


def StrongReynoldsAnalogy(Nx, Ny, Nz, x, y, z, U, V, T, Tt):
  uT  = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  vT  = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  uv  = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  vT1 = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  vTt = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  T2  = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  Tt2 = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  itr = 0.e0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_directory, Q_file)
    if file_path == os.path.join(Q_directory, "TKE.vtr") \
    or file_path == os.path.join(Q_directory, "ReynoldsStress.vtr") \
    or file_path == os.path.join(Q_directory, "SRA.vtr") \
    or file_path == os.path.join(Q_directory, "Qrms.vtr") \
    or file_path == os.path.join(Q_directory, "Qmean.vtr"):
      continue
    rho     = getScalar(file_path, Nx, Ny, Nz, 'rho')
    u, v, w = getVector(file_path, Nx, Ny, Nz, 'velocity')
    p       = getScalar(file_path, Nx, Ny, Nz, 'p')
    t       = p / (rho * Rgas)
    tt      = t + (gamma - 1.e0) / gamma * 0.5e0 * (u**2 + v**2 + w**2)
    uf  = u - U
    vf  = v - V
    Tf  = t - T
    Ttf = tt - Tt
    # calc Strong Reynolds Analogy
    SRA(uf, vf, Tf, Ttf, uT, vT, uv, vT1, vTt, T2, Tt2)
    itr += 1.e0
  uT  /= itr
  vT  /= itr
  uv  /= itr
  vT1 /= itr
  vTt /= itr
  T2  /= itr
  Tt2 /= itr
  #factor1 = 1.e0 - vTt / vT1
  #factor2 = Tt2 / (2.e0 * T2)
  print_SRA(x, y, z, uT, vT, uv, Q_directory, "SRA")


def main():
  first_path          = os.path.join(Q_directory, Q_files[0])
  Nx, Ny, Nz, x, y, z = getGrid(first_path)
  u_path  = os.path.join(Q_directory, "uF.npy")
  v_path  = os.path.join(Q_directory, "vF.npy")
  T_path  = os.path.join(Q_directory, "TF.npy")
  Tt_path = os.path.join(Q_directory, "TtF.npy")
  U  = np.load(u_path)
  V  = np.load(v_path)
  T  = np.load(T_path)
  Tt = np.load(Tt_path)
  StrongReynoldsAnalogy(Nx, Ny, Nz, x, y, z, U, V, T, Tt)


main()

