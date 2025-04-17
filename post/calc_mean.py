import numpy as np
import os
from tqdm import tqdm
from mod.mod_read import getGrid, getVector, getScalar, extract_number

Q_directory = "../../SBLI/SBLI_16delta/2.2ms/5"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
Q_files.sort(key=extract_number)

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)

Rgas  = 287.03e0
gamma = 1.4e0

rhom = np.zeros((Nz,Ny,Nx), dtype=np.float32)
um   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
vm   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
wm   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
pm   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
# Farvre average
uF   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
vF   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
wF   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
pF   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
TF   = np.zeros((Nz,Ny,Nx), dtype=np.float32)
TtF  = np.zeros((Nz,Ny,Nx), dtype=np.float32)

itr = 0
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_directory, Q_file)
  if file_path == os.path.join(Q_directory, "TKE.vtr") \
  or file_path == os.path.join(Q_directory, "ReynoldsStress.vtr") \
  or file_path == os.path.join(Q_directory, "Qmean.vtr"):
    continue
  rho     = getScalar(file_path, Nx, Ny, Nz, 'rho')
  u, v, w = getVector(file_path, Nx, Ny, Nz, 'velocity')
  p       = getScalar(file_path, Nx, Ny, Nz, 'p')
  # average
  rhom += rho
  um   += u
  vm   += v
  wm   += w
  pm   += p
  # Favre average
  uF   += rho * u
  vF   += rho * v
  wF   += rho * w
  pF   += rho * p
  TF   += rho * (p / (Rgas * rho))
  TtF  += rho * ((p / (Rgas * rho)) + 0.5e0 * (gamma - 1.e0) * (u**2 + v**2 + w**2) / gamma)
  itr  += 1

rhom /= np.float32(itr)
um   /= np.float32(itr)
vm   /= np.float32(itr)
wm   /= np.float32(itr)
pm   /= np.float32(itr)
uF   /= (np.float32(itr) * rhom)
vF   /= (np.float32(itr) * rhom)
wF   /= (np.float32(itr) * rhom)
pF   /= (np.float32(itr) * rhom)
TF   /= (np.float32(itr) * rhom)
TtF  /= (np.float32(itr) * rhom)

# save Favre average
rho_path = os.path.join(Q_directory, "rho")
u_path   = os.path.join(Q_directory, "u")
v_path   = os.path.join(Q_directory, "v")
w_path   = os.path.join(Q_directory, "w")
p_path   = os.path.join(Q_directory, "p")
uF_path  = os.path.join(Q_directory, "uF")
vF_path  = os.path.join(Q_directory, "vF")
wF_path  = os.path.join(Q_directory, "wF")
pF_path  = os.path.join(Q_directory, "pF")
TF_path  = os.path.join(Q_directory, "TF")
TtF_path = os.path.join(Q_directory, "TtF")
np.save(rho_path, rhom)
np.save(u_path,   um)
np.save(v_path,   vm)
np.save(w_path,   wm)
np.save(p_path,   pm)
np.save(uF_path,  uF)
np.save(vF_path,  vF)
np.save(wF_path,  wF)
np.save(pF_path,  pF)
np.save(TF_path,  TF)
np.save(TtF_path, TtF)

