import os
import numpy as np
import mod.mod_plot as myplt
from mod.mod_read import gridInfo

data_directory = "./np_data"
Qp_directory   = os.path.join(data_directory, "Qp")
Qrms_directory = os.path.join(data_directory, "Qrms")
os.makedirs(Qrms_directory, exist_ok = True)

x, y, z, Nx, Ny, Nz = gridInfo(data_directory)

Qp_files  = [f for f in os.listdir(Qp_directory) if f.endswith(".npy")]
num_files = len(Qp_files)

Qm = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
Q2 = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)

for Qp_file in Qp_files:
  file_path = os.path.join(Qp_directory, Qp_file)
  Qm += np.load(file_path)

Qm = Qm / float(num_files)

for Qp_file in Qp_files:
  file_path = os.path.join(Qp_directory, Qp_file)
  Q2 += (np.load(file_path) - Qm)**2

Qrms = np.sqrt(Q2 / float(num_files))

file_path = os.path.join(Qp_directory, "Qm.npy")
np.save(file_path, Qm)
file_path = os.path.join(Qrms_directory, "Qrms.npy")
np.save(file_path, Qrms)

