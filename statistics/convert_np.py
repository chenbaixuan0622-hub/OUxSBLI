import os
import numpy as np
from mod.mod_read import getGrid, getScalar, getVector
from mod.mod_turb_stat import non_dim_tbl

# search vtk files
data_directory = "../3D_solver/TBL/KEEP"
save_directory = "./np_data"
os.makedirs(save_directory, exist_ok = True)

vtk_files = [f for f in os.listdir(data_directory) if f.endswith(".vtr")]

# get grid info
first_path          = os.path.join(data_directory, vtk_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)
save_path           = os.path.join(save_directory, "x")
np.save(save_path, x)
save_path           = os.path.join(save_directory, "y")
np.save(save_path, y)
save_path           = os.path.join(save_directory, "z")
np.save(save_path, z)

xp  = np.zeros_like(x, dtype=np.float32)
yp  = np.zeros_like(y, dtype=np.float32)
zp  = np.zeros_like(z, dtype=np.float32)
xps = np.zeros_like(x, dtype=np.float32)
yps = np.zeros_like(y, dtype=np.float32)
zps = np.zeros_like(z, dtype=np.float32)

Q   = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
Qpm = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
Qp2 = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)

Q_directory = os.path.join(save_directory, "Q")
os.makedirs(Q_directory, exist_ok = True)
Qp_directory = os.path.join(save_directory, "Qp")
os.makedirs(Qp_directory, exist_ok = True)
Qvd_directory = os.path.join(save_directory, "Qvd")
os.makedirs(Qvd_directory, exist_ok = True)
Qrms_directory = os.path.join(save_directory, "Qrms")
os.makedirs(Qrms_directory, exist_ok = True)
for vtk_file in vtk_files:
  file_path  = os.path.join(data_directory, vtk_file)
  rho        = getScalar(file_path, Nx, Ny, Nz, 'rho')
  p          = getScalar(file_path, Nx, Ny, Nz, 'p')
  u, v, w    = getVector(file_path, Nx, Ny, Nz, 'velocity')
  Q[0,:,:,:] = rho[:,:,:]
  Q[1,:,:,:] = u[:,:,:]
  Q[2,:,:,:] = v[:,:,:]
  Q[3,:,:,:] = w[:,:,:]
  Q[4,:,:,:] = p[:,:,:]
  xp, yp, zp, Qp, Qvd = non_dim_tbl(Q,x,y,z)
  file_name  = os.path.splitext(file_path.split('/')[-1])[0]
  save_path  = os.path.join(Q_directory, file_name)
  np.save(save_path, Q)
  save_path  = os.path.join(Qp_directory, file_name)
  np.save(save_path, Qp)
  save_path  = os.path.join(Qvd_directory, file_name)
  np.save(save_path, Qvd)
  
  xps += xp
  yps += yp
  zps += zp
  Qpm += Qp
  Qp2 += Qp**2

num_files = len(vtk_files)

xps = xps / float(num_files)
yps = yps / float(num_files)
zps = zps / float(num_files)
Qpm = Qpm / float(num_files)
Qp2 = Qp2 / float(num_files)

del rho
del p
del u
del v
del w
del Q
del Qp
del Qvd

Qrms = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
Qrms = np.sqrt(Qp2 - Qpm**2)

file_path = os.path.join(save_directory, "xplus.npy")
np.save(file_path, xps)
file_path = os.path.join(save_directory, "yplus.npy")
np.save(file_path, yps)
file_path = os.path.join(save_directory, "zplus.npy")
np.save(file_path, zps)
file_path = os.path.join(save_directory, "Qpmean.npy")
np.save(file_path, Qpm)
file_path = os.path.join(Qrms_directory, "Qrms.npy")
np.save(file_path, Qrms)

