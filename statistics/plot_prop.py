import os
import numpy as np
import mod.mod_plot as myplt
from mod.mod_read import gridInfo

data_directory = "./np_data"
file_name      = "Qrms.npy"

x, y, z, Nx, Ny, Nz = gridInfo(data_directory)

# set grid
X, Y = np.meshgrid(x,y)

data_directory = os.path.join(data_directory, "Qrms")
file_path      = os.path.join(data_directory, file_name)

Q = np.load(file_path)

if len(Q[:,0,0,0]) == 5:
  save_path = os.path.join(data_directory, "rho.png")
  myplt.plot_scalar(X,Y,Q[0,int(0.5*Nz),:,:],save_path)
  save_path = os.path.join(data_directory, "u.png")
  myplt.plot_scalar(X,Y,Q[1,int(0.5*Nz),:,:],save_path)
  save_path = os.path.join(data_directory, "v.png")
  myplt.plot_scalar(X,Y,Q[2,int(0.5*Nz),:,:],save_path)
  save_path = os.path.join(data_directory, "w.png")
  myplt.plot_scalar(X,Y,Q[3,int(0.5*Nz),:,:],save_path)
  save_path = os.path.join(data_directory, "p.png")
  myplt.plot_scalar(X,Y,Q[4,int(0.5*Nz),:,:],save_path)
else:
  save_path = os.path.join(data_directory, "u.png")
  myplt.plot_scalar(X,Y,Q[0,int(0.5*Nz),:,:],save_path)
  save_path = os.path.join(data_directory, "v.png")
  myplt.plot_scalar(X,Y,Q[1,int(0.5*Nz),:,:],save_path)
  save_path = os.path.join(data_directory, "w.png")
  myplt.plot_scalar(X,Y,Q[2,int(0.5*Nz),:,:],save_path)

