import os
import numpy as np
import mod.mod_plot as myplt
from mod.mod_read import gridInfo

data_directory = "./np_data"
Qrms_directory = os.path.join(data_directory, "Qrms")
file_name      = "Qrms.npy"

x, y, z, Nx, Ny, Nz = gridInfo(data_directory)
file_path = os.path.join(data_directory, "yplus.npy")
yp        = np.load(file_path)

file_path = os.path.join(data_directory, "Q", "Q00100.npy")
Q         = np.load(file_path)
u0        = np.mean(Q[1,:,-1,:])
del Q

file_path = os.path.join(Qrms_directory, file_name)
Qrms      = np.load(file_path)

X, Z = np.meshgrid(x,z)

# pressure rms at wall
save_path = os.path.join(Qrms_directory, "p_at_wall.png")
myplt.plot_scalar(X,Z,Qrms[4,:,0,:],save_path)

X, Y = np.meshgrid(x,y)

save_path = os.path.join(Qrms_directory, "rho.png")
myplt.plot_scalar(X,Y,Qrms[0,int(0.5*Nz),:,:],save_path)
save_path = os.path.join(Qrms_directory, "u.png")
myplt.plot_scalar(X,Y,Qrms[1,int(0.5*Nz),:,:],save_path)
save_path = os.path.join(Qrms_directory, "v.png")
myplt.plot_scalar(X,Y,Qrms[2,int(0.5*Nz),:,:],save_path)
save_path = os.path.join(Qrms_directory, "w.png")
myplt.plot_scalar(X,Y,Qrms[3,int(0.5*Nz),:,:],save_path)
save_path = os.path.join(Qrms_directory, "p.png")
myplt.plot_scalar(X,Y,Qrms[4,int(0.5*Nz),:,:],save_path)

save_path = os.path.join(Qrms_directory, "p_at_wall.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print('# x       prms', file=f)
  for i in range(Nx):
    print(f'{x[i]:.3e}', f'{np.mean(Qrms[4,:,0,i]):.3e}', file=f)

save_path = os.path.join(Qrms_directory, "velocity_rms.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print('# yplus   urms/u0   vrms/u0   wrms/u0', file=f)
  for j in range(Ny):
    print(f'{yp[j]:.3e}', f'{np.mean(Qrms[1,:,j,:])/u0:.3e}', \
    f'{np.mean(Qrms[2,:,j,:])/u0:.3e}', f'{np.mean(Qrms[3,:,j,:])/u0:.3e}', file=f)

