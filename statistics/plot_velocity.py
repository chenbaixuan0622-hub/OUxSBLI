import os
import numpy as np
import mod.mod_plot as myplt

data_directory = "./np_data"
file_name      = "Q01000.npy"

# load grid infor
file_path = os.path.join(data_directory, "yplus.npy")
yp        = np.load(file_path)
Ny        = len(yp)

Qvd_directory = os.path.join(data_directory, 'Qvd')

file_path = os.path.join(Qvd_directory, file_name)
Q = np.load(file_path)
file_path = os.path.join(Qvd_directory, "u_y_plus.d")
with open(file_path, "w", encoding="UTF-8") as f:
  print('# yplus   uplus', file=f)
  for j in range(Ny):
    print(f'{yp[j]:.3e}', f'{np.mean(Q[0,:,j,:]):.3e}', file=f)

