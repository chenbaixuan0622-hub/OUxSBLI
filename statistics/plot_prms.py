import os
import numpy as np
import mod.mod_plot as myplt
from mod.mod_read import getGrid

Q_directory = "../../../../../media/user/HD-EDS-E/hatayama/TOS/TOS20240830"
prms_path = os.path.join(Q_directory, "prms.npy")
p_path    = os.path.join(Q_directory, "pm.npy")
prms      = np.load(prms_path)
p         = np.load(p_path)

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]

first_path       = os.path.join(Q_directory, Q_files[0])
Nx, _, _, x, _, _ = getGrid(first_path)

N1 = 400
N2 = Nx-1

delta = 2.e-3

x0 = 0.e0

p1 = np.mean(p[:,:,0])

save_path = os.path.join(Q_directory, "prms_wall.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print('# x       prms', file=f)
  for i in range(N1, N2):
    print(f'{(x[i] - x0) / delta:.3e}', f'{np.mean(prms[:,0,i]) / p1:.3e}', file=f)

