import numpy as np
import os
import matplotlib.pyplot as plt
from mod.mod_recurrence import recurrence_plot
from mod.mod_POD import make_data, make_grid


# for plot
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['font.size'] = 12

# parameter
Q_dir   = "../../../../../../media/user/HD-EDS-E/TBL/SBLI_1delta"
Lx1     = 32.e-3
Lx2     = 33.e-3
Ly1     = 1.e-3
Ly2     = 1.1e-3
Lz1     = 0.e-3 
stridex = 4
stridey = 8
endT    = 0.1e-3


# D[space=nx*ny, time]
save_path = os.path.join(Q_dir, "D.npy")
if os.path.isfile(save_path):
  D = np.load(save_path)
  _, _, _, _, _, x, y, t = make_grid(Lx1, Lx2, Ly1, Ly2, stridex, stridey, endT, Q_dir)
else:
  x, y, t, D = make_data(Lx1, Lx2, Ly2, stridex, stridey, endT, Q_dir)
  np.save(save_path, D)


u_series = (D[:100,1] - D[:100,1].mean()) / np.std(D[:100,1])
fig = plt.figure(figsize=(15, 14))
ax = fig.add_subplot(1,2,1)
ax.plot(u_series)

ax = fig.add_subplot(1,2,2)
ax.imshow(recurrence_plot(u_series[:,None]))
plt.show()

