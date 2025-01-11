import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mod.mod_Poisson import DCT_x_Dirichlet_y_Neumann


nx = 33
ny = 33
Lx = 2.e0
Ly = 2.e0
dx = Lx / (nx-1)
dy = Ly / (ny-1)
x  = np.linspace(0.5e0 * dx, Lx - 0.5e0 * dx, nx, dtype=np.float32)
y  = np.linspace(0.5e0 * dy, Ly - 0.5e0 * dy, ny, dtype=np.float32)
xc = 0.5e0 * (x[:-1] + x[1:])
yc = 0.5e0 * (y[:-1] + y[1:])

f = np.zeros((ny-1,nx-1), dtype=np.float32)
for j in range(ny-1):
  for i in range(nx-1):
    f[j,i] = -2.e0 * np.pi**2 * np.cos(np.pi * xc[i]) * np.cos(np.pi * yc[j])


# theoretical
# p[:,0]     = 0.e0
# p[:,-1]    = 0.e0
# dpdy[0,:]  = 0.e0
# dpdy[-1,:] = 0.e0
p_th = np.zeros((ny-1,nx-1), dtype=np.float32)
for j in range(ny-1):
  for i in range(nx-1):
    p_th[j,i] = np.sin(np.pi * xc[i]) * np.cos(np.pi * yc[j])


pl = np.zeros(ny-1, dtype=np.float32)
pr = np.zeros(ny-1, dtype=np.float32)
p_dct = DCT_x_Dirichlet_y_Neumann(f, dx, dy, pl, pr)

# plot
contour = [p_th, p_dct]
fig, ax = plt.subplots(1, 2, figsize=(12, 8))
X, Y = np.meshgrid(xc, yc)

for i in range(2):
  ax[i].set_aspect('equal', adjustable='box')
  im = ax[i].contourf(X, Y, contour[i], levels=100, cmap='jet')
  divider = make_axes_locatable(ax[i])
  cax = divider.append_axes('right', '5%', pad='3%')
  fig.colorbar(im, cax=cax, extendrect=True)
  fig.tight_layout()
plt.show()

