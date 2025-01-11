import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable



nx = 129
ny = 129
Lx = 2.e0
Ly = 2.e0
dx = Lx / (nx-1)
dy = Ly / (ny-1)
x  = np.linspace(-Lx/2 + 0.5e0 * dx, Lx/2 - 0.5e0 * dx, nx, dtype=np.float32)
y  = np.linspace(-Ly/2 + 0.5e0 * dy, Ly/2 - 0.5e0 * dy, ny, dtype=np.float32)
xc = 0.5e0 * (x[:-1] + x[1:])
yc = 0.5e0 * (y[:-1] + y[1:])

f = np.zeros((ny-1,nx-1), dtype=np.float32)
for j in range(ny-1):
  for i in range(nx-1):
    f[j,i] = np.exp(-10.e0*(xc[i]**2 + yc[j]**2))



# FFT for x and y direction
kx  = np.linspace(0.e0, nx-2, nx-1)
ky  = np.linspace(0.e0, ny-2, ny-1)
mwx = -2.e0 * (1.e0 - np.cos(2.e0 * np.pi * kx / (nx-1))) / dx**2
mwy = -2.e0 * (1.e0 - np.cos(2.e0 * np.pi * ky / (ny-1))) / dy**2
fhat = np.fft.fft2(f)
phat = np.zeros_like(fhat)
for j in range(ny-1):
  for i in range(nx-1):
    if mwx[i] + mwy[j] != 0.e0:
      phat[j,i] = fhat[j,i] / (mwx[i] + mwy[j])
    else:
      phat[j,i] = 0.e0
p_fft = np.real(np.fft.ifft2(phat))



# plot
fig, ax = plt.subplots(figsize=(12, 8))
X, Y = np.meshgrid(xc, yc)
ax.set_aspect('equal', adjustable='box')
im = ax.contourf(X, Y, p_fft, levels=100, cmap='jet')
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', '5%', pad='3%')
fig.colorbar(im, cax=cax, extendrect=True)
fig.tight_layout()

plt.show()

