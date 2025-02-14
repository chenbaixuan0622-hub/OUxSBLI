import numpy as np
from scipy.fftpack import dct, idct
from scipy.linalg import solve_banded, lu_factor, lu_solve
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def dct2(a):
  return dct(dct(a, axis=0, norm='ortho'), axis=1, norm='ortho')


def idct2(a):
  return idct(idct(a, axis=0, norm='ortho'), axis=1, norm='ortho')


nx = 33
ny = 33
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
    f[j,i] = -2.e0 * np.pi**2 * np.cos(np.pi * xc[i]) * np.cos(np.pi * yc[j])


# theoretical
p_th = np.zeros((ny-1,nx-1), dtype=np.float32)
for j in range(ny-1):
  for i in range(nx-1):
    p_th[j,i] = np.cos(np.pi * xc[i]) * np.cos(np.pi * yc[j])



# DCT for x and y direction
kx  = np.linspace(0.e0, nx-2, nx-1)
ky  = np.linspace(0.e0, ny-2, ny-1)
mwx = 2.e0 * (np.cos(np.pi * kx / (nx-1)) - 1.e0) / dx**2
mwy = 2.e0 * (np.cos(np.pi * ky / (ny-1)) - 1.e0) / dy**2
fhat = dct2(f)
phat = fhat / (mwx[None,:] + mwy[:,None])
phat[0,0] = 0.e0
p_dct = idct2(phat)
p_dct += np.mean(p_th) - np.mean(p_dct)


# DCT for y direction
fhat  = dct(f, axis=0, norm='ortho')
# tridiagonal matrix
n     = nx-1
upper = np.ones(n) / dx**2
diag  = -2.e0 * np.ones(n) / dx**2
lower = np.ones(n) / dx**2

diag[0]   = 1.e0          # Dirichlet boundary condition
diag[-1]  = -1.e0 / dx**2 # Neumann boundary condition
upper[1]  = 0.e0          # Dirichlet boundary condition
# for Banded matrix
upper[0]  = 0.e0
lower[-1] = 0.e0

phat = np.zeros_like(fhat)
for k in range(ny-2):
  mw = 2.e0 * (np.cos(np.pi * k / (ny-1)) - 1.e0) / dy**2
  diag[1:-2] -= mw  
  ab = np.vstack([upper, diag, lower])
  phat[k,:] = solve_banded((1,1), ab, fhat[k,:])

p_dct1 = idct(phat, axis=0, norm='ortho')



# plot
contour = [p_th, p_dct, p_dct1]
fig, ax = plt.subplots(1, 3, figsize=(12, 8))
X, Y = np.meshgrid(xc, yc)
for i in range(3):
  ax[i].set_aspect('equal', adjustable='box')
  im = ax[i].contourf(X, Y, contour[i], levels=100, cmap='jet')
  divider = make_axes_locatable(ax[i])
  cax = divider.append_axes('right', '5%', pad='3%')
  fig.colorbar(im, cax=cax, extendrect=True)
  fig.tight_layout()

plt.show()

