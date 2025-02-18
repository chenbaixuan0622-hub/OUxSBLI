import numpy as np
import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mod.mod_Poisson import Poisson
from mod.mod_BiCGStab import BiCGStab

'''
#1D problems###################################################
# parameters
n   = 30
Lx  = 1.e0
dx  = Lx / n
x   = np.linspace(0.5e0 * dx, Lx - 0.5e0 * dx, n)

f = np.ones(n)
dp1 = 1.e0
dp2 = 2.e0
test1 = Poisson(f, x)
p1 = test1.Neumann(dp1, dp2)

# theoretical
p_th = 0.5e0 * x**2 + x
p_th += p1[0] - p_th[0]

# plot
plt.plot(x, p1, color='blue')
plt.plot(x, p_th,  'o', color='black', markerfacecolor='white')
plt.show()


pl = 1.e0
pr = 2.e0
p2 = test1.Dirichlet(pl, pr)

# theoretical
p_th = 0.5e0 * x**2 + 0.5e0 * x
p_th += p2[0] - p_th[0]

# plot
plt.plot(x, p2, color='blue')
plt.plot(x, p_th,  'o', color='black', markerfacecolor='white')
plt.show()


p3 = test1.Dirichlet_Neumann(dp1, pr)

# theoretical
p_th = 0.5e0 * x**2 + x - 0.5e0

# plot
plt.plot(x, p3, color='blue')
plt.plot(x, p_th,  'o', color='black', markerfacecolor='white')
plt.show()
'''

#2D problems###################################################
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


dpx1 = np.zeros(ny-1)
dpx2 = np.zeros(ny-1)
dpy1 = np.zeros(nx-1)
dpy2 = np.zeros(nx-1)

test2 = Poisson(f, xc, yc)
p4 = test2.x_Neumann_y_Neumann(dpx1, dpx2, dpy1, dpy2)


# plot
contour = [p_th, p4]
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


px1 = np.zeros(ny-1)
px2 = np.zeros(ny-1)
p5 = test2.x_Dirichlet_y_Neumann(px1, px2, dpy1, dpy2)

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal', adjustable='box')
im = ax.contourf(X, Y, p5, levels=100, cmap='jet')
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', '5%', pad='3%')
fig.colorbar(im, cax=cax, extendrect=True)
fig.tight_layout()

plt.show()


nx = 65
ny = 65
Lx = 1.e0
Ly = 1.e0
dx = Lx / (nx-1)
dy = Ly / (ny-1)
x  = np.linspace(0.5e0 * dx, Lx - 0.5e0 * dx, nx, dtype=np.float32)
y  = np.linspace(0.5e0 * dy, Ly - 0.5e0 * dy, ny, dtype=np.float32)
xc = 0.5e0 * (x[:-1] + x[1:])
yc = 0.5e0 * (y[:-1] + y[1:])
f  = np.zeros((ny-1,nx-1), dtype=np.float32)

px1 = np.zeros(ny-1)
px2 = np.zeros(ny-1)
py1 = np.zeros(nx-1)
py2 = np.ones(nx-1)

test3 = Poisson(f, xc, yc)
p6 = test3.x_Dirichlet_y_Dirichlet(px1, px2, py1, py2)

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal', adjustable='box')
X, Y = np.meshgrid(xc, yc)
im = ax.contourf(X, Y, p6, levels=100, cmap='jet')
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', '5%', pad='3%')
fig.colorbar(im, cax=cax, extendrect=True)
fig.tight_layout()

plt.show()


p7 = np.zeros_like(p6, dtype=np.float32)
f  = np.zeros_like(p7, dtype=np.float32)
p7 = jnp.array(p7)
f  = jnp.array(f)
p7 = BiCGStab(p7, f, dx, dy)

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal', adjustable='box')
X, Y = np.meshgrid(xc, yc)
im = ax.contourf(X, Y, p7, levels=100, cmap='jet')
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', '5%', pad='3%')
fig.colorbar(im, cax=cax, extendrect=True)
fig.tight_layout()

plt.show()

