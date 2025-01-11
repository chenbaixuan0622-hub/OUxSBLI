import numpy as np
from scipy.fftpack import dct, idct
from scipy.linalg import solve_banded
import matplotlib.pyplot as plt


# parameters
n   = 30
Lx  = 1.e0
dx  = Lx / n
rhs = np.ones(n)


# boundary conditions
dpl = 1.e0
dpr = 2.e0
rhs[0]  = rhs[0]  + dpl / dx
rhs[-1] = rhs[-1] - dpr / dx


# DCT
rhshat = dct(rhs, type=2, norm='ortho') # DCT-II
kx     = np.linspace(0.e0, n-1, n)
a      = np.pi * kx / n
mw     = 2.e0 * (np.cos(a) - 1.e0) / dx**2
phat   = np.zeros_like(rhshat)
phat[1:] = rhshat[1:] / mw[1:]
phat[0]  = 0.e0
p_dct = idct(phat, type=2, norm='ortho')


# theoretical
x    = np.linspace(0.5e0 * dx, Lx - 0.5e0 * dx, n)
p_th = 0.5e0 * x**2 + x
p_th += p_dct[0] - p_th[0]


# tridiagonal matrix
upper = np.ones(n) / dx**2
diag  = -2.e0 * np.ones(n) / dx**2
lower = np.ones(n) / dx**2

diag[0]   = 1.e0          # Dirichlet boundary condition
diag[-1]  = -1.e0 / dx**2 # Neumann boundary condition
upper[1]  = 0.e0 # Dirichlet boundary condition
# for Banded matrix
upper[0]  = 0.e0
lower[-1] = 0.e0

ab   = np.vstack([upper, diag, lower])
p_td = solve_banded((1,1), ab, rhs)
p_td += p_dct[0] - p_td[0]

# plot
plt.plot(x, p_dct, color='blue')
plt.plot(x, p_td,  color='red', linestyle='dotted')
plt.plot(x, p_th,  'o', color='black', markerfacecolor='white')
plt.show()

