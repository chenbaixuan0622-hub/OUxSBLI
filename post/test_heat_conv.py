import numpy as np
from scipy.linalg import lu_factor, lu_solve
import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt
import time
#from mod.mod_matrix import Poisson_coef, Dirichlet_bc, Neumann_bc, check_matrix
from mod.mod_matrix_jax import Poisson_coef, Dirichlet_bc


Lx = 1.e0
Ly = 1.e0
Lz = 1.e0

nx = 9#65
ny = 9#65

x = np.linspace(0.e0, Lx, nx, dtype=np.float32)
y = np.linspace(0.e0, Ly, ny, dtype=np.float32)

dx = Lx / float(nx-1)
dy = Ly / float(ny-1)

# source term
f = np.zeros((ny,nx), dtype=np.float32)
f[-1,:] = 1.e0

# LU decomposition
a, A = Poisson_coef(x, y)
f[1:-1,1:-1] /= a

B = Dirichlet_bc(f, nx, ny, 0.0e0, 0.0e0, 0.0e0, 1.e0)

lu_factor, piv = lu_factor(A)
T = lu_solve((lu_factor, piv), B).reshape((ny,nx))

X, Y = np.meshgrid(x, y)
plt.contourf(X, Y, T, np.linspace(0, 1, 100), cmap='jet', extend='both')
plt.show()

