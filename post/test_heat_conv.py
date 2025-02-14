import numpy as np
from scipy.linalg import lu_factor, lu_solve
import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt
import time
from mod.mod_Poisson import BiCGStab2D
from mod.mod_matrix import Poisson_coef, Dirichlet_bc, Neumann_bc, check_matrix, Poisson_cyclic


def set_init(nx, ny, nz):
  T = np.zeros((nz,ny,nx), dtype=np.float32)
  T[:,-1,:] = 1.e0
  return T


'''
@jax.jit
def set_bc(T):
  T = T.at[1:-1,1:-1,0].set(0.e0)
  T = T.at[1:-1,1:-1,-1].set(0.e0)
  T = T.at[1:-1,0,:].set(0.e0)
  T = T.at[1:-1,-1,:].set(1.e0)
  T = T.at[0,:,:].set(T[-2,:,:])
  T = T.at[-1,:,:].set(T[1,:,:])
  return T
'''


@jax.jit
def set_bc(T):
  T = T.at[1:-1,0].set(0.e0)
  T = T.at[1:-1,-1].set(0.e0)
  T = T.at[0,:].set(0.e0)
  T = T.at[-1,:].set(1.e0)
  return T


Lx = 1.e0
Ly = 1.e0
Lz = 1.e0

nx = 65
ny = 65
nz = 65

x = np.linspace(0.e0, Lx, nx, dtype=np.float32)
y = np.linspace(0.e0, Ly, ny, dtype=np.float32)
z = np.linspace(0.e0, Lz, nz, dtype=np.float32)

'''
x = np.zeros(nx, dtype=np.float32)
y = np.zeros(ny, dtype=np.float32)

dx = Lx / float(nx-1)
dy = Ly / float(ny-1)

for i in range(nx//2):
  x[i+1] = x[i] + dx
for i in range(nx//2,nx-1):
  x[i+1] = x[i] + 1.01e0 * dx

for j in range(ny//2):
  y[j+1] = y[j] + dy
for j in range(ny//2,ny-1):
  y[j+1] = y[j] + 1.01e0 * dy
'''

'''
T = set_init(nx, ny, nz)

plt.imshow(T[1,:,:], cmap='jet', origin='lower')
plt.show()

f = np.zeros_like(T, dtype=np.float32)
err_tol = 1e-16

# BiCGStab
start = time.time()

Tj  = jnp.array(T, dtype=jnp.float32)
fj  = jnp.array(f, dtype=jnp.float32)
dxj = jnp.array(-x[:-1] + x[1:], dtype=jnp.float32)
dyj = jnp.array(-y[:-1] + y[1:], dtype=jnp.float32)
Tj, err = BiCGStab2D(Tj[0,:,:], fj[0,:,:], dxj, dyj, err_tol, set_bc)
T = np.array(Tj)
jax.block_until_ready(Tj)

end = time.time()

print("err = ", err, " time = ", end - start)

plt.imshow(T[:,:], cmap='jet', origin='lower')
plt.show()
'''

# source term
f = np.zeros((nz,ny,nx), dtype=np.float32)
f[:,-1,:] = 1.e0

'''
# LU decomposition
a, A = Poisson_coef(x, y)
f[1:-1,1:-1] /= a

B = Dirichlet_bc(f, nx, ny, 0.0e0, 0.0e0, 0.0e0, 1.e0)
#A, B = Neumann_bc(A, f, nx, ny)

check_matrix(nx, ny, A, B)

lu_factor, piv = lu_factor(A)
T = lu_solve((lu_factor, piv), B).reshape((ny,nx))
'''

T = Poisson_cyclic(x, y, z, f)

X, Y = np.meshgrid(x, y)
plt.contourf(X, Y, T[2,:,:], np.linspace(0, 1, 100), cmap='jet', extend='both')
plt.show()

