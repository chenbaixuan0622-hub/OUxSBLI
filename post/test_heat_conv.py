import numpy as np
import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt
import time
from numba import njit
from mod.mod_Poisson_ue import BiCGStab
#from mod.mod_Poisson import BiCGStab


def set_init(nx, ny, nz):
  T = np.zeros((nz,ny,nx), dtype=np.float32)
  T[:,-1,:] = 1.e0
  return T


@njit(cache=True, nogil=True)
def set_bc(T):
  T[:,:,0]  = 0.e0
  T[:,:,-1] = 0.e0
  T[:,0,:]  = 0.e0
  T[0,:,:]  = 0.e0
  T[-1,:,:] = 0.e0
  T[:,-1,:] = 1.e0
  return T


@jax.jit
def set_bc_jax(T):
  T = T.at[:,:,0].set(0.e0)
  T = T.at[:,:,-1].set(0.e0)
  T = T.at[:,0,:].set(0.e0)
  T = T.at[0,:,:].set(T[-2,:,:])
  T = T.at[-1,:,:].set(T[1,:,:])
  T = T.at[:,-1,:].set(1.e0)
  return T

Lx = 1.e0
Ly = 1.e0
Lz = 1.e0

dx = 0.01e0
dy = 0.01e0
dz = 0.01e0

x = np.linspace(0.e0, Lx, int(Lx / dx))
y = np.linspace(0.e0, Ly, int(Ly / dy))
z = np.linspace(0.e0, Lz, int(Lz / dz))

T = set_init(len(x), len(y), len(z))

plt.imshow(T[-1,:,:], cmap='jet', origin='lower')
plt.show()

f = np.zeros_like(T, dtype=np.float32)
dx2dy2 = dx**2 * dy**2
dy2dz2 = dy**2 * dz**2
dz2dx2 = dz**2 * dx**2
err_tol = 1e-16
start = time.time()
#T, err = BiCGStab(T, f, dx2dy2, dy2dz2, dz2dx2, err_tol, set_bc)
#Tj = jax.device_put(T, jax.devices()[0])
#fj = jax.device_put(f, jax.devices()[0])
Tj = jnp.array(T)
fj = jnp.array(f)
dx = jnp.array(-x[:-1] + x[1:])
dy = jnp.array(-y[:-1] + y[1:])
dz = jnp.array(-z[:-1] + z[1:])
Tj, err = BiCGStab(Tj, fj, dx, dy, dz, err_tol, set_bc_jax)
T = np.array(Tj)
jax.block_until_ready(Tj)
end = time.time()

print("err = ", err, " time = ", end - start)

plt.imshow(T[-10,:,:], cmap='jet', origin='lower')
plt.show()


