import numpy as np
import jax
import jax.numpy as jnp
from functools import partial

def set_grid(nx, ny, nz, Lx, Ly, Lz):
  L  = 2.e0 * jnp.pi
  dx = L / jnp.float32(nx-2) 

  x = jnp.zeros(nx, dtype=jnp.float32)
  x[1:-1] = jnp.linspace(0, L, nx-2, dtype=jnp.float32)
  x[0]  = x[1]  - dx
  x[-1] = x[-2] + dx

  y = x
  z = y

  dxdy = jnp.zeros((nz-2,ny-2,nx-2), dtype=jnp.float32)
  dydz = jnp.zeros((nz-2,ny-2,nx-2), dtype=jnp.float32)
  dzdx = jnp.zeros((nz-2,ny-2,nx-2), dtype=jnp.float32)
  dxdy = dx**2
  dydz = dxdy
  dzdx = dydz
  J    = 1.e0 / (dx * dydz)
  return x, y, z, dxdy, dydz, dzdx, J


def set_init(nx, ny, nz, gamma, rho0, M0, x, y, z):
  Q = jnp.zeros((nz, ny, nx, 5), dtype=jnp.float32)
  for k in range(nz):
    for j in range(ny):
      for i in range(nx):
        Q[k,j,i,0] =  rho0
        Q[k,j,i,1] =  rho0 * M0 * jnp.sin(x[i]) * jnp.cos(y[j]) * jnp.cos(z[k])
        Q[k,j,i,2] = -rho0 * M0 * jnp.cos(x[i]) * jnp.sin(y[j]) * jnp.cos(z[k])
        Q[k,j,i,3] = 0.e0
        Q[k,j,i,4] = (1.e0 / gamma + 0.0625e0 * rho0 * (M0**2) * \
        (jnp.cos(2.e0 * x[i]) + jnp.cos(2.e0 * y[j])) * (jnp.cos(2.e0 * z[k]) + 2.e0) \
        ) / (gamma - 1.e0) + 0.5e0 * (Q[k,j,i,1]**2 + Q[k,j,i,2]**2) / Q[k,j,i,0]

  Q[1:-1,1:-1,0,:] = Q[1:-1,1:-1,-2,:]
  Q[1:-1,1:-1,-1,:] = Q[1:-1,1:-1,1,:]
  Q[1:-1,0,1:-1,:] = Q[1:-1,-2,1:-1,:]
  Q[1:-1,-1,1:-1,:] = Q[1:-1,1,1:-1,:]

  Q[1:-1,0,0,:] = Q[1:-1,-2,-2,:]
  Q[1:-1,0,-1,:] = Q[1:-1,-2,1,:]
  Q[1:-1,-1,0,:] = Q[1:-1,1,-2,:]
  Q[1:-1,-1,-1,:] = Q[1:-1,1,1,:]

  Q[0,1:-1,1:-1,:] = Q[-2,1:-1,1:-1,:]
  Q[-1,1:-1,1:-1,:] = Q[1,1:-1,1:-1,:]
  return Q


@jax.jit
def set_bc(Q):
  Q = Q.at[1:-1,1:-1,0,:].set(Q[1:-1,1:-1,-2,:])
  Q = Q.at[1:-1,1:-1,-1,:].set(Q[1:-1,1:-1,1,:])
  Q = Q.at[1:-1,0,1:-1,:].set(Q[1:-1,-2,1:-1,:])
  Q = Q.at[1:-1,-1,1:-1,:].set(Q[1:-1,1,1:-1,:])

  Q = Q.at[1:-1,0,0,:].set(Q[1:-1,-2,-2,:])
  Q = Q.at[1:-1,0,-1,:].set(Q[1:-1,-2,1,:])
  Q = Q.at[1:-1,-1,0,:].set(Q[1:-1,1,-2,:])
  Q = Q.at[1:-1,-1,-1,:].set(Q[1:-1,1,1,:])

  Q = Q.at[0,1:-1,1:-1,:].set(Q[-2,1:-1,1:-1,:])
  Q = Q.at[-1,1:-1,1:-1,:].set(Q[1,1:-1,1:-1,:])
  return Q

