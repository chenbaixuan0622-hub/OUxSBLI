import jax
import jax.numpy as jnp
from functools import partial

@partial(jax.jit, static_argnums=(0, 1, 2, 3))
def calc_E(nx, ny, nz, gamma, rho, u, v, w, p):
  E = jnp.zeros((nz-2,ny-2,nx-1,5), dtype=jnp.float32)
  for k in range(1,nz-1):
    for j in range(1,ny-1):
      for i in range(nx-1):
        C              = 0.25e0 * (rho[k,j,i] + rho[k,j,i+1]) * (u[k,j,i] + u[k,j,i+1])
        E[k-1,j-1,i,0] = C
        E[k-1,j-1,i,1] = 0.5e0 * C * (u[k,j,i] + u[k,j,i+1]) + 0.5e0 * (p[k,j,i] + p[k,j,i+1])
        E[k-1,j-1,i,2] = 0.5e0 * C * (v[k,j,i] + v[k,j,i+1])
        E[k-1,j-1,i,3] = 0.5e0 * C * (w[k,j,i] + w[k,j,i+1])
        E[k-1,j-1,i,4] = 0.5e0 * C * (p[k,j,i] / rho[k,j,i] + p[k,j,i+1] / rho[k,j,i+1]) / (gamma - 1.e0) \
                   + 0.5e0 * (u[k,j,i] * p[k,j,i+1] + u[k,j,i+1] * p[k,j,i]) \
                   + 0.5e0 * C * (u[k,j,i] * u[k,j,i+1] + v[k,j,i] * v[k,j,i+1] + w[k,j,i] * w[k,j,i+1])
  return E


@partial(jax.jit, static_argnums=(0, 1, 2, 3))
def calc_F(nx, ny, nz, gamma, rho, u, v, w, p):
  F = jnp.zeros((nz-2,ny-1,nx-2,5), dtype=jnp.float32)
  for k in range(1,nz-1):
    for j in range(ny-1):
      for i in range(1,nx-1):
        C              = 0.25e0 * (rho[k,j,i] + rho[k,j+1,i]) * (v[k,j,i] + v[k,j+1,i])
        F[k-1,j,i-1,0] = C
        F[k-1,j,i-1,1] = 0.5e0 * C * (u[k,j,i] + u[k,j+1,i])
        F[k-1,j,i-1,2] = 0.5e0 * C * (v[k,j,i] + v[k,j+1,i]) + 0.5e0 * (p[k,j,i] + p[k,j+1,i])
        F[k-1,j,i-1,3] = 0.5e0 * C * (w[k,j,i] + w[k,j+1,i])
        F[k-1,j,i-1,4] = 0.5e0 * C * (p[k,j,i] / rho[k,j,i] + p[k,j+1,i] / rho[k,j+1,i]) / (gamma - 1.e0) \
                       + 0.5e0 * (v[k,j,i] * p[k,j+1,i] + v[k,j+1,i] * p[k,j,i]) \
                       + 0.5e0 * C * (u[k,j,i] * u[k,j+1,i] + v[k,j,i] * v[k,j+1,i] + w[k,j,i] * w[k,j+1,i])
  return F


@partial(jax.jit, static_argnums=(0, 1, 2, 3))
def calc_G(nx, ny, nz, gamma, rho, u, v, w, p):
  G = jnp.zeros((nz-1,ny-2,nx-2), dtype=jnp.float32)
  for k in range(nz-1):
    for j in range(1,ny-1):
      for i in range(1,nx-1):
        C              = 0.25e0 * (rho[k,j,i] + rho[k+1,j,i]) * (w[k,j,i] + w[k+1,j,i])
        G[k,j-1,i-1,0] = C
        G[k,j-1,i-1,1] = 0.5e0 * C * (u[k,j,i] + u[k+1,j,i])
        G[k,j-1,i-1,2] = 0.5e0 * C * (v[k,j,i] + v[k+1,j,i])
        G[k,j-1,i-1,3] = 0.5e0 * C * (w[k,j,i] + w[k+1,j,i]) + 0.5e0 * (p[k,j,i] + p[k+1,j,i])
        G[k,j-1,i-1,4] = 0.5e0 * C * (p[k,j,i] / rho[k,j,i] + p[k+1,j,i] / rho[k+1,j,i]) / (gamma - 1.e0) \
                       + 0.5e0 * (w[k,j,i] * p[k+1,j,i] + w[k+1,j,i] * p[k,j,i]) \
                       + 0.5e0 * c * (u[k,j,i] * u[k+1,j,i] + v[k,j,i] * v[k+1,j,i] + w[k,j,i] * w[k+1,j,i])
  return G

