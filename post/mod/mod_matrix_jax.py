import jax
import jax.numpy as jnp
from functools import partial


@partial(jax.jit, static_argnums=(0,1))
def create_matrix(nx, ny, a1, a2, a3, a4):
  # A1, A2, A3, A4 : scalar or 1d array [(nx-2) * (ny-2)]
  n    = nx * ny
  diag = jnp.ones(n, dtype=jnp.float32)

  A1 = jnp.zeros((ny,nx), dtype=jnp.float32)
  A2 = jnp.zeros((ny,nx), dtype=jnp.float32)
  A3 = jnp.zeros((ny,nx), dtype=jnp.float32)
  A4 = jnp.zeros((ny,nx), dtype=jnp.float32)
  
  A1 = A1.at[1:-1,1:-1].set(a1)
  A2 = A2.at[1:-1,1:-1].set(a2)
  A3 = A3.at[1:-1,1:-1].set(a3)
  A4 = A4.at[1:-1,1:-1].set(a4)

  A1 = A1.flatten()
  A2 = A2.flatten()
  A3 = A3.flatten()
  A4 = A4.flatten()

  x1 = A1[1:]
  x2 = A2[:-1]

  y1 = A3[nx:]
  y2 = A4[:-nx]

  A = jnp.diag(diag) + jnp.diag(x1, -1) + jnp.diag(x2, 1) + jnp.diag(y1, -nx) + jnp.diag(y2, nx)
  return A


@partial(jax.jit, static_argnums=(1,2))
def Dirichlet_bc(B, nx, ny, x0, x1, y0, y1):
  # B : 2d array
  B = B.flatten()
  B = B.at[::nx].set(x0)
  B = B.at[nx-1::nx].set(x1)
  B = B.at[:nx].set(y0)
  B = B.at[-nx:].set(y1)
  return B


@partial(jax.jit, static_argnums=(1,2))
def wall_Neumann_bc_A(A, nx, ny):
  for i in range(1,nx-1):
    A = A.at[i,i+nx].set(-1.e0)
  return A


@partial(jax.jit, static_argnums=(1,2))
def wall_Neumann_bc_B(B, nx, ny, x0, x1, y1):
  # B : 2d array
  B = B.flatten()
  B = B.at[::nx].set(x0)
  B = B.at[nx-1::nx].set(x1)
  B = B.at[:nx].set(0.e0)
  B = B.at[-nx:].set(y1)
  return B


@jax.jit
def Poisson_coef(x, y, mz=0.e0):
  dx = -x[:-1] + x[1:]
  dy = -y[:-1] + y[1:]

  Ax = 0.5e0 * (dx[:-1] + dx[1:]) * dx[:-1] * dx[1:]
  Ay = 0.5e0 * (dy[:-1] + dy[1:]) * dy[:-1] * dy[1:]
  A  = - (dx[None,:-1] + dx[None,1]) / Ax[None,:] \
       - (dy[:-1,None] + dy[1,None]) / Ay[:,None] + mz
  A1 = dx[None,1:]  / (A * Ax[None,:]) 
  A2 = dx[None,:-1] / (A * Ax[None,:])
  A3 = dy[1:,None]  / (A * Ay[:,None])
  A4 = dy[:-1,None] / (A * Ay[:,None])
  return A, create_matrix(len(x), len(y), A1, A2, A3, A4)

