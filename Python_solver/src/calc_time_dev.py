import jax
import jax.numpy as jnp
from functools import partial
from calc_conv import calc_E, calc_F, calc_G


@partial(jax.jit, static_argnums=(0, 1, 2, 3, 4))
def calc_R(nx, ny, nz, gamma, dt, dx, dy, dz, J, Q):
  rho = Q[:,:,:,0] * J
  u   = Q[:,:,:,1] / Q[:,:,:,0]
  v   = Q[:,:,:,2] / Q[:,:,:,0]
  w   = Q[:,:,:,3] / Q[:,:,:,0]
  p   = (gamma - 1.e0) * (Q[:,:,:,4] * J - 0.5e0 * rho * (u**2 + v**2 + w**2))
  E   = calc_E(nx, ny, nz, gamma, rho, u, v, w, p)
  F   = calc_F(nx, ny, nz, gamma, rho, u, v, w, p)
  G   = calc_G(nx, ny, nz, gamma, rho, u, v, w, p)
  for k in range(nz-1):
    for j in range(ny-1):
      for i in range(nx-1):
        for l in range(5):
          R[k,j,i,l]  = dt * (dy[j] * dz[k] * (-E[k,j,i,l] + E[k,j,i+1,l]) \
                            + dz[k] * dx[i] * (-F[k,j,i,l] + F[k,j+1,i,l]) \
                            + dx[i] * dy[j] * (-G[k,j,i,l] + G[k+1,j,i,l]))
  return R


@jax.jit
def Runge_Kutta(i, xs):
  nx, ny, nz, dxdy, dydz, dzdx, J, Q = xs
  Qs = jnp.zeros_like(Q, dtype=jnp.float64)
  # 1st step
  R  = calc_R(nx, ny, nz, gamma, dt, dxdy, dydz, dzdx, J, Q)
  Qs = Qs.at[1:-1,1:-1,1:-1,:].set(Q[1:-1,1:-1,1:-1,:] - R)
  Qs = set_bc(Qs)
  # 2nd step
  R  = calc_R(nx, ny, nz, gamma, dt, dxdy, dydz, dzdx, J, Qs)
  Qs = Qs.at[1:-1,1:-1,1:-1,:].set((0.25e0 * (3.e0 * Q[1:-1,1:-1,1:-1,:] + Qs[1:-1,1:-1,1:-1,:] - R)))
  Qs = set_bc(Qs)
  # 3rd step
  R  = calc_R(nx, ny, nz, gamma, dt, dxdy, dydz, dzdx, J, Qs)
  Q  = Q.at[1:-1,1:-1,1:-1,:].set((Q[1:-1,1:-1,1:-1,:] + 2.e0 * Qs[1:-1,1:-1,1:-1,:] - 2.e0 * R) / 3.e0)
  Q  = set_bc(Q)
  return (nx, ny, nz, dxdy, dydz, dzdx, J, Q)

