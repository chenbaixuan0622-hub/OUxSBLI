import jax
import jax.numpy as jnp
from jax import lax
from functools import partial


@jax.jit
def Laplacian(A1, A2, A3, A4, A5, A6, p):
  # p[nz, ny, nx]
  return A1 * p[1:-1,1:-1,:-2] + A2 * p[1:-1,1:-1,2:] + \
         A3 * p[1:-1,:-2,1:-1] + A4 * p[1:-1,2:,1:-1] + \
         A5 * p[:-2,1:-1,1:-1] + A6 * p[2:,1:-1,1:-1] + \
         p[1:-1,1:-1,1:-1]


@partial(jax.jit, static_argnums=(5, 6))
def BiCGStab(p, f, dx, dy, dz, err_tol, set_bc):
  A  = - 1.e0 / (dx[None,None,:-1] * dx[None,None,1:]) \
       - 1.e0 / (dy[None,:-1,None] * dy[None,1:,None]) \
       - 1.e0 / (dz[:-1,None,None] * dz[1:,None,None])
  A1 = 2.e0 / (dx[None,None,:-1] * (dx[None,None,:-1] + dx[None,None,1:])) / A
  A2 = 2.e0 / (dx[None,None,1:]  * (dx[None,None,:-1] + dx[None,None,1:])) / A
  A3 = 2.e0 / (dy[None,:-1,None] * (dy[None,:-1,None] + dy[None,1:,None])) / A
  A4 = 2.e0 / (dy[None,1:,None]  * (dy[None,:-1,None] + dy[None,1:,None])) / A
  A5 = 2.e0 / (dz[:-1,None,None] * (dz[:-1,None,None] + dz[1:,None,None])) / A
  A6 = 2.e0 / (dz[1:,None,None]  * (dz[:-1,None,None] + dz[1:,None,None])) / A
  B  = f[1:-1,1:-1,1:-1] / A

  r0 = jnp.zeros_like(p, dtype=jnp.float32)
  d  = jnp.zeros_like(p, dtype=jnp.float32)
  r  = jnp.zeros_like(p, dtype=jnp.float32)

  r0 = r0.at[1:-1,1:-1,1:-1].set(B - Laplacian(A1, A2, A3, A4, A5, A6, p))

  print(r0)

  r = r0
  d = r0
  err_r = 1.e0

  def condition_fun(x):
    A1, A2, A3, A4, A5, A6, p, r, r0, d, err_r, err_tol = x
    return err_r > err_tol

  def body_fun(x):
    A1, A2, A3, A4, A5, A6, p, r, r0, d, err_r, err_tol = x
    Ad = jnp.zeros_like(p, dtype=jnp.float32)
    Ad = Ad.at[1:-1,1:-1,1:-1].set(Laplacian(A1, A2, A3, A4, A5, A6, d))
    alpha = jnp.sum(r0 * r) / jnp.sum(r0 * Ad)
    s     = r - alpha * Ad
    set_bc(s)
    As = jnp.zeros_like(p, dtype=jnp.float32)
    As = As.at[1:-1,1:-1,1:-1].set(Laplacian(A1, A2, A3, A4, A5, A6, s))
    w     = jnp.sum(As * s) / jnp.sum(As * As)
    dp    = alpha * d + w * s
    p    += dp
    set_bc(p)
    rs    = r
    r     = s - w * As
    beta  = alpha * jnp.sum(r0 * r) / (w * jnp.sum(r0 * rs))
    d     = r + beta * (d - w * Ad)
    set_bc(d)
    err_r = jnp.sqrt(jnp.sum(dp**2) / jnp.sum(p**2))
    return (A1, A2, A3, A4, A5, A6, p, r, r0, d, err_r, err_tol)

  x0 = (A1, A2, A3, A4, A5, A6, p, r, r0, d, err_r, err_tol)
  x = lax.while_loop(condition_fun, body_fun, x0)
  A1, A2, A3, A4, A5, A6, p, r, r0, d, err_r, err_tol = x
  return p, err_r

