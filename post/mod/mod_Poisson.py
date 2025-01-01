import jax
import jax.numpy as jnp
from jax import lax
from functools import partial


@jax.jit
def Laplacian2D(A1, A2, A3, A4, p):
  # p[nz, ny, nx]
  return A1 * p[1:-1,:-2] + A2 * p[1:-1,2:] + \
         A3 * p[:-2,1:-1] + A4 * p[2:,1:-1] + \
         p[1:-1,1:-1]


@jax.jit
def Laplacian3D(A1, A2, A3, A4, A5, A6, p):
  # p[nz, ny, nx]
  return A1 * p[1:-1,1:-1,:-2] + A2 * p[1:-1,1:-1,2:] + \
         A3 * p[1:-1,:-2,1:-1] + A4 * p[1:-1,2:,1:-1] + \
         A5 * (p[:-2,1:-1,1:-1] + p[2:,1:-1,1:-1]) + \
         p[1:-1,1:-1,1:-1]


@jax.jit
def safe_divide(numerator, denominator):
  safe_denom = jnp.where(denominator < 1e-30, 1.e0, denominator)
  result = numerator / safe_denom
  return jnp.where(denominator < 1e-30, 0.e0, result)


@partial(jax.jit, static_argnums=(4, 5))
def BiCGStab2D(p, f, dx, dy, err_tol, set_bc):
  Ax = 0.5e0 * (dx[:-1] + dx[1:]) * dx[:-1] * dx[1:]
  Ay = 0.5e0 * (dy[:-1] + dy[1:]) * dy[:-1] * dy[1:]
  A  = - (dx[None,:-1] + dx[None,1]) / Ax[None,:] \
       - (dy[:-1,None] + dy[1,None]) / Ay[:,None]
  A1 = dx[None,1:]  / (A * Ax[None,:]) 
  A2 = dx[None,:-1] / (A * Ax[None,:])
  A3 = dy[1:,None]  / (A * Ay[:,None])
  A4 = dy[:-1,None] / (A * Ay[:,None])
  B  = f[1:-1,1:-1] / A

  r0 = jnp.zeros_like(p)
  r0 = r0.at[1:-1,1:-1].set(B - Laplacian2D(A1, A2, A3, A4, p))

  r = r0
  d = r0
  err_r = 1.e0

  def condition_fun(x):
    A1, A2, A3, A4, p, r, r0, d, err_r, err_tol = x
    return err_r > err_tol

  def body_fun(x):
    A1, A2, A3, A4, p, r, r0, d, err_r, err_tol = x
    Ad = jnp.zeros_like(p)
    Ad = Ad.at[1:-1,1:-1].set(Laplacian2D(A1, A2, A3, A4, d))
    alpha = safe_divide(jnp.sum(r0 * r), jnp.sum(r0 * Ad))
    s     = r - alpha * Ad
    set_bc(s)
    As = jnp.zeros_like(p)
    As = As.at[1:-1,1:-1].set(Laplacian2D(A1, A2, A3, A4, s))
    w     = safe_divide(jnp.sum(As * s), jnp.sum(As * As))
    dp    = alpha * d + w * s
    p    += dp
    set_bc(p)
    rs    = r
    r     = s - w * As
    beta  = alpha * safe_divide(jnp.sum(r0 * r), w * jnp.sum(r0 * rs))
    d     = r + beta * (d - w * Ad)
    set_bc(d)
    err_r = jnp.sqrt(jnp.sum(dp**2) / jnp.sum(p**2))
    return (A1, A2, A3, A4, p, r, r0, d, err_r, err_tol)

  x0 = (A1, A2, A3, A4, p, r, r0, d, err_r, err_tol)
  x  = lax.while_loop(condition_fun, body_fun, x0)
  A1, A2, A3, A4, p, r, r0, d, err_r, err_tol = x
  return p, err_r


@partial(jax.jit, static_argnums=(5, 6))
def BiCGStab3D(p, f, dx, dy, dz, err_tol, set_bc):
  Ax = 0.5e0 * (dx[:-1] + dx[1:]) * dx[:-1] * dx[1:]
  Ay = 0.5e0 * (dy[:-1] + dy[1:]) * dy[:-1] * dy[1:]
  Az = 0.5e0 * (dz[:-1] + dz[1:]) * dz[:-1] * dz[1:]
  A  = - (dx[None,None,:-1] + dx[None,None,1]) / Ax[None,None,:] \
       - (dy[None,:-1,None] + dy[None,1,None]) / Ay[None,:,None] \
       - (dz[:-1,None,None] + dz[1,None,None]) / Az[:,None,None]
  A1 = dx[None,None,1:]  / (A * Ax[None,None,:]) 
  A2 = dx[None,None,:-1] / (A * Ax[None,None,:])
  A3 = dy[None,1:,None]  / (A * Ay[None,:,None])
  A4 = dy[None,:-1,None] / (A * Ay[None,:,None])
  A5 = dz[1:,None,None]  / (A * Az[:,None,None])
  A6 = dz[:-1,None,None] / (A * Az[:,None,None])
  B  = f[1:-1,1:-1,1:-1] / A

  r0 = jnp.zeros_like(p)
  r0 = r0.at[1:-1,1:-1,1:-1].set(B - Laplacian3D(A1, A2, A3, A4, A5, A6, p))

  r = r0
  d = r0
  err_r = 1.e0

  def condition_fun(x):
    A1, A2, A3, A4, A5, A6, p, r, r0, d, err_r, err_tol = x
    return err_r > err_tol

  def body_fun(x):
    A1, A2, A3, A4, A5, A6, p, r, r0, d, err_r, err_tol = x
    Ad = jnp.zeros_like(p)
    Ad = Ad.at[1:-1,1:-1,1:-1].set(Laplacian3D(A1, A2, A3, A4, A5, A6, d))
    alpha = jnp.sum(r0 * r) / jnp.sum(r0 * Ad)
    s     = r - alpha * Ad
    set_bc(s)
    As = jnp.zeros_like(p)
    As = As.at[1:-1,1:-1,1:-1].set(Laplacian3D(A1, A2, A3, A4, A5, A6, s))
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

