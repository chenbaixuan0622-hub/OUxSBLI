import jax
import jax.numpy as jnp


@jax.jit
def set_bc(p):
  p = p.at[1:-1,0].set(0.e0)
  p = p.at[1:-1,-1].set(0.e0)
  p = p.at[0,:].set(0.e0)
  p = p.at[-1,:].set(1.e0)
  return p


@jax.jit
def Poisson(A1, A2, p):
  return A1 * (p[1:-1,:-2] + p[1:-1,2:]) + \
         A2 * (p[:-2,1:-1] + p[2:,1:-1]) + p[1:-1,1:-1]


@jax.jit
def BiCGStab(p, f, dx, dy, mu=0.e0):
  A  = jnp.float32(-2.e0 * (dx**2 + dy**2) / (dx**2 * dy**2))
  A1 = jnp.float32(1.e0 / (dx**2 * A))
  A2 = jnp.float32(1.e0 / (dy**2 * A))
  B  = f / A
  p  = set_bc(p)
  r0 = jnp.zeros_like(p, dtype=jnp.float32)
  r0 = r0.at[1:-1,1:-1].set(B[1:-1,1:-1] - Poisson(A1, A2, p))
  r  = r0
  d  = r0
  err = 1.e0
  err_tol = 1e-6

  @jax.jit
  def safe_divide(a, b):
    return a / b
    '''
    if jnp.abs(b) < 1.e-30:
      return 0.e0
    else:
      return a / b
    '''
  
  @jax.jit
  def condition(x):
    err = x[-1]
    return err > err_tol

  @jax.jit
  def body(x):
    p, d, r0, r, A1, A2, err = x
    Ad = jnp.zeros_like(p, dtype=jnp.float32)
    Ad = Ad.at[1:-1,1:-1].set(Poisson(A1, A2, d))
    alpha = safe_divide(jnp.sum(r0 * r), jnp.sum(r0 * Ad))
    s  = r - alpha * Ad
    s  = set_bc(s)
    As = jnp.zeros_like(p, dtype=jnp.float32)
    As = As.at[1:-1,1:-1].set(Poisson(A1, A2, s))
    w  = safe_divide(jnp.sum(As * s), jnp.sum(As * As))
    dp = alpha * d + w * s
    p += dp
    p  = set_bc(p)
    rs = r
    r  = s - w * As
    beta = alpha * safe_divide(jnp.sum(r0 * r), (w * jnp.sum(r0 * rs)))
    d  = r + beta * (d - w * Ad)
    d  = set_bc(d)
    err = jnp.sqrt(safe_divide(jnp.sum(dp**2), jnp.sum(p**2)))
    jax.debug.print('{}', err)
    return (p, d, r0, r, A1, A2, err)

  x = (p, d, r0, r, A1, A2, err)
  x = jax.lax.while_loop(condition, body, x)
  p = x[0]
  return p

