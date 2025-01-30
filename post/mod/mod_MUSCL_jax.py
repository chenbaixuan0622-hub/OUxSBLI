import jax
import jax.numpy as jnp
from functools import partial


class MUSCL:
  
  @staticmethod
  @jax.jit
  def __minmod2(x, y):
    sgn = jnp.sign(x)
    return sgn * jnp.maximum(jnp.minimum(jnp.abs(x), sgn * y), 0.e0)
 

  @staticmethod
  @jax.jit
  def __minmod3(x, y, z):
    sgn = jnp.sign(x)
    return sgn * jnp.maximum(jnp.minimum(jnp.abs(x), jnp.minimum(sgn * y, sgn * z)), 0.e0)

  @staticmethod
  @jax.jit
  def __d33(d1, d2, d3):
    da = MUSCL.__minmod3(d1, 2.e0 * d2, 2.e0 * d3)
    db = MUSCL.__minmod3(d2, 2.e0 * d1, 2.e0 * d3)
    dc = MUSCL.__minmod3(d3, 2.e0 * d1, 2.e0 * d2)
    return da - 2.e0 * db + dc

  @staticmethod
  @partial(jax.jit, static_argnums=())
  def MUSCL3rd(x):
    k  = 1.e0 / 3.e0
    b  = (3.e0 - k) / (1.e0 - k)
    d  = -x[:-1] + x[1:]
    d1 = MUSCL.__minmod2(d[0], b * d[1])
    d2 = MUSCL.__minmod2(d[1], b * d[0])
    d3 = MUSCL.__minmod2(d[2], b * d[1])
    d4 = MUSCL.__minmod2(d[1], b * d[2])
    xl = x[1] + 0.25e0 * ((1.e0 - k) * d1 + (1.e0 + k) * d2)
    xr = x[2] - 0.25e0 * ((1.e0 - k) * d3 + (1.e0 + k) * d4)
    return xl, xr


  @staticmethod
  @partial(jax.jit, static_argnums=())
  def MUSCL4th(x):
    d  = -x[:-1] + x[1:]
    d1 = d[1] - MUSCL.__d33(d[0], d[1], d[2]) / 6.e0
    d2 = d[2] - MUSCL.__d33(d[1], d[2], d[3]) / 6.e0
    d3 = d[3] - MUSCL.__d33(d[2], d[3], d[4]) / 6.e0
    dl = MUSCL.__minmod2(d1, 4.e0 * d2)
    dr = MUSCL.__minmod2(d2, 4.e0 * d1)
    xl = x[2] + (dl + 2.e0 * dr) / 6.e0
    dl = MUSCL.__minmod2(d2, 4.e0 * d3)
    dr = MUSCL.__minmod2(d3, 4.e0 * d2)
    xr = x[3] - (dr + 2.e0 * dl) / 6.e0
    return xl, xr


  @staticmethod
  @partial(jax.jit, static_argnums=(0, 1, 2))
  def MUSCL4th_x(nx, ny, nz, u):
    us = jnp.zeros((nz-2,ny-2,nx-1))
    ul, ur = MUSCL.MUSCL4th(u)
    return ul, ur


  '''
  @staticmethod
  @partial(jax.jit, static_argnums=(0, 1, 2, 3, 4))
  def __cond3(cond1, cond2, fun1, fun2, fun3, x):
    def fun23(x):
      return jax.lax.cond(cond2, fun2, fun3, x)
    return jax.lax.cond(cond1, fun1, fun23, x)


  @staticmethod
  @partial(jax.jit, static_argnums=(0))
  def MUSCL_x(nx, u, i, j, k):
    @jax.jit
    def fun1(x):
      u, i, j, k = x
      ul, ur = MUSCL.MUSCL4th(u[k,j,i-2:i+3])
      return (ul, ur)

    @jax.jit
    def fun2(x):
      u, i, j, k = x
      ul, ur = MUSCL.MUSCL3rd(u[k,j,i-1:i+2])
      return (ul, ur)

    @jax.jit
    def fun3(x):
      u, i, j, k = x
      um = 0.5e0 * (u[k,j,i] + u[k,j,i+1]) 
      return (um, um)

    x = MUSCL.__cond3(2 <= i and i <= nx-3, 1 <= i and i <= nx-2, fun1, fun2, fun3, (u, i, j, k))
    ul, ur = x
    return ul, ur


  @staticmethod
  @partial(jax.jit, static_argnums=(0))
  def MUSCL_y(ny, v, i, j, k):
    def fun1(x):
      v, i, j, k = x
      vl, vr = MUSCL.MUSCL4th(v[k,j-2:j+3,i])
      return (vl, vr)

    def fun2(x):
      v, i, j, k = x
      vl, vr = MUSCL.MUSCL3rd(v[k,j-1:j+1,i])
      return (vl, vr)

    def fun3(x):
      w, i, j, k = x
      vm = 0.5e0 * (v[k,j,i] + v[k,j+1,i])
      return (vm, vm)

    x = MUSCL.__cond3(2 <= j and j <= ny-3, 1 <= j and j <= ny-2, fun1, fun2, fun3, (v, i, j, k))
    vl, vr = x
    return vl, vr


  @staticmethod
  @partial(jax.jit, static_argnums=(0))
  def MUSCL_z(nz, w, i, j, k):
    def fun1(x):
      w, i, j, k = x
      wl, wr = MUSCL.MUSCL4th(w[k-2:k+3,j,i])
      return (wl, wr)

    def fun2(x):
      w, i, j, k = x
      wl, wr = MUSCL.MUSCL3rd(w[k-1:k+1,j,i])
      return (wl, wr)

    def fun3(x):
      w, i, j, k = x
      wm = 0.5e0 * (w[k,j,i] + w[k+1,j,i])
      return (wm, wm)

    x = MUSCL.__cond3(2 <= k and k <= nz-3, 1 <= k and k <= nz-2, fun1, fun2, fun3, (w, i, j, k))
    wl, wr = x
    return wl, wr


  @staticmethod
  @partial(jax.jit, static_argnums=(0, 1, 2))
  def u_staggered(nx, ny, nz, u):
    us = jnp.zeros((nz-2,ny-2,nx-1))
    def loop_k(k, x):
      nx, u, us = x
      def loop_j(j, x):
        nx, u, us, k = x
        def loop_i(i, x):
          nx, u, us, k, j = x
          ul, ur = 0.5e0, 0.5e0 #MUSCL.__MUSCL_x(nx, u, i, j, k)
          us = us.at[k-1,j-1,i].set(0.5e0 * (ul + ur))
          return (nx, u, us, k, j)
        us = jax.lax.fori_loop(0, nx-1, loop_i, (nx, u, us, k, j))
        return (nx, u, us, k)
      us = jax.lax.fori_loop(1, ny-1, loop_j, (nx, u, us, k))
      return (nx, u, us)
    us = jax.lax.fori_loop(1, nz-1, loop_k, (nx, u, us))
    return us


  @staticmethod
  @partial(jax.jit, static_argnums=(0, 1, 2))
  def v_staggered(nx, ny, nz, v):
    vs = jnp.zeros((nz-2,ny-1,nx-2), dtype=jnp.float32)
    def loop_k(k, x):
      ny, v, vs = x
      def loop_j(j, x):
        ny, v, vs, k = x
        def loop_i(i, x):
          ny, v, vs, k, j = x
          vl, vr = MUSCL.__MUSCL_y(ny, v, i, j, k)
          vs = vs.at[k-1,j,i-1].set(0.5e0 * (vl + vr))
          return vs
        vs = jax.lax.fori_loop(1, nx-1, loop_i, (ny, v, vs, k, j))
        return vs
      vs = jax.lax.fori_loop(0, ny-1, loop_j, (ny, v, vs, k))
      return vs
    vs = jax.lax.fori_loop(1, nz-1, loop_k, (ny, v, vs))
    return vs


  @staticmethod
  @partial(jax.jit, static_argnums=(0, 1, 2))
  def w_staggered(nx, ny, nz, w):
    ws = jnp.zeros((nz-1,ny-2,nx-2), dtype=jnp.float32)
    def loop_k(k, x):
      nz, w, ws = x
      def loop_j(j, x):
        nz, w, ws, k = x
        def loop_i(i, x):
          nz, w, ws, k, j = x
          wl, wr = MUSCL.__MUSCL_z(nz, w, i, j, k)
          ws = ws.at[k,j-1,i-1].set(0.5e0 * (wl + wr))
          return ws
        ws = jax.lax.fori_loop(1, nx-1, loop_i, (nz, w, ws, k, j))
        return ws
      ws = jax.lax.fori_loop(1, ny-1, loop_j, (nz, w, ws, k))
      return ws
    ws = jax.lax.fori_loop(0, nz-1, loop_k, (nz, w, ws))
    return ws
  '''
