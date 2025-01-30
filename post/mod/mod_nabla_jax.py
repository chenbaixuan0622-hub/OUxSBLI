import jax
import jax.numpy as jnp
from functools import partial
from dataclasses import dataclass
from mod_MUSCL_jax import MUSCL


@jax.tree_util.register_pytree_node_class
@dataclass(frozen=True)
class Vector(MUSCL):
  nx: int
  ny: int
  nz: int
  dx: float
  dy: float
  dz: float

  def tree_flatten(self):
    return (self.nx, self.ny, self.nz, self.dx, self.dy, self.dz), ()

  @classmethod
  def tree_unflatten(cls, static, dynamic):
    return cls(*dynamic)

  @jax.jit
  def def_print(self):
    jax.debug.print('{}', self.nx)

  @jax.jit
  def divergence(self, U):
    u   = jnp.zeros((nz-2,ny-2,nx-1), dtype=jnp.float32)
    v   = jnp.zeros((nz-2,ny-1,nx-2), dtype=jnp.float32)
    w   = jnp.zeros((nz-1,ny-2,nx-2), dtype=jnp.float32)
    div = jnp.zeros((nz,  ny,  nx),   dtype=jnp.float32)

    @jax.jit
    def loop_k(k, u):
      @jax.jit
      def loop_j(j, u):
        @jax.jit
        def loop_i(i, u):
          u = u.at[k,j,i].set(1.e0)
          return u
          u = jax.lax.fori_loop(1, nx-1, loop_i, u)
        return u
        u = jax.lax.fori_loop(1, ny-1, loop_j, u)
      return u
    u = jax.lax.fori_loop(1, nz-1, loop_k, u)
    
    return div


nx = 64
ny = 64
nz = 64
dx = 1.e0 / nx
dy = 1.e0 / ny
dz = 1.e0 / nz
U  = jnp.zeros((nz, ny, nx), dtype=jnp.float32)
nabla = Vector(nx, ny, nz, dx, dy, dz)
nabla.def_print()
nabla.divergence(U)

