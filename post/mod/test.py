import jax
import jax.numpy as jnp
from mod_MUSCL_jax import MUSCL


a6 = jnp.ones(6)
al, ar = MUSCL.MUSCL4th(a6)
print(al, ar)

a4 = jnp.ones(4)
al, ar = MUSCL.MUSCL3rd(a4)
print(al, ar)

'''
x = 3.e0

@jax.jit
def fun1(x):
  return x

@jax.jit
def fun2(x):
  return 2.e0 * x

@jax.jit
def fun3(x):
  return 3.e0 * x


y = MUSCL.cond3(3 <= x and x <= 10, 2 <= x and x <= 11, fun1, fun2, fun3, x)

print(y)
'''

nx = 31
ny = 32
nz = 33
u  = jnp.ones((nz,   ny,   nx),   dtype=jnp.float32)
us = jnp.ones((nz-2, ny-2, nx-1), dtype=jnp.float32)

def loop(i, x):
  u, us = x
  u4 = u[1:-1,1:-1,i-1:i+2]#.transpose((2,1,0))
  us = us.at[:,:,i].set(MUSCL.MUSCL3rd(u4))
  return (u, us)

x = jax.lax.fori_loop(1, nx-2, loop, (u, us))
u, us = x
#ul, ur = MUSCL.MUSCL3rd(u)
print(jnp.shape(us))

