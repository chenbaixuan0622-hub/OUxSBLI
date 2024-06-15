import jax
import jax.numpy as jnp
from jax import lax
from functools import partial

@jax.jit
def unit_vector(r1,r2):
  R = jnp.sqrt(r1**2 + r2**2)
  r = jnp.array([r1 / R, r2 / R])
  return r

@jax.jit
def orthogonal_unit_vector(x1,y1,x2,y2):
  r = jnp.array([x2 - x1, y2 - y1])
  R = jnp.sqrt(r[0]**2 + r[1]**2)
  n = jnp.array([-r[1] / R, r[0] / R])
  return n

@jax.jit
def body_velocity_corr(i,xs):
  n, u1u2, u12, u22 = xs
  u1u2 = ((n - 1.e0) * u1u2 + u1[i] * u2[i]) / n
  u12  = ((n - 1.e0) * u12  + u1[i]**2)      / n
  u22  = ((n - 1.e0) * u22  + u2[i]**2)      / n
  n += 1.e0
  return (n, u1u2, u12, u22)

@partial(jax.jit, static_argnames=['length'])
def velocity_corr(length,u1,u2):
  n = 1.e0
  u1u2 = 0.e0
  u12 = 0.e0
  u22 = 0.e0
  xs = (n, u1u2, u12, u22)
  lax.fori_loop(0,length,body_velocity_corr,xs)
  return u1u2 / (jnp.sqrt(u12) * jnp.sqrt(u22))

# longitudinal velocity correlation function
@partial(jax.jit, static_argnames=['length'])
def longitudinal_corr(length,x1,y1,u1,v1,x2,y2,u2,v2):
  r = unit_vector(x2 - x1, y2 - y1)
  U1 = r[0] * u1 + r[1] * v1
  U2 = r[0] * u2 + r[1] * v2
  R11 = velocity_corr(length,U1,U2)
  return R11

# lateral velociy correlation function
@partial(jax.jit, static_argnames=['length'])
def lateral_corr(length,x1,y1,u1,v1,x2,y2,u2,v2):
  n = orthogonal_unit_vector(x1,y1,x2,y2)
  U1 = n[0] * u1 + n[1] * v1
  U2 = n[0] * u2 + n[1] * v2 
  R22 = velocity_corr(length,U1,U2)
  return R22

# longitudinal and lateral correlation
@partial(jax.jit, static_argnames=['span', 'length'])
def longitudinal_and_lateral_corr(span,length,Nx1,Ny1,x,y,us,vs):
  '''
  R11 = jnp.zeros(length-1)
  R22 = jnp.zeros(length-1)
  R11_span = jnp.zeros((span,length-1))
  R22_span = jnp.zeros((span,length-1))
  for i in range(length-1):
    I = Nx1 + i + 1
    for j in range(span):
      R11_span[j,i] = longitudinal_corr(length,x[Nx1],y[Ny1],us[:,j,0],vs[:,j,0],x[I],y[Ny1],us[:,j,i+1],vs[:,j,i+1])
      R22_span[j,i] = lateral_corr(length,x[Nx1],y[Ny1],us[:,j,0],vs[:,j,0],x[I],y[Ny1],us[:,j,i+1],vs[:,j,i+1])
    R11[i] = jnp.mean(R11_span[:,i])
    R22[i] = jnp.mean(R22_span[:,i])
  '''
  def calc_corr(i):
    I = Nx1 + i + 1
    R11_span = jax.vmap(lambda j: \
    longitudinal_corr(length,x[Nx1],y[Ny1],us[:,j,0],vs[:,j,0],x[I],y[Ny1],us[:,j,i+1],vs[:,j,i+1]))(jnp.arange(span))
    R22_span = jax.vmap(lambda j: \
    lateral_corr(length,x[Nx1],y[Ny1],us[:,j,0],vs[:,j,0],x[I],y[Ny1],us[:,j,i+1],vs[:,j,i+1]))(jnp.arange(span))
    return jnp.mean(R11_span), jnp.mean(R22_span)

  R11, R22 = jax.vmap(calc_corr)(jnp.arange(length - 1))
  return R11, R22

# integral scale
@jax.jit
def body_integral_scale(i,xs):
  x, L11, L22 = xs
  dx = -x[i] + x[i+1]
  L11 += 0.5e0 * (R11[i] + R11[i+1]) * dx
  L22 += 0.5e0 * (R22[i] + R22[i+1]) * dx
  return (x, L11, L22)

@jax.jit
def integral_scale(x,R11,R22):
  L11 = 0.e0
  L22 = 0.e0
  lax.fori_loop(0,len(x)-1,body_integral_scale,xs)
  return L11, L22

