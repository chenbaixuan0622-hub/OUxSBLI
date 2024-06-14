import numpy as np
from numba import jit

@jit(nopython=True, cache=True, fastmath=True)
def unit_vector(r1,r2):
  R = np.sqrt(r1**2 + r2**2)
  r = np.array([r1 / R, r2 / R])
  return r

@jit(nopython=True, cache=True, fastmath=True)
def orthogonal_unit_vector(x1,y1,x2,y2):
  r = np.array([x2 - x1, y2 - y1])
  R = np.sqrt(r[0]**2 + r[1]**2)
  n = np.array([-r[1] / R, r[0] / R])
  return n

@jit(nopython=True, cache=True, fastmath=True)
def velocity_corr(u1,u2):
  n = 1.e0
  u1u2 = 0.e0
  u12 = 0.e0
  u22 = 0.e0
  for i in range(u1.size):
    u1u2 = ((n - 1.e0) * u1u2 + u1[i] * u2[i]) / n
    u12  = ((n - 1.e0) * u12  + u1[i]**2)      / n
    u22  = ((n - 1.e0) * u22  + u2[i]**2)      / n
    n += 1.e0
  return u1u2 / (np.sqrt(u12) * np.sqrt(u22))

# longitudinal velocity correlation function
@jit(nopython=True, cache=True, fastmath=True)
def longitudinal_corr(x1,y1,u1,v1,x2,y2,u2,v2):
  r = unit_vector(x2 - x1, y2 - y1)
  U1 = r[0] * u1 + r[1] * v1
  U2 = r[0] * u2 + r[1] * v2
  R11 = velocity_corr(U1,U2)
  return R11

# lateral velociy correlation function
@jit(nopython=True, cache=True, fastmath=True)
def lateral_corr(x1,y1,u1,v1,x2,y2,u2,v2):
  n = orthogonal_unit_vector(x1,y1,x2,y2)
  U1 = n[0] * u1 + n[1] * v1
  U2 = n[0] * u2 + n[1] * v2 
  R22 = velocity_corr(U1,U2)
  return R22

# longitudinal and lateral correlation
#@jit(nopython=True, cache=True, fastmath=True)
def longitudinal_and_lateral_corr(span,length,Nx1,Ny1,x,y,us,vs):
  R11 = np.zeros(length-1)
  R22 = np.zeros(length-1)
  R11_span = np.zeros((span,length-1))
  R22_span = np.zeros((span,length-1))
  for i in range(length-1):
    I = Nx1 + i + 1
    for j in range(span):
      print(x[Nx1],x[I])
      R11_span[j,i] = longitudinal_corr(x[Nx1],y[Ny1],us[:,j,0],vs[:,j,0],x[I],y[Ny1],us[:,j,i+1],vs[:,j,i+1])
      R22_span[j,i] = lateral_corr(x[Nx1],y[Ny1],us[:,j,0],vs[:,j,0],x[I],y[Ny1],us[:,j,i+1],vs[:,j,i+1])
    R11[i] = np.mean(R11_span[:,i])
    R22[i] = np.mean(R22_span[:,i])
  return R11, R22

# integral scale
@jit(nopython=True, cache=True, fastmath=True)
def integral_scale(x,R11,R22):
  L11 = 0.e0
  L22 = 0.e0
  for i in range(len(x)-1):
    dx = -x[i] + x[i+1]
    L11 += 0.5e0 * (R11[i] + R11[i+1]) * dx 
    L22 += 0.5e0 * (R22[i] + R22[i+1]) * dx
  return L11, L22

