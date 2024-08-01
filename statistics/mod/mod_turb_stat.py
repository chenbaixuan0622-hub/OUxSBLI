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

@jit(nopython=True, cache=True, fastmath=True)
def blt(y,u):
  d  = 0.e0
  u0 = np.mean(u[:,-1,:])
  for j in range(len(y)):
    if np.mean(u[:,j,:]) >= 0.99e0 * u0:
      d = y[j] - (-y[j-1] + y[j]) * (np.mean(u[:,j,:]) - 0.99e0 * u0) \
      / (-np.mean(u[:,j-1,:]) + np.mean(u[:,j,:]) + 1.e-20)
  return d

@jit(nopython=True, cache=True, fastmath=True)
def Sutherland(T):
  mu0  = 1.716e-5
  T0   = 273.2e0
  S    = 111.e0
  nu   = mu0 * ((T0 + S) / (np.mean(T) + S)) * (np.mean(T) / T0)**1.5
  return nu

# calc yplus, uplus
@jit(nopython=True, cache=True, fastmath=True)
def non_dim_tbl(Q,x,y,z):
  # Q[rho,u,v,w,p]
  R    = 287.03e0
  rhow = np.mean(Q[0,:,0,:])
  nu   = Sutherland(Q[4,:,0,:] / (R * Q[0,:,0,:])) / rhow
  dudy = np.mean(-Q[1,:,0,:] + Q[1,:,1,:]) / (-y[0] + y[1])
  tw   = rhow * nu * dudy
  ut   = np.sqrt(tw / rhow)
  nx   = len(x)
  ny   = len(y)
  nz   = len(z)
  Qp   = np.zeros((5,nz,ny,nx), dtype=np.float32)
  Qvd  = np.zeros((3,nz,ny,nx), dtype=np.float32)
  xp   = ut * x[:] / nu
  yp   = ut * y[:] / nu
  zp   = ut * z[:] / nu
  Qp[0,:,0,:] = 1.e0
  Qp[1,:,0,:] = Q[1,:,0,:] / ut
  Qp[2,:,0,:] = Q[2,:,0,:] / ut
  Qp[3,:,0,:] = Q[3,:,0,:] / ut
  Qp[4,:,0,:] = Q[4,:,0,:] / (rhow * ut**2)
  for k in range(nz):
    for j in range(1,ny):
      for i in range(nx):
        # van Driest transformation
        uvd = Q[1,k,j-1,i] + np.sqrt(Q[0,k,j,i] / rhow) * (-Q[1,k,j-1,i] + Q[1,k,j,i])
        vvd = Q[2,k,j-1,i] + np.sqrt(Q[0,k,j,i] / rhow) * (-Q[2,k,j-1,i] + Q[2,k,j,i])
        wvd = Q[3,k,j-1,i] + np.sqrt(Q[0,k,j,i] / rhow) * (-Q[3,k,j-1,i] + Q[3,k,j,i])
        Qp[0,k,j,i]  = Q[0,k,j,i] / rhow
        Qp[1,k,j,i]  = Q[1,k,j,i] / ut
        Qp[2,k,j,i]  = Q[2,k,j,i] / ut
        Qp[3,k,j,i]  = Q[3,k,j,i] / ut
        Qp[4,k,j,i]  = Q[4,k,j,i] / (rhow * ut**2)
        Qvd[0,k,j,i] = uvd / ut
        Qvd[1,k,j,i] = vvd / ut
        Qvd[2,k,j,i] = wvd / ut
  return xp, yp, zp, Qp, Qvd

# calc point-wise turbulent kinetic energy
@jit(nopython=True, cache=True, fastmath=True)
def TKE(u,v,w,U,V,W,tke):
  u2 = 0.e0
  v2 = 0.e0
  w2 = 0.e0
  for t in range(len(u)):
    u2 += (u[t] - U)**2
    v2 += (v[t] - V)**2
    w2 += (w[t] - W)**2
  tke += 0.5e0 * (u2 + v2 + w2)
  return tke

