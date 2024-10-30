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

@jit(nopython=True, cache=True, fastmath=True, nogil=True)
def Sutherland(T):
  mu0 = 1.716e-5
  T0  = 273.2e0
  S   = 111.e0
  mu  = mu0 * ((T0 + S) / (T + S)) * (T / T0)**1.5
  return mu

@jit(nopython=True, cache=True, fastmath=True, nogil=True)
def non_dim_tbl(Q, x, y, z):
  # Q[rho,u,v,w,p]
  R   = 287.03e0
  nx  = len(x)
  ny  = len(y)
  nz  = len(z)
  nu  = np.zeros((nz,nx), dtype=np.float32)
  tw  = np.zeros((nz,nx), dtype=np.float32)
  ut  = np.zeros((nz,nx), dtype=np.float32)
  uvd = np.zeros(ny, dtype=np.float32)
  twm = 0.e0
  utm = 0.e0
  itr = 0.e0
  for k in range(nz):
    for i in range(nx):
      rhow    = Q[0,k,0,i]
      mu      = Sutherland(Q[4,k,0,i] / (R * rhow))
      nu[k,i] = mu / rhow
      dudy    = (-Q[1,k,0,i] + Q[1,k,1,i]) / (-y[0] + y[1])
      tw[k,i] = mu * dudy
      if dudy >= 0.e0:
        ut[k,i] = np.sqrt(tw[k,i] / rhow)
        twm += tw[k,i]
        utm += ut[k,i]
        itr += 1.e0
      else:
        ut[k,i] = -np.sqrt(-tw[k,i] / rhow)
  twm /= itr
  utm /= itr
  for j in range(1,ny-1):
    # van Driest transformation
    uvd[j] = uvd[j-1] + np.sqrt(np.mean(Q[0,:,j,:] / Q[0,:,0,:])) * 0.5e0 * np.mean(-Q[1,:,j-1,:] + Q[1,:,j+1,:])
  uvd[ny-1] = uvd[ny-2] + np.sqrt(np.mean(Q[0,:,ny-1,:] / Q[0,:,0,:])) * np.mean(-Q[1,:,ny-2,:] + Q[1,:,ny-1,:])
  yp  = utm * y[:] / np.mean(nu)
  up  = uvd[:] / utm
  return yp, tw, ut, twm, utm, up

@jit(nopython=True, cache=True, fastmath=True, nogil=True)
def tau_2d(Q, x, y, z):
  # Q[rho,u,v,w,p]
  R    = 287.03e0
  nx   = len(x)
  ny   = len(y)
  nz   = len(z)
  tw   = np.zeros((nz, nx), dtype=np.float32)
  ut   = np.zeros((nz, nx), dtype=np.float32)
  for k in range(nz):
    for i in range(nx):
      rhow    = Q[0,k,0,i]
      mu      = Sutherland(Q[4,k,0,i] / (R * Q[0,k,0,i]))
      dudy    = (-Q[1,k,0,i] + Q[1,k,2,i]) / (-y[0] + y[2])
      tw[k,i] = mu * np.abs(dudy)
      ut[k,i] = np.sqrt(tw[k,i] / rhow)
  return tw, ut

# calc point-wise turbulent kinetic energy
@jit(nopython=True, cache=True, fastmath=True, nogil=True)
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

@jit(nopython=True, cache=True, fastmath=True, nogil=True)
def tau(Nx, Ny, Nz, dx, y, dz, Rgas, rho, u, v, w, p):
  mu  = 1.716e-5 * ((273.2e0 + 111.e0) / (p / (rho * Rgas) + 111.e0)) * ((p / (rho * Rgas)) / 273.2e0)**1.5
  mux = 0.5e0 * (mu[1:-1,1:-1,:-1] + mu[1:-1,1:-1,1:])
  muy = 0.5e0 * (mu[1:-1,:-1,1:-1] + mu[1:-1,1:,1:-1])
  muz = 0.5e0 * (mu[:-1,1:-1,1:-1] + mu[1:,1:-1,1:-1])
  txx = np.zeros((Nz-2,Ny-2,Nx-1), dtype=np.float32)
  tyx = np.zeros((Nz-2,Ny-2,Nx-1), dtype=np.float32)
  tzx = np.zeros((Nz-2,Ny-2,Nx-1), dtype=np.float32)
  txy = np.zeros((Nz-2,Ny-1,Nx-2), dtype=np.float32)
  tyy = np.zeros((Nz-2,Ny-1,Nx-2), dtype=np.float32)
  tzy = np.zeros((Nz-2,Ny-1,Nx-2), dtype=np.float32)
  txz = np.zeros((Nz-1,Ny-2,Nx-2), dtype=np.float32)
  tyz = np.zeros((Nz-1,Ny-2,Nx-2), dtype=np.float32)
  tzz = np.zeros((Nz-1,Ny-2,Nx-2), dtype=np.float32)
  for k in range(1,Nz-1):
    for j in range(1,Ny-1):
      for i in range(Nx-1):
        uy = 0.50e0 * (-u[k,j-1,i] + u[k,j+1,i] - u[k,j-1,i+1] + u[k,j+1,i+1]) / (-y[j-1] + y[j+1])
        vy = 0.50e0 * (-v[k,j-1,i] + v[k,j+1,i] - v[k,j-1,i+1] + v[k,j+1,i+1]) / (-y[j-1] + y[j+1])
        wz = 0.25e0 * (-w[k-1,j,i] + w[k+1,j,i] - w[k-1,j,i+1] + w[k+1,j,i+1]) / dz
        uz = 0.25e0 * (-u[k-1,j,i] + u[k+1,j,i] - u[k-1,j,i+1] + u[k+1,j,i+1]) / dz
        txx[k-1,j-1,i] = 2.e0 * mux[k-1,j-1,i] * (2.e0 * (-u[k,j,i] + u[k,j,i+1]) / dx - vy - wz) / 3.e0
        tyx[k-1,j-1,i] = mux[k-1,j-1,i] * ((-v[k,j,i] + v[k,j,i+1]) / dx + uy)
        tzx[k-1,j-1,i] = mux[k-1,j-1,i] * ((-w[k,j,i] + w[k,j,i+1]) / dx + uz)
  for k in range(1,Nz-1):
    for j in range(Ny-1):
      for i in range(1,Nx-1):
        vz = 0.25e0 * (-v[k-1,j,i] + v[k+1,j,i] - v[k-1,j+1,i] + v[k+1,j+1,i]) / dz
        wz = 0.25e0 * (-w[k-1,j,i] + w[k+1,j,i] - w[k-1,j+1,i] + w[k+1,j+1,i]) / dz
        ux = 0.25e0 * (-u[k,j,i-1] + u[k,j,i+1] - u[k,j+1,i-1] + u[k,j+1,i+1]) / dx
        vx = 0.25e0 * (-v[k,j,i-1] + v[k,j,i+1] - v[k,j+1,i-1] + v[k,j+1,i+1]) / dx
        txy[k-1,j,i-1] = muy[k-1,j,i-1] * ((-u[k,j,i] + u[k,j+1,i]) / (-y[j] + y[j+1]) + vx)
        tyy[k-1,j,i-1] = 2.e0 * muy[k-1,j,i-1] * (2.e0 * (-v[k,j,i] + v[k,j+1,i]) / (-y[j] + y[j+1]) - wz - ux) / 3.e0
        tzy[k-1,j,i-1] = muy[k-1,j,i-1] * ((-w[k,j,i] + w[k,j+1,i]) / (-y[j] + y[j+1]) + vz)
  for k in range(Nz-1):
    for j in range(1,Ny-1):
      for i in range(1,Nx-1):
        wx = 0.25e0 * (-w[k,j,i-1] + w[k,j,i+1] - w[k+1,j,i-1] + w[k+1,j,i+1]) / dx
        ux = 0.25e0 * (-u[k,j,i-1] + u[k,j,i+1] - u[k+1,j,i-1] + u[k+1,j,i+1]) / dx
        vy = 0.50e0 * (-v[k,j-1,i] + v[k,j+1,i] - v[k+1,j-1,i] + v[k+1,j+1,i]) / (-y[j-1] + y[j+1])
        wy = 0.50e0 * (-w[k,j-1,i] + w[k,j+1,i] - w[k+1,j-1,i] + w[k+1,j+1,i]) / (-y[j-1] + y[j+1])
        txz[k,j-1,i-1] = muz[k,j-1,i-1] * ((-u[k,j,i] + u[k+1,j,i]) / dz + wx)
        tyz[k,j-1,i-1] = muz[k,j-1,i-1] * ((-v[k,j,i] + v[k+1,j,i]) / dz + wy)
        tzz[k,j-1,i-1] = 2.e0 * muz[k,j-1,i-1] * (2.e0 * (-w[k,j,i] + w[k+1,j,i]) / dz - ux - vy) / 3.e0
  return txx, tyx, tzx, txy, tyy, tzy, txz, tyz, tzz

