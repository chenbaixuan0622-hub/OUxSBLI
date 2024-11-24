import numpy as np
from numba import njit
from scipy.interpolate import interp1d
from mod.mod_MUSCL import MUSCL3rd, MUSCL4th


@njit(cache=True, fastmath=True, nogil=True)
def divergence(x, y, z, u, v, w):
  Nx = len(x)
  Ny = len(y)
  Nz = len(z)
  div = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  if Ny >= 3:
    for k in range(1,Nz-1):
      for j in range(1,Ny-1):
        for i in range(1,Nx-1):
          dx = -x[i-1] + x[i+1]
          dy = -y[j-1] + y[j+1]
          dz = -z[k-1] + z[k+1]
          div[k,j,i] = (-u[k,j,i-1] + u[k,j,i+1]) / dx + \
                       (-v[k,j-1,i] + v[k,j+1,i]) / dy + \
                       (-w[k-1,j,i] + w[k+1,j,i]) / dz
    div[:,:,0]  = div[:,:,1]
    div[:,:,-1] = div[:,:,-2]
    div[:,0,:]  = div[:,1,:]
    div[:,-1,:] = div[:,-2,:]
    div[0,:,:]  = div[1,:,:]
    div[-1,:,:] = div[-2,:,:]
  else:
    for k in range(1,Nz-1):
      for j in range(Ny-1):
        for i in range(1,Nx-1):
          dx = -x[i-1] + x[i+1]
          dy = -y[j] + y[j+1]
          dz = -z[k-1] + z[k+1]
          div[k,j,i] = (-u[k,j,i-1] + u[k,j,i+1]) / dx + \
                       (-v[k,j,i]   + v[k,j+1,i]) / dy + \
                       (-w[k-1,j,i] + w[k+1,j,i]) / dz
    div[:,:,0]  = div[:,:,1]
    div[:,:,-1] = div[:,:,-2]
    div[0,:,:]  = div[1,:,:]
    div[-1,:,:] = div[-2,:,:]
    div[:,1,:]  = div[:,0,:]
  return div


@njit(cache=True, fastmath=True, nogil=True)
def Ducros(x, y, z, u, v, w):
  Nx = len(x)
  Ny = len(y)
  Nz = len(z)
  eps = 1.e-16
  sensor = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  for k in range(1,Nz-1):
    for j in range(1,Ny-1):
      for i in range(1,Nx-1):
        dx = -x[i-1] + x[i+1]
        dy = -y[j-1] + y[j+1]
        dz = -z[k-1] + z[k+1]
        dudx = (-u[k,j,i-1] + u[k,j,i+1]) / dx
        dvdy = (-v[k,j-1,i] + v[k,j+1,i]) / dy
        dwdz = (-w[k-1,j,i] + w[k+1,j,i]) / dz
        dwdy = (-w[k,j-1,i] + w[k,j+1,i]) / dy
        dvdz = (-v[k-1,j,i] + v[k+1,j,i]) / dz
        dudz = (-u[k-1,j,i] + u[k+1,j,i]) / dz
        dwdx = (-w[k,j,i-1] + w[k,j,i+1]) / dx
        dvdx = (-v[k,j,i-1] + v[k,j,i+1]) / dx
        dudy = (-u[k,j-1,i] + u[k,j+1,i]) / dy
        div  = dudx + dvdy + dwdz
        rot1 = dwdy - dvdz 
        rot2 = dudz - dwdx
        rot3 = dvdx - dudy
        sensor[k,j,i] = div**2 / (rot1**2 + rot2**2 + rot3**2 + eps)
  sensor[:,:,0]  = sensor[:,:,1]
  sensor[:,:,-1] = sensor[:,:,-2]
  sensor[:,0,:]  = sensor[:,1,:]
  sensor[:,-1,:] = sensor[:,-2,:]
  sensor[0,:,:]  = sensor[1,:,:]
  sensor[-1,:,:] = sensor[-2,:,:]
  return sensor


@njit(cache=True, fastmath=True, nogil=True)
def interp(x, a):
  Nx = len(x)
  xi = np.linspace(x[0], x[-1], 2*Nx-1)
  ai = np.zeros(2*Nx-1, dtype=np.float32)
  return xi, ai


@njit(cache=True, fastmath=True, nogil=True)
def MUSCL(xi, x, a):
  Nx = len(x)
  ai = np.zeros_like(xi)
  ai[0:-1:2] = a
  for i in range(Nx-1):
    if 2 <= i and i <= Nx-4:
      a6 = a[i-2:i+3]
      al, ar = MUSCL4th(a6)
    elif 1 <= i and i <= Nx-3:
      a4 = a[i-1:i+2]
      al, ar = MUSCL3rd(a4)
    elif i == 0:
      al, ar = a[0], a[0]
    else:
      al, ar = a[-1], a[-1]
    ai[2*i+1] = 0.5e0 * (al + ar)
  return ai


def reflected_shock(x, y, z, ny, u, v, w):
  Nz  = len(z)
  ls  = np.zeros(Nz, dtype=np.float32)
  lsi = np.zeros(Nz, dtype=np.float32)
  div = divergence(x, y, z, u, v, w)
  Nx  = len(x)
  xi  = np.linspace(x[0], x[-1], 2*Nx-1)
  for k in range(Nz):
    ls[k]  = x[np.argmin(div[k,ny,:])]
    #divi   = MUSCL(xi, x, div[k,ny,:])
    f      = interp1d(x, div[k,ny,:], kind='linear')
    divi   = f(xi)
    #divi   = np.interp(xi, x, div[k,ny,:])
    lsi[k] = xi[np.argmin(divi)]
  return ls, lsi

