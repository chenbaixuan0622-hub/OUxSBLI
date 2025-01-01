import numpy as np
from numba import njit


@njit(cache=True, nogil=True)
def Neumann_bc2D(phi):
  phi[1:-1,0]  = phi[1:-1,1]
  phi[1:-1,-1] = phi[1:-1,-2]
  phi[0,:]     = phi[1,:]
  phi[-1,:]    = phi[-2,:]
  return phi


@njit(cache=True, nogil=True)
def Neumann_bc3D(phi):
  phi[1:-1,1:-1,0]  = phi[1:-1,1:-1,1]
  phi[1:-1,1:-1,-1] = phi[1:-1,1:-1,-2]
  phi[1:-1,0,:]     = phi[1:-1,1,:]
  phi[1:-1,-1,:]    = phi[1:-1,-2,:]
  phi[0,:,:]        = phi[-2,:,:]
  phi[:-1,:,:]      = phi[1,:,:]
  return phi


@njit(cache=True, nogil=True)
def df(f, h1, h2):
  return (-h2**2 * f[0] + (h2**2 - h1**2) * f[1] + h1**2 * f[2]) / (h1 * h2 * (h1 + h2))


class Scalar:
  def __init__(self, p, x, y, z=None):
    self.p  = p
    self.x  = x
    self.y  = y
    self.z  = z
    self.dx = -x[:-1] + x[1:]
    self.dy = -y[:-1] + y[1:]
    if self.z is not None:
      self.dz = -z[:-1] + z[1:]
      self.nz = len(z)

  @staticmethod
  @njit(cache=True, nogil=True)
  def gradient2D(p, x, y, dx, dy):
    grad = np.zeros((2,len(y),len(x)))
    for j in range(1,len(y)-1):
      for i in range(1,len(x)-1):
        grad[0,j,i] = df(p[j,i-1:i+1], dx[i-1], dx[i])
        grad[1,j,i] = df(p[j-1:j+1,i], dy[j-1], dy[j])
    return grad

  @staticmethod
  @njit(cache=True, nogil=True)
  def gradient3D(p, x, y, z, dx, dy, dz):
    grad = np.zeros((3,len(z),len(y),len(x)))
    for k in range(1,len(z)-1):
      for j in range(1,len(y)-1):
        for i in range(1,len(x)-1):
          grad[0,k,j,i] = df(p[k,j,i-1:i+1], dx[i-1], dx[i])
          grad[1,k,j,i] = df(p[k,j-1:j+1,i], dy[j-1], dy[j])
          grad[2,k,j,i] = df(p[k-1:k+1,j,i], dz[k-1], dz[k])
    return grad

  def gradient(self):
    if self.z is None:
      return self.gradient2D(self.p, self.x, self.y, self.dx, self.dy)
    else:
      return self.gradient3D(self.p, self.x, self.y, self.z, self.dx, self.dy, self.dz)


class Vector:
  def __init__(self, u, x, y, z=None):
    self.u  = u
    self.x  = x
    self.y  = y
    self.z  = z
    self.dx = -x[:-1] + x[1:]
    self.dy = -y[:-1] + y[1:]
    self.nx = len(x)
    self.ny = len(y)
    if self.z is not None:
      self.dz = -z[:-1] + z[1:]
      self.nz = len(z)
  
  @staticmethod
  @njit(cache=True, nogil=True)
  def divergence2D(u, nx, ny, dx, dy):
    div = np.zeros((ny,nx), dtype=np.float32)
    for j in range(1,ny-1):
      for i in range(1,nx-1):
        div[j,i] = df(u[0,j,i-1:i+1], dx[i-1], dx[i]) \
                 + df(u[1,j-1:j+1,i], dy[j-1], dy[j])
    return div

  @staticmethod
  @njit(cache=True, nogil=True)
  def divergence3D(u, nx, ny, nz, dx, dy, dz):
    div = np.zeros((nz,ny,nx), dtype=np.float32)
    for k in range(1,nz-1):
      for j in range(1,ny-1):
        for i in range(1,nx-1):
          div[k,j,i] = df(u[0,k,j,i-1:i+1], dx[i-1], dx[i]) \
                     + df(u[1,k,j-1:j+1,i], dy[j-1], dy[j]) \
                     + df(u[2,k-1:k+1,j,i], dz[k-1], dz[k])
    return div

  @staticmethod
  @njit(cache=True, nogil=True)
  def rotation2D(u, nx, ny, dx, dy):
    rot = np.zeros((ny,nx), dtype=np.float32)
    for j in range(1,ny-1):
      for i in range(1,nx-1):
        rot[j,i] = df(u[1,j,i-1:i+1], dx[i-1], dx[i]) - df(u[0,j-1:j+1,i], dy[j-1], dy[j])
    return rot

  @staticmethod
  @njit(cache=True, nogil=True)
  def rotation3D(u, nx, ny, nz, dx, dy, dz):
    rot = np.zeros((3,nz,ny,nx), dtype=np.float32)
    for k in range(1,nz-1):
      for j in range(1,ny-1):
        for i in range(1,nx-1):
          rot[0,k,j,i] = df(u[2,k,j-1:j+1,i], dy[j-1], dy[j]) - df(u[1,k-1:k+1,j,i], dz[k-1], dz[k])
          rot[1,k,j,i] = df(u[0,k-1:k+1,j,i], dz[k-1], dz[k]) - df(u[2,k,j,i-1:i+1], dx[i-1], dx[i])
          rot[2,k,j,i] = df(u[1,k,j,i-1:i+1], dx[i-1], dx[i]) - df(u[0,k,j-1:j+1,i], dy[j-1], dy[j])
    return rot

  def divergence(self):
    if self.z is None:
      return self.divergence2D(self.u, self.nx, self.ny, self.dx, self.dy)
    else:
      return self.divergence3D(self.u, self.nx, self.ny, self.nz, self.dx, self.dy, self.dz)

  def rotation(self):
    if self.z is None:
      return self.rotation2D(self.u, self.nx, self.ny, self.dx, self.dy)
    else:
      return self.rotation3D(self.u, self.nx, self.ny, self.nz, self.dx, self.dy, self.dz)

