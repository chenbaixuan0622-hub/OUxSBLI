import numpy as np
from numba import njit
from mod.mod_MUSCL import u_staggered, v_staggered, w_staggered


class bc:
  def __init__(self, v):
    self.dim = v.ndim
    if self.dim == 2:
      self.v = v[np.newaxis,:,:]
    elif self.dim == 3:
      self.v = v
    elif self.dim != 2 and self.dim != 3:
      raise Exception('Invalid dimension')

  def periodic(self, x=None, y=None, z=None):
    v = self.v
    
    if x is True:
      v[:,:,0]  = v[:,:,-2]
      v[:,:,-1] = v[:,:,1]

    if y is True:
      v[:,0,:]  = v[:,-2,:]
      v[:,-1,:] = v[:,1,:]

    if z is True and self.dim == 3:
      v[0,:,:]  = v[-2,:,:]
      v[-1,:,:] = v[1,:,:]

    if self.dim == 3:
      return v
    else:
      return v[0,:,:]

  def Neumann(self, x1=None, x2=None, y1=None, y2=None, z1=None, z2=None):
    v = self.v

    if x1 is True:
      v[:,:,0]  = v[:,:,1]
    if x2 is True:
      v[:,:,-1] = v[:,:,-2]
    if y1 is True:
      v[:,0,:]  = v[:,1,:]
    if y2 is True:
      v[:,-1,:] = v[:,-2,:]
    if z1 is True and self.dim == 3:
      v[0,:,:]  = v[1,:,:]
    if z2 is True and self.dim == 3:
      v[-1,:,:] = v[-2,:,:]

    if self.dim == 3:
      return v
    else:
      return v[0,:,:]

  def Dirichlet(self, x1=None, x2=None, y1=None, y2=None, z1=None, z2=None):
    v = self.v

    if self.dim == 2:
      x1 = x1[np.newaxis,:]
      x2 = x2[np.newaxis,:]
      y1 = y1[np.newaxis,:]
      y2 = y2[np.newaxis,:]

    if x1 is not None:
      v[:,:,0]  = x1
    if x2 is not None:
      v[:,:,-1] = x2
    if y1 is not None:
      v[:,0,:]  = y1
    if y2 is not None:
      v[:,-1,:] = y2
    if z1 is not None and self.dim == 3:
      v[0,:,:]  = z1
    if z2 is not None and self.dim == 3:
      v[-1,:,:] = z2

    if self.dim == 3:
      return v
    else:
      return v[0,:,:]


@njit(cache=True, nogil=True)
def df(f, h1, h2):
  return (-h2**2 * f[0] + (h2**2 - h1**2) * f[1] + h1**2 * f[2]) / (h1 * h2 * (h1 + h2))


class Scalar():
  def __init__(self, p, x, y, z=None):
    self.p  = p
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
  def gradient2D(p, nx, ny, dx, dy):
    grad = np.zeros((2,ny,nx), dtype=np.float32)
    for j in range(1,ny-1):
      for i in range(1,nx-1):
        grad[0,j,i] = df(p[j,i-1:i+1], dx[i-1], dx[i])
        grad[1,j,i] = df(p[j-1:j+1,i], dy[j-1], dy[j])
    return grad

  @staticmethod
  @njit(cache=True, nogil=True)
  def gradient3D(p, nx, ny, nz, dx, dy, dz, periodic_z=False):
    grad = np.zeros((3,nz,ny,nx), dtype=np.float32)
    for k in range(1,nz-1):
      for j in range(1,ny-1):
        for i in range(1,nx-1):
          grad[0,k,j,i] = df(p[k,j,i-1:i+1], dx[i-1], dx[i])
          grad[1,k,j,i] = df(p[k,j-1:j+1,i], dy[j-1], dy[j])
          grad[2,k,j,i] = df(p[k-1:k+1,j,i], dz[k-1], dz[k])
    if periodic_z is True:
      for j in range(1,ny-1):
        for i in range(1,nx-1):
          grad[0,0,j,i]  = df(p[0,j,i-1:i+1], dx[i-1], dx[i])
          grad[1,0,j,i]  = df(p[0,j-1:j+1,i], dy[j-1], dy[j])
          grad[2,0,j,i]  = df(np.array([p[-1,j,i], p[1,j,i]]), dz[-1], dz[0])
          grad[0,-1,j,i] = df(p[-1,j,i-1:i+1], dx[i-1], dx[i])
          grad[1,-1,j,i] = df(p[-1,j-1:j+1,i], dy[j-1], dy[j])
          grad[2,-1,j,i] = df(np.array([p[-2,j,i], p[0,j,i]]), dz[-1], dz[0])
    return grad

  @staticmethod
  @njit(cache=True, nogil=True)
  def gradient_MUSCL(px, py, pz, nx, ny, nz, dx, dy, dz):
    grad = np.zeros((3,nz,ny,nx), dtype=np.float32)
    for k in range(1,nz-1):
      for j in range(1,ny-1):
        for i in range(1,nx-1):
          grad[0,k,j,i] = (-px[k-1,j-1,i-1] + px[k-1,j-1,i]) / (0.5e0 * (dx[i-1] + dx[i]))
          grad[1,k,j,i] = (-py[k-1,j-1,i-1] + py[k-1,j,i-1]) / (0.5e0 * (dy[j-1] + dy[j]))
          grad[2,k,j,i] = (-pz[k-1,j-1,i-1] + pz[k,j-1,i-1]) / (0.5e0 * (dz[k-1] + dz[k]))
    return grad

  def gradient(self, periodic_z=False, TVD=False):
    if TVD is False:
      if self.z is None:
        return self.gradient2D(self.p, self.nx, self.ny, self.dx, self.dy)
      else:
        return self.gradient3D(self.p, self.nx, self.ny, self.nz, self.dx, self.dy, self.dz, periodic_z)
    else:
      px = u_staggered(self.nx, self.ny, self.nz, self.p)
      py = v_staggered(self.nx, self.ny, self.nz, self.p)
      pz = w_staggered(self.nx, self.ny, self.nz, self.p)
      return self.gradient_MUSCL(px, py, pz, self.nx, self.ny, self.nz, self.dx, self.dy, self.dz)


class Vector():
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
  def divergence3D(u, nx, ny, nz, dx, dy, dz, periodic_z=False):
    div = np.zeros((nz,ny,nx), dtype=np.float32)
    for k in range(1,nz-1):
      for j in range(1,ny-1):
        for i in range(1,nx-1):
          div[k,j,i] = df(u[0,k,j,i-1:i+1], dx[i-1], dx[i]) \
                     + df(u[1,k,j-1:j+1,i], dy[j-1], dy[j]) \
                     + df(u[2,k-1:k+1,j,i], dz[k-1], dz[k])
    if periodic_z is True:
      for j in range(1,ny-1):
        for i in range(1,nx-1):
          div[0,j,i]  = df(u[0,0,j,i-1:i+1],  dx[i-1], dx[i]) \
                      + df(u[1,0,j-1:j+1,i],  dy[j-1], dy[j]) \
                      + df(np.array([u[2,-1,j,i], u[2,0,j,i],  u[2,1,j,i]]),  dz[-1], dz[0])
          div[-1,j,i] = df(u[0,-1,j,i-1:i+1], dx[i-1], dx[i]) \
                      + df(u[1,-1,j-1:j+1,i], dy[j-1], dy[j]) \
                      + df(np.array([u[2,-2,j,i], u[2,-1,j,i], u[2,0,j,i]]),  dz[-1], dz[0])
    return div

  @staticmethod
  @njit(cache=True, nogil=True)
  def divergence_MUSCL(nx, ny, nz, dx, dy, dz, us, vs, ws):
    # us[nz-2,ny-2,nx-1], vs[nz-2,ny-1,nx-2], ws[nz-1,ny-2,nx-2]
    div = np.zeros((nz,ny,nx), dtype=np.float32)
    div[1:-1,1:-1,1:-1] = (-us[:,:,:-1] + us[:,:,1:]) / (0.5e0 * (dx[None,None,:-1] + dx[None,None,1:])) \
                        + (-vs[:,:-1,:] + vs[:,1:,:]) / (0.5e0 * (dy[None,:-1,None] + dy[None,1:,None])) \
                        + (-ws[:-1,:,:] + ws[1:,:,:]) / (0.5e0 * (dz[:-1,None,None] + dz[1:,None,None]))
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
  def rotation3D(u, nx, ny, nz, dx, dy, dz, periodic_z=False):
    rot = np.zeros((3,nz,ny,nx), dtype=np.float32)
    for k in range(1,nz-1):
      for j in range(1,ny-1):
        for i in range(1,nx-1):
          rot[0,k,j,i] = df(u[2,k,j-1:j+1,i], dy[j-1], dy[j]) - df(u[1,k-1:k+1,j,i], dz[k-1], dz[k])
          rot[1,k,j,i] = df(u[0,k-1:k+1,j,i], dz[k-1], dz[k]) - df(u[2,k,j,i-1:i+1], dx[i-1], dx[i])
          rot[2,k,j,i] = df(u[1,k,j,i-1:i+1], dx[i-1], dx[i]) - df(u[0,k,j-1:j+1,i], dy[j-1], dy[j])
    if periodic_z is True:
      for j in range(1,ny-1):
        for i in range(1,nx-1):
          U1 = np.array([u[0,-1,j,i], u[0,0,j,i], u[0,1,j,i]])
          U2 = np.array([u[1,-1,j,i], u[1,0,j,i], u[1,1,j,i]])
          rot[0,0,j,i]  = df(u[2,0,j-1:j+1,i], dy[j-1], dy[j]) - df(U2,               dz[-1],  dz[0])
          rot[1,0,j,i]  = df(U1,               dz[-1],  dz[0]) - df(u[2,0,j,i-1:i+1], dx[i-1], dx[i])
          rot[2,0,j,i]  = df(u[1,0,j,i-1:i+1], dx[i-1], dx[i]) - df(u[0,0,j-1:j+1,i], dy[j-1], dy[j])

          U1 = np.array([u[0,-2,j,i], u[0,-1,j,i], u[0,0,j,i]])
          U2 = np.array([u[1,-2,j,i], u[1,-1,j,i], u[1,0,j,i]])
          rot[0,-1,j,i] = df(u[2,-1,j-1:j+1,i], dy[j-1], dy[j]) - df(U2,                dz[-1],  dz[0])
          rot[1,-1,j,i] = df(U1,                dz[-1],  dz[0]) - df(u[2,-1,j,i-1:i+1], dx[i-1], dx[i])
          rot[2,-1,j,i] = df(u[1,-1,j,i-1:i+1], dx[i-1], dx[i]) - df(u[0,-1,j-1:j+1,i], dy[j-1], dy[j])
    return rot

  @staticmethod
  @njit(cache=True, nogil=True)
  def symmetry3D(u, nx, ny, nz, dx, dy, dz, periodic_z=False):
    D2 = np.zeros((nz,ny,nx), dtype=np.float32)
    for k in range(1,nz-1):
      for j in range(1,ny-1):
        for i in range(1,nx-1):
          D11 = df(u[0,k,j,i-1:i+1], dx[i-1], dx[i])
          D22 = df(u[1,k,j-1:j+1,i], dy[j-1], dy[j])
          D33 = df(u[2,k-1:k+1,j,i], dz[k-1], dz[k])
          D12 = 0.5e0 * (df(u[0,k,j-1:j+1,i], dy[j-1], dy[j]) + df(u[1,k,j,i-1:i+1], dx[i-1], dx[i]))
          D23 = 0.5e0 * (df(u[1,k-1:k+1,j,i], dz[k-1], dz[k]) + df(u[2,k,j-1:j+1,i], dy[j-1], dy[j]))
          D31 = 0.5e0 * (df(u[2,k,j,i-1:i+1], dx[i-1], dx[i]) + df(u[0,k-1:k+1,j,i], dz[k-1], dz[k]))
          D2[k,j,i] = 2.e0 * (D11**2 + D22**2 + D33**2) + 4.e0 * (D12**2 + D23**2 + D31**2)
    '''
    if periodic_z is True:
      for j in range(1,ny-1):
        for i in range(1,nx-1):
          D11 = df(u[0,0,j,i-1:i+1], dx[i-1], dx[i])
          D22 = df(u[1,0,j-1:j+1,i], dy[j-1], dy[j])
          U1  = np.array([u[0,-1,j,i], u[0,1,j,i]])
          U2  = np.array([u[1,-1,j,i], u[1,1,j,i]])
          U3  = np.array([u[2,-1,j,i], u[2,1,j,i]])
          D33 = df(U3, dz[-1], dz[0])
          D12 = 0.5e0 * (df(u[0,0,j-1:j+1,i], dy[j-1], dy[j]) + df(u[1,0,j,i-1:i+1], dx[i-1], dx[i]))
          D23 = 0.5e0 * (df(U2, dz[-1], dz[0])                + df(u[2,0,j-1:j+1,i], dy[j-1], dy[j]))
          D31 = 0.5e0 * (df(u[2,0,j,i-1:i+1], dx[i-1], dx[i]) + df(U1, dz[-1], dz[0]))
          D2[0,j,i] = 2.e0 * (D11**2 + D22**2 + D33**2) + 4.e0 * (D12**2 + D23**2 + D31**2)

          D11 = df(u[0,-1,j,i-1:i+1], dx[i-1], dx[i])
          D22 = df(u[1,-1,j-1:j+1,i], dy[j-1], dy[j])
          U1  = np.array([u[0,-2,j,i], u[0,0,j,i]])
          U2  = np.array([u[1,-2,j,i], u[1,0,j,i]])
          U3  = np.array([u[2,-2,j,i], u[2,0,j,i]])
          D33 = df(U3, dz[-1], dz[0])
          D12 = 0.5e0 * (df(u[0,-1,j-1:j+1,i], dy[j-1], dy[j]) + df(u[1,-1,j,i-1:i+1], dx[i-1], dx[i]))
          D23 = 0.5e0 * (df(U2, dz[-1], dz[0])                 + df(u[2,-1,j-1:j+1,i], dy[j-1], dy[j]))
          D31 = 0.5e0 * (df(u[2,-1,j,i-1:i+1], dx[i-1], dx[i]) + df(U1, dz[-1], dz[0]))
          D2[-1,j,i] = 2.e0 * (D11**2 + D22**2 + D33**2) + 4.e0 * (D12**2 + D23**2 + D31**2)
    '''
    return D2

  @staticmethod
  @njit(cache=True, nogil=True)
  def safe_divide(nx, ny, a, b):
    c = np.zeros((ny,nx), dtype=np.float32)
    for j in range(ny):
      for i in range(nx):
        if b[j,i] != 0.e0:
          c[j,i] = a[j,i] / b[j,i]
    return c

  def divergence(self, periodic_z=None, TVD=False):
    if TVD is False:
      if self.z is None:
        return self.divergence2D(self.u, self.nx, self.ny, self.dx, self.dy)
      else:
        return self.divergence3D(self.u, self.nx, self.ny, self.nz, self.dx, self.dy, self.dz, periodic_z)
    else:
      us = u_staggered(self.nx, self.ny, self.nz, self.u[0,:,:,:])
      vs = v_staggered(self.nx, self.ny, self.nz, self.u[1,:,:,:])
      ws = w_staggered(self.nx, self.ny, self.nz, self.u[2,:,:,:])
      return self.divergence_MUSCL(self.nx, self.ny, self.nz, self.dx, self.dy, self.dz, us, vs, ws)

  def rotation(self, periodic_z=None):
    if self.z is None:
      return self.rotation2D(self.u, self.nx, self.ny, self.dx, self.dy)
    else:
      return self.rotation3D(self.u, self.nx, self.ny, self.nz, self.dx, self.dy, self.dz, periodic_z)

  def Fcriterion(self):
    div  = self.divergence3D(self.u, self.nx, self.ny, self.nz, self.dx, self.dy, self.dz)
    W    = self.rotation3D(self.u, self.nx, self.ny, self.nz, self.dx, self.dy, self.dz)
    W2   = np.mean(W[0,1:-1,:,1:-1]**2 + W[1,1:-1,:,1:-1]**2 + W[2,1:-1,:,1:-1]**2, axis=0)
    div2 = np.mean(div[1:-1,:,1:-1]**2, axis=0)
    return self.safe_divide(self.nx-2, self.ny, div2 - W2, div2 + W2)

