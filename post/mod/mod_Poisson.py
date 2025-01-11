import numpy as np
from scipy.fftpack import dct, idct, dst, idst
from scipy.linalg import solve_banded
import matplotlib.pyplot as plt


class Poisson:
  def __init__(self, f, x, y=None, z=None):
    self.f  = f
    self.x  = x
    self.y  = y
    self.z  = z
    self.dx = -x[0] + x[1]
    self.nx = len(x)
    if self.nx != f.shape[-1]:
      raise Exception('Invalid length')
    if self.y is not None:
      self.dy = -y[0] + y[1]
      self.ny = len(y)
      if self.ny != f.shape[-2]:
        raise Exception('Invalid length')
    if self.z is not None:
      self.dz = -z[0] + z[1]
      self.nz = len(z)
      if self.nz != f.shape[-3]:
        raise Exception('Invalid length')


  def Dirichlet(self, p1, p2, f_hat=None, my=0.e0):
    if f_hat is None:
      f = self.f
    else:
      f = f_hat
    upper =         np.ones(self.nx) / self.dx**2
    diag  = -2.e0 * np.ones(self.nx) / self.dx**2 + my
    lower =         np.ones(self.nx) / self.dx**2
    upper[0]  = 0.e0
    lower[-1] = 0.e0
    diag[0]   = -3.e0 / self.dx**2
    diag[-1]  = -3.e0 / self.dx**2
    f[0]      = f[0]  - 2.e0 * p1 / self.dx**2
    f[-1]     = f[-1] - 2.e0 * p2 / self.dx**2
    ab = np.stack([upper, diag, lower])
    return solve_banded((1,1), ab, f)


  def Neumann(self, dp1, dp2, f_hat=None, my=0.e0):
    if f_hat is None:
      f = self.f
    else:
      f = f_hat
    f[0]  += dp1 / self.dx
    f[-1] -= dp2 / self.dx
    fhat = dct(f, norm='ortho')
    k  = np.linspace(0.e0, self.nx-1, self.nx)
    mx = 2.e0 * (np.cos(np.pi * k / self.nx) - 1.e0) / self.dx**2
    p_hat     = np.zeros_like(fhat)
    p_hat[1:] = fhat[1:] / (mx[1:] + my)
    p_hat[0]  = 0.e0
    return idct(p_hat, norm='ortho')


  def Dirichlet_Neumann(self, dp1, p2, f_hat=None, my=0.e0):
    if f_hat is None:
      f = self.f
    else:
      f = f_hat
    f[0]  += dp1 / self.dx
    f[-1] -= 2.e0 * p2 / self.dx**2
    fhat = dct(f, type=4, norm='ortho')
    k  = np.linspace(0.e0, self.nx-1, self.nx) + 0.5e0
    mx = 2.e0 * (np.cos(np.pi * k / self.nx) - 1.e0) / self.dx**2
    p_hat = fhat / (mx + my)
    return idct(p_hat, type=4, norm='ortho')


  # y Dirichlet
  def __y_Dirichlet(self, f, v1, v2, py1, py2, Poisson_1D, mz=0.e0):
    f[0,:]  -= 2.e0 * py1 / self.dy**2
    f[-1,:] -= 2.e0 * py2 / self.dy**2
    f_hat   = dst(f,  type=2, norm='ortho', axis=0)
    v1_hat  = dst(v1, type=2, norm='ortho')
    v2_hat  = dst(v2, type=2, norm='ortho')
    p_hat   = np.zeros_like(f_hat)
    for k in range(self.ny):
      my = 2.e0 * (np.cos(np.pi * (k+1) / self.ny) - 1.e0) / self.dy**2 + mz
      p_hat[k,:] = Poisson_1D(v1[k], v2[k], f_hat[k,:], my)
    return dst(p_hat, type=3, axis=0, norm='ortho')


  # y Neumann
  def __y_Neumann(self, f, v1, v2, dpy1, dpy2, Poisson_1D, mz=0.e0):
    f[0,:]  += dpy1 / self.dy
    f[-1,:] -= dpy2 / self.dy
    f_hat   = dct(f, axis=0, norm='ortho')
    v1_hat  = dct(v1, norm='ortho')
    v2_hat  = dct(v2, norm='ortho')
    p_hat   = np.zeros_like(f_hat)
    for k in range(self.ny):
      my = 2.e0 * (np.cos(np.pi * k / self.ny) - 1.e0) / self.dy**2 + mz
      p_hat[k,:] = Poisson_1D(v1[k], v2[k], f_hat[k,:], my)
    return idct(p_hat, axis=0, norm='ortho')
  

  @staticmethod
  def __z_periodic(f, nz, dz, v1, v2, v3, v4, Poisson_2D):
    f_hat  = np.fft.fft(f,  axis=0)
    v1_hat = np.fft.fft(v1, axis=0)
    v2_hat = np.fft.fft(v2, axis=0)
    v3_hat = np.fft.fft(v3, axis=0)
    v4_hat = np.fft.fft(v4, axis=0)
    p_hat  = np.zeros_like(f_hat)
    for k in range(nz):
      mw = -2.e0 * (1.e0 - np.cos(2.e0 * np.pi * k / nz)) / dz**2
      p_hat[k,:,:] = Poisson_2D(v1_hat[k,:], v2_hat[k,:], v3_hat[k,:], v4_hat[k,:], f_hat[k,:,:], mw)
    return np.real(np.fft.ifft(p_hat, axis=0))


  # 2D Poisson equation
  def x_Neumann_y_Dirichlet(self, v1, v2, v3, v4, f_hat=None, mw=None):
    if f_hat is None and mw is None:
      if len(v1) == self.ny and len(v2) == self.ny and len(v3) == self.nx and len(v4) == self.nx:
        return self.__y_Dirichlet(self.f, v1, v2, v3, v4, self.Neumann)
      else:
        raise Exception('Invalid length')
    elif f_hat is not None and mw is not None:
      return self.__y_Dirichlet(f_hat, v1, v2, v3, v4, self.Neumann, mw)
    else:
      raise Exception('Invalid arguments')


  def x_Dirichlet_y_Dirichlet(self, v1, v2, v3, v4, f_hat=None, mw=None):
    if f_hat is None and mw is None:
      if len(v1) == self.ny and len(v2) == self.ny and len(v3) == self.nx and len(v4) == self.nx:
        return self.__y_Dirichlet(self.f, v1, v2, v3, v4, self.Dirichlet)
      else:
        raise Exception('Invalid length')
    elif f_hat is not None and mw is not None:
      return self.__y_Dirichlet(f_hat, v1, v2, v3, v4, self.Dirichlet, mw)
    else:
      raise Exception('Invalid arguments')


  def x_Neumann_y_Neumann(self, v1, v2, v3, v4, f_hat=None, mw=None):
    if f_hat is None and mw is None:
      if len(v1) == self.ny and len(v2) == self.ny and len(v3) == self.nx and len(v4) == self.nx:
        return self.__y_Neumann(self.f, v1, v2, v3, v4, self.Neumann)
      else:
        raise Exception('Invalid length')
    elif f_hat is not None and mw is not None:
      return self.__y_Neumann(f_hat, v1, v2, v3, v4, self.Neumann, mw)
    else:
      raise Exception('Invalid arguments')


  def x_Dirichlet_y_Neumann(self, v1, v2, v3, v4, f_hat=None, mw=None):
    if f_hat is None and mw is None:
      if len(v1) == self.ny and len(v2) == self.ny and len(v3) == self.nx and len(v4) == self.nx:
        return self.__y_Neumann(self.f, v1, v2, v3, v4, self.Dirichlet)
      else:
        raise Exception('Invalid length')
    elif f_hat is not None and mw is not None:
      return self.__y_Neumann(f_hat, v1, v2, v3, v4, self.Dirichlet, mw)
    else:
      raise Exception('Invalid arguments')


  # 3D Poisson equation z periodic
  def x_Neumann_y_Dirichlet_z_periodic(self, dp1x, dp2x, p1y, p2y):
    return self.__z_periodic(self.f, self.nz, self.dz, dp1x, dp2x, p1y, p2y, self.x_Neumann_y_Dirichlet)


  def x_Dirichlet_y_Dirichlet_z_periodic(self, dp1x, dp2x, p1y, p2y):
    return self.__z_periodic(self.f, self.nz, self.dz, dp1x, dp2x, p1y, p2y, self.x_Dirichlet_y_Dirichlet)

  def x_Neumann_y_Neumann_z_periodic(self, dp1x, dp2x, dp1y, dp2y):
    if np.shape(dp1x) == np.shape(self.f[:,:,0]) and np.shape(dp2x) == np.shape(self.f[:,:,0]) \
      and np.shape(dp1y) == np.shape(self.f[:,0,:]) and np.shape(dp2y) == np.shape(self.f[:,0,:]):
      return self.__z_periodic(self.f, self.nz, self.dz, dp1x, dp2x, dp1y, dp2y, self.x_Neumann_y_Neumann)
    else:
      raise Exception('Invalid shape')


  def x_Dirichlet_y_Neumann_z_periodic(self, p1x, p2x, dp1y, dp2y):
    if np.shape(p1x) == np.shape(self.f[:,:,0]) and np.shape(p2x) == np.shape(self.f[:,:,0]) \
      and np.shape(dp1y) == np.shape(self.f[:,0,:]) and np.shape(dp2y) == np.shape(self.f[:,0,:]):
      return self.__z_periodic(self.f, self.nz, self.dz, p1x, p2x, dp1y, dp2y, self.x_Dirichlet_y_Neumann)
    else:
      raise Exception('Invalid shape')

