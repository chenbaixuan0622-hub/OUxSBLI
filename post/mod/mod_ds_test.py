import numpy as np
from numba import njit
import matplotlib.pyplot as plt


class RK:
  @staticmethod
  @njit(cache=True, nogil=True)
  def _RK44(RHS, params, nt, dt, x0):
    x = np.zeros((nt+1,len(x0)))
    x[0,:] = x0[:]
    Dx = np.zeros(3)
    for i in range(nt):
      dx = RHS(params, x[i,:])
      x[i+1,:] = x[i,:] + 0.5e0 * dt * dx[:]
      Dx[:] = dx[:]
      dx = RHS(params, x[i+1,:])
      x[i+1,:] = x[i,:] + 0.5e0 * dt * dx[:]
      Dx[:] += 2.e0 * dx[:]
      dx = RHS(params, x[i+1,:])
      x[i+1,:] = x[i,:] + dt * dx[:]
      Dx[:] += 2.e0 * dx[:]
      dx = RHS(params, x[i+1,:])
      Dx[:] += dx[:]
      x[i+1,:] = x[i,:] + dt * Dx[:] / 6.e0
    return x

  @staticmethod
  @njit(cache=True, nogil=True)
  def _error(RHS, params, dim, k1, k2, a11, a12, a21, a22, x, dt):
    err = np.empty(2 * dim)
    f1 = k1 - RHS(params, x + dt * (a11 * k1 + a12 * k2))
    f2 = k2 - RHS(params, x + dt * (a21 * k1 + a22 * k2))
    for i in range(dim):
      err[i]     = f1[i]
      err[dim+i] = f2[i]
    return err

  @staticmethod
  @njit(cache=True, nogil=True)
  def _GaussStep(RHS, Jacobian, error, params, dt, x, damping=1.e0, max_itr=100, tol=1e-12):
    sqrt3 = np.sqrt(3.e0)
    dim   = x.shape[0]
    k  = RHS(params, x)
    c1 = 0.5e0 - sqrt3 / 6.e0
    c2 = 0.5e0 + sqrt3 / 6.e0
    x1 = x + dt * c1 * k
    x2 = x + dt * c2 * k
    k1 = RHS(params, x1)
    k2 = RHS(params, x2)
    a11 = 0.25e0
    a12 = 0.25e0 - sqrt3 / 6.e0
    a21 = 0.25e0 + sqrt3 / 6.e0
    a22 = 0.25e0
    err = error(RHS, params, dim, k1, k2, a11, a12, a21, a22, x, dt)
    itr = 0
    while np.linalg.norm(err) > tol and itr < max_itr:
      itr += 1
      J1 = Jacobian(params, x + dt * (a11 * k1 + a12 * k2))
      J2 = Jacobian(params, x + dt * (a21 * k1 + a22 * k2))
      J = np.zeros((2*dim,2*dim))
      I = np.eye(dim)
      for i in range(dim):
        for j in range(dim):
          J[i,j]         = I[i,j] - dt * a11 * J1[i,j]
          J[i,j+dim]     = -dt * a12 * J1[i,j]
          J[i+dim,j]     = -dt * a21 * J2[i,j]
          J[i+dim,j+dim] = I[i,j] - dt * a22 * J2[i,j]
      delta = np.linalg.solve(J, err)
      k1 -= damping * delta[:dim]
      k2 -= damping * delta[dim:]
      err = error(RHS, params, dim, k1, k2, a11, a12, a21, a22, x, dt)
    return x + 0.5e0 * dt * (k1 + k2)

  @staticmethod
  @njit(cache=True, nogil=True)
  def _GaussRK2(GaussStep, RHS, Jacobian, error, params, nt, dt, x0):
    x = np.zeros((nt+1,len(x0)))
    x[0,:] = x0[:]
    for i in range(nt):
      x[i+1,:] = GaussStep(RHS, Jacobian, error, params, dt, x[i,:])
    return x



class lorenz(RK):
  def __init__(self, s=10.e0, r=28.e0, b=8.e0/3.e0):
    self.s = s
    self.r = r
    self.b = b

  @staticmethod
  @njit(cache=True, nogil=True)
  def _RHS(params, x):
    s, r, b = params
    dx = s * (x[1] - x[0])
    dy = r * x[0] - x[1] - x[0] * x[2]
    dz = x[0] * x[1] - b * x[2]
    return np.array([dx, dy, dz])

  @staticmethod
  @njit(cache=True, nogil=True)
  def _Jacobian(params, x):
    s, r, b = params
    J = np.array([[   -s,           s,  0.e0], 
                  [r - x[2],    -1.e0, -x[0]], 
                  [    x[1],     x[0],    -b]])
    return J

  def RK44(self, nt=10000, dt=0.001e0, x0=0.e0, y0=1.e0, z0=1.05e0):
    params  = (self.s, self.r, self.b)
    x = self._RK44(self._RHS, params, nt, dt, np.array([x0, y0, z0]))
    return x[:,0], x[:,1], x[:,2]

  def GaussRK2(self, nt=10000, dt=0.001e0, x0=0.e0, y0=1.e0, z0=1.05e0):
    params  = (self.s, self.r, self.b)
    x = self._GaussRK2(self._GaussStep, self._RHS, self._Jacobian, self._error, 
                       params, nt, dt, np.array([x0, y0, z0]))
    return x[:,0], x[:,1], x[:,2]


@njit(cache=True, nogil=True)
def logistic_map(nt, x, r):
  data = np.zeros(n, dtype=np.float64)
  for i in range(n):
    x = r * x * (1.e0 - x)
    data[i] = x
  return data


@njit(cache=True, nogil=True)
def coupling_system(nt, x0, y0, bxy, byx, gx=3.7e0, gy=3.72e0):
  x = np.zeros(nt+1)
  y = np.zeros(nt+1)
  x[0] = x0
  y[0] = y0
  for i in range(nt):
    x[i+1] = x[i] * (gx - (gx - byx) * x[i] - byx * y[i]) + np.random.normal(0.e0, 0.01e0)
    y[i+1] = y[i] * (gy - (gy - bxy) * y[i] - bxy * x[i]) + np.random.normal(0.e0, 0.01e0)
  return x, y


@njit(cache=True, nogil=True)
def discrete_logistic_map(nt, x0, y0, z0, bxy, byz, bxz, gx=3.7e0, gy=3.72e0, gz=3.78e0):
  x = np.zeros(nt+1)
  y = np.zeros(nt+1)
  z = np.zeros(nt+1)
  x[0] = x0
  y[0] = y0
  z[0] = z0
  for i in range(nt):
    x[i+1] = gx * x[i] * (1.e0 - x[i]) + np.random.normal(0.e0, 0.005e0)
    y[i+1] = gy * y[i] * (1.e0 - (1.e0 - bxy / gy) * y[i] - bxy * x[i] / gy) + np.random.normal(0.e0, 0.005e0)
    z[i+1] = gz * z[i] * (1.e0 - (1.e0 - (bxz + byz) / gz) * z[i] - bxz * x[i] / gz - byz * y[i] / gz) + np.random.normal(0.e0, 0.005e0)
  return x, y, z


def coupled_Lorenz(x, t, sets):
  # Lorenz system is time-invariant, no need to use the input t
  n_x, h, c, A, sigma = sets
  r = np.zeros(x.shape)
  noise = np.random.normal(loc=0, scale=sigma, size=3*n_x)
  for k in range(n_x):
    r[k*3+0]= -10*(x[k*3+0]-x[k*3+1] +c*A[k]@(x[1::3]-x[k*3+1])) + noise[0+k*3]
    r[k*3+1]= 28*(1+h[k])*x[k*3+0]-x[k*3+1] - x[k*3+0]*x[k*3+2] + noise[1+k*3]
    r[k*3+2]= -8/3*x[k*3+2]+x[k*3+0]*x[k*3+1] + noise[2+k*3]
  return r


@njit(cache=True, nogil=True)
def coupled_lorenz9(nt=1000, dt=0.01e0, a=10.e0, b=28.e0, c=8.e0/3.e0, bxy=0.3e0, bxz=0.2e0):
  byz = bxy
  x1 = np.zeros(nt+1)
  x2 = np.zeros(nt+1)
  x3 = np.zeros(nt+1)
  y1 = np.zeros(nt+1)
  y2 = np.zeros(nt+1)
  y3 = np.zeros(nt+1)
  z1 = np.zeros(nt+1)
  z2 = np.zeros(nt+1)
  z3 = np.zeros(nt+1)
  x1[0] = 0.e0
  x2[0] = 1.e0
  x3[0] = 1.05e0
  y1[0] = 0.e0
  y2[0] = 1.e0
  y3[0] = 1.05e0
  z1[0] = 0.e0
  z2[0] = 1.e0
  z3[0] = 1.05e0
  for i in range(nt):
    dot_x1 = a * (x2[i] - x1[i]) + np.random.normal(0.e0, 1.e0) 
    dot_x2 = x1[i] * (b - x3[i]) - x2[i]
    dot_x3 = x1[i] * x2[i] - c * x3[i]
    dot_y1 = a * (y2[i] - y1[i]) + bxy * x1[i] + np.random.normal(0.e0, 1.e0)
    dot_y2 = y1[i] * (b - y3[i]) - y2[i]
    dot_y3 = y1[i] * y2[i] - c * y3[i]
    dot_z1 = a * (z2[i] - z1[i]) + bxz * x1[i] + byz * y1[i] + np.random.normal(0.e0, 1.e0)
    dot_z2 = z1[i] * (b - z3[i]) - z2[i]
    dot_z3 = z1[i] * z2[i] - c * z3[i]
    x1[i+1] = x1[i] + dot_x1 * dt
    x2[i+1] = x2[i] + dot_x2 * dt
    x3[i+1] = x3[i] + dot_x3 * dt
    y1[i+1] = y1[i] + dot_y1 * dt
    y2[i+1] = y2[i] + dot_y2 * dt
    y3[i+1] = y3[i] + dot_y3 * dt
    z1[i+1] = z1[i] + dot_z1 * dt
    z2[i+1] = z2[i] + dot_z2 * dt
    z3[i+1] = z3[i] + dot_z3 * dt
  return x1, x2, x3, y1, y2, y3, z1, z2, z3


def gen_data_norm(n_x, A, n_data, c = 0.3, sigma_dyn = 0, sigma_obs = 0, dt = 0.01):
  h  = 2*(np.random.rand(n_x)-0.5)*0.06
    
  setting = (n_x, dt, n_data, h, c, A, sigma_dyn)
  u0 = np.array([7.432487609628195, 10.02071718705213, 29.62297428638419])
  x0 = u0 + np.random.random((n_x,3))*0.3

  # X = gen_data(setting, x0)
  X = RK4(coupled_Lorenz, x0.flatten(), setting)
    
  q1 = np.percentile(X, 25, interpolation='midpoint', axis=0)
  q2 = np.percentile(X, 50, interpolation='midpoint', axis=0)
  q3 = np.percentile(X, 75, interpolation='midpoint', axis=0)

  X = (X - q2)/(q3 - q1)
    
  if sigma_obs != 0:
    noise_obs = np.random.normal(loc=0, scale=sigma_obs, size=X.shape)
    X += noise_obs
  return X

