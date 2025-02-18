import numpy as np
from numba import njit
import matplotlib.pyplot as plt


def RK4(func, X0, sets):
  """
  Runge Kutta 4 solver.
  """
  n_x, dt, n_data, h, c, A, sigma = sets
  func_sets = (n_x, h, c, A, sigma)
  X  = np.zeros([n_data, len(X0)])
  X[0] = X0
  ti = 0
  for i in range(n_data-1):
    k1 = func(X[i], ti, func_sets)
    k2 = func(X[i] + dt/2. * k1, ti + dt/2., func_sets)
    k3 = func(X[i] + dt/2. * k2, ti + dt/2., func_sets)
    k4 = func(X[i] + dt    * k3, ti + dt, func_sets)
    X[i+1] = X[i] + dt / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
    ti += dt
  return X


@njit(cache=True, nogil=True)
def lorenz(nt=1000, dt=0.01e0, x0=0.e0, y0=1.e0, z0=1.05e0, s=10.e0, r=28.e0, b=8.e0/3.e0):
  x = np.zeros(nt+1)
  y = np.zeros(nt+1)
  z = np.zeros(nt+1)
  x[0] = x0
  y[0] = y0
  z[0] = z0
  for i in range(nt):
    dot_x = s * (y[i] - x[i])
    dot_y = r * x[i] - y[i] - x[i] * z[i]
    dot_z = x[i] * y[i] - b * z[i]
    x[i+1] = x[i] + dot_x * dt
    y[i+1] = y[i] + dot_y * dt
    z[i+1] = z[i] + dot_z * dt
  return x, y, z


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

