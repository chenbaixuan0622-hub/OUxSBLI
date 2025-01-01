import numpy as np
from scipy.sparse import diags
from scipy.linalg import lu_factor, lu_solve
import scipy.optimize as op
from scipy.optimize import approx_fprime
import optuna
import os
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
from tqdm import tqdm
from mod.mod_plot import set_Params
from mod.mod_read import getGrid, getVector, extract_number
from mod.mod_POD import make_data, make_grid
from mod.mod_matrix import Poisson_coef, Dirichlet_bc, Neumann_bc, wall_Neumann_bc_A, wall_Neumann_bc_B
from mod.mod_nabla import Vector, Scalar
from mod.mod_Helmholtz import plot_Helmholtz_decomposition


set_Params()

# parameter
Q_dir   = "../3D_solver/TBL/data"
Lx1     = 24.e-3#4.e-3
Lx2     = 40.e-3#20.e-3
Ly1     = 0.e-3
Ly2     = 10.e-3
stridex = 4
endT    = 0.1e-3


def get_and_plot_data(Q_dir, Lx1, Lx2, Ly1, Ly2, stridex, endT):
  Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)
  
  Nx, Ny, Nz, X, Y, Z = getGrid(os.path.join(Q_dir, Q_files[0]))
  nx1 = int(Lx1 / X[-1] * Nx)
  nx2 = int(Lx2 / X[-1] * Nx)
  indicesx = np.arange(nx1, nx2, stridex)
  x = X[indicesx]

  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir, Q_file)
    U, V, _ = getVector(file_path, Nx, Ny, Nz, 'velocity')
    u = U[0,:,indicesx].T
    v = V[0,:,indicesx].T

  y  = np.linspace(Ly1, Ly2, 100)
  Ui = np.zeros((2,len(y),len(x)), dtype=np.float32)
  for i in range(len(x)):
    Ui[0,:,i] = np.interp(y, Y, u[:,i])
    Ui[1,:,i] = np.interp(y, Y, v[:,i])

  crange = np.linspace(0, 500, 100)
  X, Y = np.meshgrid(x*1e3, y*1e3)
  plt.contourf(X, Y, Ui[0,:,:], crange, cmap='jet', extend='both')
  plt.show()
  plt.close()
  return Ui, x, y


def Helmholtz_decomposition(x, y, U):
  div = divergence(U, x, y)
  a, A = Poisson_coef(x, y)
  factor, piv = lu_factor(A)
  f = np.zeros_like(div)
  f[1:-1,1:-1] = div[1:-1,1:-1] / a

  nx, ny = len(x), len(y)
  # boundary conditions
  x0 = np.zeros(ny, dtype=np.float32)
  x1 = np.zeros(ny, dtype=np.float32)
  y0 = np.zeros(nx, dtype=np.float32)
  y1 = np.zeros(nx, dtype=np.float32)
  B = Dirichlet_bc(f, nx, ny, x0, x1, y0, y1)

  phi = lu_solve((factor, piv), B)
  phi = phi.reshape((ny,nx))
  return phi

  
def newton(f, x0, maxiter):
  x     = x0
  xbest = x0
  fbest = f(x0)
  for itr in range(maxiter):
    f_val = f(x)
    grad  = approx_fprime(x, f, 1.e-10)

    def safe_divide(numerator, denominator):
      safe_denom = np.where(denominator < 1e-30, 1.e0, denominator)
      result = numerator / safe_denom
      return np.where(denominator < 1e-30, 0.e0, result)

    tol = 1e-30
    x_new = x - safe_divide(f_val, grad)
    if f_val < tol or np.mean(x - x_new) < tol:
      break
    if f_val < fbest:
      fbest = f_val
      xbest = x
    x = x_new
    print("itr: ", itr, " residual: ", f_val)
  return xbest


def Helmholtz_decomposition_optim(x, y, U):
  nx, ny = len(x), len(y)
  div  = Vector(U, x, y).divergence()
  a, A = Poisson_coef(x, y)
  A    = wall_Neumann_bc_A(A, nx, ny)
  factor, piv = lu_factor(A)
  f = np.zeros_like(div)
  f[1:-1,1:-1] = div[1:-1,1:-1] / a

  # initial values
  x0 = np.zeros(ny, dtype=np.float32)
  x1 = np.zeros(ny, dtype=np.float32)
  #y0 = np.zeros(nx, dtype=np.float32)
  y1 = np.zeros(nx, dtype=np.float32)
  #xbc0 = np.concatenate([x0, x1, y0, y1])
  xbc0 = np.concatenate([x0, x1, y1])

  def residual(xbc):
    #x0, x1, y0, y1 = xbc[:ny], xbc[ny:ny+ny], xbc[2*ny:2*ny+nx], xbc[2*ny+nx:]
    #B = Dirichlet_bc(f, nx, ny, x0, x1, y0, y1)
    x0, x1, y1 = xbc[:ny], xbc[ny:ny+ny], xbc[2*ny:]
    B = wall_Neumann_bc_B(f, nx, ny, x0, x1, y1)
    
    phi = lu_solve((factor, piv), B)
    phi = phi.reshape((nx,ny))
    gradPhi    = Scalar(phi, x, y).gradient()
    rotA       = U - gradPhi
    rotgradPhi = Vector(gradPhi, x, y).rotation()
    divrotA    = Vector(rotA, x, y).divergence()

    res = np.mean(abs(rotgradPhi)) + np.mean(abs(divrotA))
    return res


  print(residual(xbc0))
  xbc = newton(residual, xbc0, maxiter=100)
 
  '''
  # Optuna
  def objective(trial):
    xbc = np.array([trial.suggest_float(f"xbc{i}", -1, 1) for i in range(len(xbc0))])
    return residual(xbc)

  study = optuna.create_study(direction='minimize')
  study.optimize(objective, n_trials=100)
  xbc = np.array([study.best_params[f"x{i}"] for i in range(len(xbc0))])
  '''

  print(residual(xbc))
  #x0, x1, y0, y1 = xbc[:ny], xbc[ny:ny+ny], xbc[2*ny:2*ny+nx], xbc[2*ny+nx:]
  #B = Dirichlet_bc(f, nx, ny, x0, x1, y0, y1)
  x0, x1, y1 = xbc[:ny], xbc[ny:ny+ny], xbc[2*ny:]
  B = wall_Neumann_bc_B(f, nx, ny, x0, x1, y1)
  phi = lu_solve((factor, piv), B)
  phi = phi.reshape((ny,nx))
  return phi


def main(Q_dir, Lx1, Lx2, Ly1, Ly2, stridex, endT):
  U, x, y = get_and_plot_data(Q_dir, Lx1, Lx2, Ly1, Ly2, stridex, endT)
  div = Vector(U, x, y).divergence()
  phi = Helmholtz_decomposition_optim(x, y, U)
  plot_Helmholtz_decomposition(U, phi, x, y)


main(Q_dir, Lx1, Lx2, Ly1, Ly2, stridex, endT)

