import numpy as np
from scipy.sparse import diags
from scipy.linalg import lu_factor, lu_solve
from numba import njit
import matplotlib.pyplot as plt


@njit(cache=True, nogil=True)
def create_matrix(nx, ny, a1, a2, a3, a4):
  # A1, A2, A3, A4 : scalar or 1d array [(nx-2) * (ny-2)]
  n    = nx * ny
  diag = np.ones(n, dtype=np.float32)

  A1 = np.zeros((ny,nx), dtype=np.float32)
  A2 = np.zeros((ny,nx), dtype=np.float32)
  A3 = np.zeros((ny,nx), dtype=np.float32)
  A4 = np.zeros((ny,nx), dtype=np.float32)
  
  A1[1:-1,1:-1] = a1
  A2[1:-1,1:-1] = a2
  A3[1:-1,1:-1] = a3
  A4[1:-1,1:-1] = a4

  A1 = A1.flatten()
  A2 = A2.flatten()
  A3 = A3.flatten()
  A4 = A4.flatten()

  x1 = A1[1:]
  x2 = A2[:-1]

  y1 = A3[nx:]
  y2 = A4[:-nx]

  A = np.diag(diag) + np.diag(x1, -1) + np.diag(x2, 1) + np.diag(y1, -nx) + np.diag(y2, nx)
  return A


@njit(cache=True, nogil=True)
def Dirichlet_bc(B, nx, ny, x0, x1, y0, y1):
  # B : 2d array
  B = B.flatten()
  B[::nx]     = x0
  B[nx-1::nx] = x1
  B[:nx]      = y0
  B[-nx:]     = y1
  return B


@njit(cache=True, nogil=True)
def Neumann_bc(A, B, nx, ny):
  # B : 2d array
  B = B.flatten()
  B[::nx]     = 0.e0
  B[nx-1::nx] = 0.e0
  B[:nx]      = 0.e0
  B[-nx:]     = 0.e0
  for i in range(nx*ny-nx+1,nx*ny-1):
    A[i,i-nx] = -1.e0
  for i in range(nx*ny-1):
    if i % nx == 0:
      A[i,i+1] = -1.e0
  for i in range(1,nx*ny):
    if i % nx == nx-1:
      A[i,i-1] = -1.e0
  return A, B


@njit(cache=True, nogil=True)
def wall_Neumann_bc_A(A, nx, ny):
  for i in range(1,nx-1):
    A[i,i+nx] = -1.e0
  return A


@njit(cache=True, nogil=True)
def wall_Neumann_bc_B(B, nx, ny, x0, x1, y1):
  # B : 2d array
  B = B.flatten()
  B[::nx]     = x0
  B[nx-1::nx] = x1
  B[:nx]      = 0.e0
  B[-nx:]     = y1
  return B


@njit(cache=True, nogil=True)
def Poisson_coef(x, y, mz=0.e0):
  dx = -x[:-1] + x[1:]
  dy = -y[:-1] + y[1:]

  Ax = 0.5e0 * (dx[:-1] + dx[1:]) * dx[:-1] * dx[1:]
  Ay = 0.5e0 * (dy[:-1] + dy[1:]) * dy[:-1] * dy[1:]
  A  = - (dx[None,:-1] + dx[None,1]) / Ax[None,:] \
       - (dy[:-1,None] + dy[1,None]) / Ay[:,None] + mz
  A1 = dx[None,1:]  / (A * Ax[None,:]) 
  A2 = dx[None,:-1] / (A * Ax[None,:])
  A3 = dy[1:,None]  / (A * Ay[:,None])
  A4 = dy[:-1,None] / (A * Ay[:,None])
  return A, create_matrix(len(x), len(y), A1, A2, A3, A4)


def check_matrix(nx, ny, A, B=None):
  if B is None:
    plt.imshow(A, cmap='jet', interpolation='none')
    plt.show()
    plt.close()
  else:
    AB = np.append(A, B.reshape(nx*ny,-1), axis=-1)

    plt.imshow(AB, cmap='jet', interpolation='none')
    plt.show()
    plt.close()

