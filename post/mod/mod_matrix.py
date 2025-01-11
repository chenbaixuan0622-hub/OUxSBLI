import numpy as np
from scipy.sparse import diags
from scipy.linalg import lu_factor, lu_solve
import matplotlib.pyplot as plt


def create_matrix(nx, ny, A1, A2, A3, A4):
  # A1, A2, A3, A4 : scalar or 1d array
  n    = nx * ny
  diag = np.ones(n, dtype=np.float32)

  x1 = np.zeros(n-1,  dtype=np.float32)
  x2 = np.zeros(n-1,  dtype=np.float32)
  y1 = np.zeros(n-nx, dtype=np.float32)
  y2 = np.zeros(n-nx, dtype=np.float32)

  idx1 = np.arange(1, n,    dtype=np.int32)
  idx2 = np.arange(0, n-1,  dtype=np.int32)
  idy1 = np.arange(0, n-nx, dtype=np.int32)
  idy2 = np.arange(0, n-nx, dtype=np.int32)
  # boundary condition for x = 0
  idx1[np.arange(1,n)    % nx == 0] = 0
  idx2[np.arange(0,n-1)  % nx == 0] = 0
  idy1[np.arange(0,n-nx) % nx == 0] = 0
  idy2[np.arange(0,n-nx) % nx == 0] = 0
  # boundary condition for x = -1
  idx1[np.arange(1,n)    % nx == nx-1] = 0
  idx2[np.arange(0,n-1)  % nx == nx-1] = 0
  idy1[np.arange(0,n-nx) % nx == nx-1] = 0
  idy2[np.arange(0,n-nx) % nx == nx-1] = 0
  # boundary condition for y = 0
  idx1[:nx-1] = 0
  idx2[:nx]   = 0
  idy2[:nx]   = 0
  # boundary condition for y = -1
  idx1[n-nx-1:]  = 0
  idx2[n-nx:]    = 0
  idy1[n-nx-nx:] = 0

  idx1 = np.delete(idx1, np.where(idx1 == 0))
  idx2 = np.delete(idx2, np.where(idx2 == 0))
  idy1 = np.delete(idy1, np.where(idy1 == 0))
  idy2 = np.delete(idy2, np.where(idy2 == 0))

  x1[idx1-1] = A1
  x2[idx2]   = A2

  y1[idy1]   = A3
  y2[idy2]   = A4

  diagonals = [diag, x1, x2, y1, y2]
  offsets   = [0, -1, 1, -nx, nx]

  A = diags(diagonals, offsets, format="csr")
  return A


def Dirichlet_bc(B, nx, ny, x0, x1, y0, y1):
  # B : 2d array
  B = B.flatten()
  B[::nx]     = x0
  B[nx-1::nx] = x1
  B[:nx]      = y0
  B[-nx:]     = y1
  return B


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


def wall_Neumann_bc_A(A, nx, ny):
  for i in range(1,nx-1):
    A[i,i+nx] = -1.e0
  return A


def wall_Neumann_bc_B(B, nx, ny, x0, x1, y1):
  # B : 2d array
  B = B.flatten()
  B[::nx]     = x0
  B[nx-1::nx] = x1
  B[:nx]      = 0.e0
  B[-nx:]     = y1
  return B


def Poisson_coef(x, y):
  dx = -x[:-1] + x[1:]
  dy = -y[:-1] + y[1:]

  Ax = 0.5e0 * (dx[:-1] + dx[1:]) * dx[:-1] * dx[1:]
  Ay = 0.5e0 * (dy[:-1] + dy[1:]) * dy[:-1] * dy[1:]
  A  = - (dx[None,:-1] + dx[None,1]) / Ax[None,:] \
       - (dy[:-1,None] + dy[1,None]) / Ay[:,None]
  A1 = dx[None,1:]  / (A * Ax[None,:]) 
  A2 = dx[None,:-1] / (A * Ax[None,:])
  A3 = dy[1:,None]  / (A * Ay[:,None])
  A4 = dy[:-1,None] / (A * Ay[:,None])
  A1 = A1.flatten()
  A2 = A2.flatten()
  A3 = A3.flatten()
  A4 = A4.flatten()
  return A, create_matrix(len(x), len(y), A1, A2, A3, A4).toarray()


def Poisson_coef_FFT(x, y, kz):
  dx = -x[:-1] + x[1:]
  dy = -y[:-1] + y[1:]

  Ax = 0.5e0 * (dx[:-1] + dx[1:]) * dx[:-1] * dx[1:]
  Ay = 0.5e0 * (dy[:-1] + dy[1:]) * dy[:-1] * dy[1:]
  A  = - (dx[None,:-1] + dx[None,1]) / Ax[None,:] \
       - (dy[:-1,None] + dy[1,None]) / Ay[:,None] - kz**2
  A1 = dx[None,1:]  / (A * Ax[None,:]) 
  A2 = dx[None,:-1] / (A * Ax[None,:])
  A3 = dy[1:,None]  / (A * Ay[:,None])
  A4 = dy[:-1,None] / (A * Ay[:,None])
  A1 = A1.flatten()
  A2 = A2.flatten()
  A3 = A3.flatten()
  A4 = A4.flatten()
  return A, create_matrix(len(x), len(y), A1, A2, A3, A4).toarray()


def Poisson_cyclic(x, y, z, f):
  nx    = len(x)
  ny    = len(y)
  nz    = len(z)
  f_hat = np.fft.fft(f, axis=0)
  p_hat = np.zeros_like(f, dtype=np.complex64)
  for k in range(nz):
    kz = 2.e0 * np.pi * k / nz if k <= nz // 2 else 2.e0 * np.pi * (k - nz) / nz
    a, A = Poisson_coef_FFT(x, y, kz)
    
    f_hat[k,1:-1,1:-1] /= a
    B = f_hat[k,:,:].flatten()

    factor, piv  = lu_factor(A)
    p_hat[k,:,:] = lu_solve((factor, piv), B).reshape(ny,nx)
  
  phi = np.real(np.fft.ifft(p_hat, axis=0))
  return phi


def check_matrix(nx, ny, A, B):
  AB = np.append(A, B.reshape(nx*ny,-1), axis=-1)

  plt.imshow(AB, cmap='jet', interpolation='none')
  plt.show()
  plt.close()

