import numpy as np
from numba import njit


@njit(cache=True, nogil=True)
def Laplacian(A1, A2, A3, p):
  # p[nz, ny, nx]
  return A1 * (p[1:-1,1:-1,:-2] + p[1:-1,1:-1,2:]) + \
         A2 * (p[1:-1,:-2,1:-1] + p[1:-1,2:,1:-1]) + \
         A3 * (p[:-2,1:-1,1:-1] + p[2:,1:-1,1:-1]) + \
         p[1:-1,1:-1,1:-1]


@njit(cache=True, nogil=True)
def BiCGStab(p, f, dx2dy2, dy2dz2, dz2dx2, err_tol, set_bc):
  A1 = dy2dz2 / (-2.e0 * (dy2dz2 + dz2dx2 + dx2dy2))
  A2 = dz2dx2 / (-2.e0 * (dy2dz2 + dz2dx2 + dx2dy2))
  A3 = dx2dy2 / (-2.e0 * (dy2dz2 + dz2dx2 + dx2dy2))
  B  = f / (-2.e0 * (dy2dz2 + dz2dx2 + dx2dy2))

  r0 = np.zeros_like(p, dtype=np.float32)
  d  = np.zeros_like(p, dtype=np.float32)
  r  = np.zeros_like(p, dtype=np.float32)
  Ad = np.zeros_like(p, dtype=np.float32)
  s  = np.zeros_like(p, dtype=np.float32)
  As = np.zeros_like(p, dtype=np.float32)
  rs = np.zeros_like(p, dtype=np.float32)
  dp = np.zeros_like(p, dtype=np.float32)

  r0[1:-1,1:-1,1:-1] = B[1:-1,1:-1,1:-1] - Laplacian(A1, A2, A3, p)

  r = r0
  d = r0
  err_r = 1.e0

  while err_r > err_tol:
    Ad[1:-1,1:-1,1:-1] = Laplacian(A1, A2, A3, d)
    alpha = np.sum(r0 * r) / np.sum(r0 * Ad)
    s     = r - alpha * Ad
    set_bc(s)
    As[1:-1,1:-1,1:-1] = Laplacian(A1, A2, A3, s)
    w     = np.sum(As * s) / np.sum(As * As)
    dp    = alpha * d + w * s
    p    += dp
    set_bc(p)
    rs    = r
    r     = s - w * As
    beta  = alpha * np.sum(r0 * r) / (w * np.sum(r0 * rs))
    d     = r + beta * (d - w * Ad)
    set_bc(d)
    err_r = np.sqrt(np.sum(dp**2) / np.sum(p**2))
  return p, err_r

