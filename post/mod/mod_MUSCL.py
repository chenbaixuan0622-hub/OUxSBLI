import numpy as np
from numba import njit


@njit(cache=True, fastmath=True, nogil=True)
def minmod2(x, y):
  sgn = np.sign(x)
  return sgn * max(min(abs(x), sgn * y), 0.e0)

  
@njit(cache=True, fastmath=True, nogil=True)
def minmod3(x, y, z):
  sgn = np.sign(x)
  return sgn * max(min(abs(x), sgn * y, sgn * z), 0.e0)


@njit(cache=True, fastmath=True, nogil=True)
def d33(d1, d2, d3):
  da = minmod3(d1, 2.e0 * d2, 2.e0 * d3)
  db = minmod3(d2, 2.e0 * d1, 2.e0 * d3)
  dc = minmod3(d3, 2.e0 * d1, 2.e0 * d2)
  return da - 2.e0 * db + dc


@njit(cache=True, fastmath=True, nogil=True)
def MUSCL3rd(x):
  k = 1.e0 / 3.e0
  b = (3.e0 - k) / (1.e0 - k)
  d = -x[:-2] + x[1:]
  d1 = minmod2(d[0], b * d[1])
  d2 = minmod2(d[1], b * d[0])
  d3 = minmod2(d[2], b * d[1])
  d4 = minmod2(d[1], b * d[2])
  xl = x[1] + 0.25e0 * ((1.e0 - k) * d1 + (1.e0 + k) * d2)
  xr = x[2] - 0.25e0 * ((1.e0 - k) * d3 + (1.e0 + k) * d4)
  return xl, xr


@njit(cache=True, fastmath=True, nogil=True)
def MUSCL4th(x):
  d  = -x[:-1] + x[1:]
  d1 = d[1] - d33(d[0], d[1], d[2]) / 6.e0
  d2 = d[2] - d33(d[1], d[2], d[3]) / 6.e0
  d3 = d[3] - d33(d[2], d[3], d[4]) / 6.e0
  dl = minmod2(d1, 4.e0 * d2)
  dr = minmod2(d2, 4.e0 * d1)
  xl = x[2] + (dl + 2.e0 * dr) / 6.e0
  dl = minmod2(d2, 4.e0 * d3)
  dr = minmod2(d3, 4.e0 * d2)
  xr = x[3] - (dr + 2.e0 * dl) / 6.e0
  return xl, xr

  
@njit(cache=True, fastmath=True, nogil=True)
def u_staggered(nx, ny, nz, u):
  us = np.zeros((nz-2,ny-2,nx-1), dtype=np.float32)
  for k in range(1,nz-1):
    for j in range(1,ny-1):
      for i in range(nx-1):
        if 2 <= i and i <= nx-4:
          ul, ur = MUSCL4th(u[k,j,i-2:i+3])
        elif 1 <= i and i <= nx-3:
          ul, ur = MUSCL3rd(u[k,j,i-1:i+2])
        else:
          um = 0.5e0 * (u[k,j,i] + u[k,j,i+1])
          ul, ur = um, um
        us[k-1,j-1,i] = 0.5e0 * (ul + ur)
  return us


@njit(cache=True, fastmath=True, nogil=True)
def v_staggered(nx, ny, nz, v):
  vs = np.zeros((nz-2,ny-1,nx-2), dtype=np.float32)
  for k in range(1,nz-1):
    for j in range(ny-1):
      for i in range(1,nx-1):
        if 2 <= j and j <= ny-4:
          vl, vr = MUSCL4th(v[k,j-2:j+3,i])
        elif 1 <= j and j <= ny-3:
          vl, vr = MUSCL3rd(v[k,j-1:j+2,i])
        else:
          vm = 0.5e0 * (v[k,j,i] + v[k,j+1,i])
          vl, vr = vm, vm
        vs[k-1,j,i-1] = 0.5e0 * (vl + vr)
  return vs


@njit(cache=True, fastmath=True, nogil=True)
def w_staggered(nx, ny, nz, w):
  ws = np.zeros((nz-1,ny-2,nx-2), dtype=np.float32)
  for k in range(nz-1):
    for j in range(1,ny-1):
      for i in range(1,nx-1):
        if 2 <= k and k <= nz-4:
          wl, wr = MUSCL4th(w[k-2:k+3,j,i])
        elif 1 <= k and k <= nz-3:
          wl, wr = MUSCL3rd(w[k-1:k+2,j,i])
        else:
          wm = 0.5e0 * (w[k,j,i] + w[k+1,j,i])
          wl, wr = wm, wm
        ws[k,j-1,i-1] = 0.5e0 * (wl + wr)
  return ws

