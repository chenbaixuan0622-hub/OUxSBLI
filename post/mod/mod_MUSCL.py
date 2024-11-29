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
  d  = -x[:-2] + x[1:]
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


