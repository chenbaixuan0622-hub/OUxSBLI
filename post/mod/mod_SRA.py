import numpy as np
from numba import njit


@njit(cache=True, fastmath=True, nogil=True)
def corr_x(Nx, ny1, ny2, Nz, uf, vf, wf, uu, vv, ww):
    # x corr
    for k in range(Nz):
      for i in range(Nx):
        uu[k,0,i] += uf[k,ny1,0] * uf[k,ny1,i]
        uu[k,1,i] += uf[k,ny2,0] * uf[k,ny2,i]
        vv[k,0,i] += vf[k,ny1,0] * vf[k,ny1,i]
        vv[k,1,i] += vf[k,ny2,0] * vf[k,ny2,i]
        ww[k,0,i] += wf[k,ny1,0] * wf[k,ny1,i]
        ww[k,1,i] += wf[k,ny2,0] * wf[k,ny2,i]


@njit(cache=True, fastmath=True, nogil=True)
def corr_z(Nx, ny1, ny2, Nz, uf, vf, wf, uu, vv, ww):
    # z corr
    for k in range(Nz):
      for i in range(Nx):
        uu[k,0,i] += uf[0,ny1,i] * uf[k,ny1,i]
        uu[k,1,i] += uf[0,ny2,i] * uf[k,ny2,i]
        vv[k,0,i] += vf[0,ny1,i] * vf[k,ny1,i]
        vv[k,1,i] += vf[0,ny2,i] * vf[k,ny2,i]
        ww[k,0,i] += wf[0,ny1,i] * wf[k,ny1,i]
        ww[k,1,i] += wf[0,ny2,i] * wf[k,ny2,i]


@njit(cache=True, fastmath=True, nogil=True)
def SRA(Nx, Ny, Nz, uF, vF, TF, TtF, uu, vv, uv, TT, TtTt, uT, vT, vTt):
  for k in range(Nz):
    for j in range(Ny-1):
      for i in range(Nx):
        uu[k,j,i]   +=  uF[k,j+1,i]**2
        vv[k,j,i]   +=  vF[k,j+1,i]**2
        uv[k,j,i]   +=  uF[k,j+1,i] *  vF[k,j+1,i]
        TT[k,j,i]   +=  TF[k,j+1,i]**2
        TtTt[k,j,i] += TtF[k,j+1,i]**2
        uT[k,j,i]   +=  uF[k,j+1,i] *  TF[k,j+1,i]
        vT[k,j,i]   +=  vF[k,j+1,i] *  TF[k,j+1,i]
        vTt[k,j,i]  +=  vF[k,j+1,i] * TtF[k,j+1,i]

