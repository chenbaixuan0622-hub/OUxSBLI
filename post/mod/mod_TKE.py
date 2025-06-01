import numpy as np
import vtk
import os
from numba import jit, njit, prange
from mod.mod_turb_stat import non_dim_tbl, Sutherland, tau


@njit(cache=True, fastmath=True, nogil=True)
def BudgetTerms(x, y, z, rho, u, v, w, p, Rho, U, V, W, P, UF, VF, WF, Pro, T, PI, D, eps):
  nx = len(x)
  ny = len(y)
  nz = len(z)
  dx = -x[0] + x[1]
  dz = -z[0] + z[1]
  Rgas = 287.03
  # TKE Budget
  Pro += Production(nx, ny, nz, dx, y, dz, rho, u, v, w, UF, VF, WF)
  T   += Transport(nx, ny, nz, dx, y, dz, rho, u, v, w, UF, VF, WF)
  PI  += Pressure(nx, ny, nz, dx, y, dz, p, u, v, w, P, UF, VF, WF)
  D   += MolecularDiffusion(nx, ny, nz, dx, y, dz, Rgas, rho, u, v, w, p, U, V, W)
  eps += Dissipation(nx, ny, nz, dx, y, dz, Rgas, rho, u, v, w, p, Rho, U, V, W, P, UF, VF, WF)


@njit(cache=True, fastmath=True, nogil=True)
def Production(nx, ny, nz, dx, y, dz, rho, u, v, w, U, V, W):
  # rho, u, v, w : instantaneus
  # U, V, W      : Favre average
  # Favere decompose
  uF = u - U
  vF = v - V
  wF = w - W
  P  = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  for k in range(1,nz-1):
    for j in range(1,ny-1):
      for i in range(1,nx-1):
        uu = uF[k,j,i]**2
        uv = uF[k,j,i] * vF[k,j,i]
        uw = uF[k,j,i] * wF[k,j,i]
        vv = vF[k,j,i]**2
        vw = vF[k,j,i] * wF[k,j,i]
        ww = wF[k,j,i]**2
        P[k-1,j-1,i-1] = rho[k,j,i] * ( \
                        + uu * (-U[k,j,i-1] + U[k,j,i+1]) / (2.e0 * dx) \
                        + uv * (-V[k,j,i-1] + V[k,j,i+1]) / (2.e0 * dx) \
                        + uw * (-W[k,j,i-1] + W[k,j,i+1]) / (2.e0 * dx) \
                        + uv * (-U[k,j-1,i] + U[k,j+1,i]) / (-y[j-1] + y[j+1]) \
                        + vv * (-V[k,j-1,i] + V[k,j+1,i]) / (-y[j-1] + y[j+1]) \
                        + vw * (-W[k,j-1,i] + W[k,j+1,i]) / (-y[j-1] + y[j+1]) \
                        + uw * (-U[k-1,j,i] + U[k+1,j,i]) / (2.e0 * dz) \
                        + vw * (-V[k-1,j,i] + V[k+1,j,i]) / (2.e0 * dz) \
                        + ww * (-W[k-1,j,i] + W[k+1,j,i]) / (2.e0 * dz))
  return -P


@njit(cache=True, fastmath=True, nogil=True)
def Transport(nx, ny, nz, dx, y, dz, rho, u, v, w, U, V, W):
  # rho, u, v, w : instantaneus
  # U, V, W      : Favre average
  uK = np.zeros((nz-2,ny-2,nx-1), dtype=np.float32)
  vK = np.zeros((nz-2,ny-1,nx-2), dtype=np.float32)
  wK = np.zeros((nz-1,ny-2,nx-2), dtype=np.float32)
  T  = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  K  = 0.5e0 * ((u - U)**2 + (v - V)**2 + (w - W)**2)
  for k in range(1,nz-1):
    for j in range(1,ny-1):
      for i in range(nx-1):
        uf = u[k,j,i] - U[k,j,i] + u[k,j,i+1] - U[k,j,i+1]
        uK[k-1,j-1,i] = 0.125e0 * (rho[k,j,i] + rho[k,j,i+1]) \
                        * uf * (K[k,j,i] + K[k,j,i+1])
  
  for k in range(1,nz-1):
    for j in range(ny-1):
      for i in range(1,nx-1):
        vf = v[k,j,i] - V[k,j,i] + v[k,j+1,i] - V[k,j+1,i]
        vK[k-1,j,i-1] = 0.125e0 * (rho[k,j,i] + rho[k,j+1,i]) \
                        * vf * (K[k,j,i] + K[k,j+1,i])

  for k in range(nz-1):
    for j in range(1,ny-1):
      for i in range(1,nx-1):
        wf = w[k,j,i] - W[k,j,i] + w[k+1,j,i] - W[k+1,j,i]
        uK[k-1,j-1,i] = 0.125e0 * (rho[k,j,i] + rho[k+1,j,i]) \
                        * wf * (K[k,j,i] + K[k,j,i+1])

  for k in range(nz-2):
    for j in range(ny-2):
      for i in range(nx-2):
        T[k,j,i] = (-uK[k,j,i] + uK[k,j,i+1]) / dx \
                 + (-vK[k,j,i] + vK[k,j+1,i]) / (-y[j] + y[j+1]) \
                 + (-wK[k,j,i] + wK[k+1,j,i]) / dz
  return -T


@njit(cache=True, fastmath=True, nogil=True)
def Pressure(nx, ny, nz, dx, y, dz, p, u, v, w, P, U, V, W):
  # p, u, v, w : instantaneus
  # P          : Reynolds average
  # U, V, W    : Favre average
  pu  = np.zeros((nz-2,ny-2,nx-1), dtype=np.float32)
  pv  = np.zeros((nz-2,ny-1,nx-2), dtype=np.float32)
  pw  = np.zeros((nz-1,ny-2,nx-2), dtype=np.float32)
  PDl = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  PDf = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  for k in range(1,nz-1):
    for j in range(1,ny-1):
      for i in range(nx-1):#centered
        pf = p[k,j,i] - P[k,j,i] + p[k,j,i+1] - P[k,j,i+1]
        uf = u[k,j,i] - U[k,j,i] + u[k,j,i+1] - U[k,j,i+1]
        # Pressure Dilatation
        pu[k-1,j-1,i] = 0.25e0 * pf * uf

  for k in range(1,nz-1):
    for j in range(ny-1):#centered
      for i in range(1,nx-1):
        pf = p[k,j,i] - P[k,j,i] + p[k,j+1,i] - P[k,j+1,i]
        vf = v[k,j,i] - V[k,j,i] + v[k,j+1,i] - V[k,j+1,i]
        # Pressure Dilatation
        pv[k-1,j,i-1] = 0.25e0 * pf * vf

  for k in range(nz-1):#centered
    for j in range(1,ny-1):
      for i in range(1,nx-1):
        pf = p[k,j,i] - P[k,j,i] + p[k+1,j,i] - P[k+1,j,i]
        wf = w[k,j,i] - W[k,j,i] + w[k+1,j,i] - W[k+1,j,i]
        # Pressure Dilatation
        pw[k,j-1,i-1] = 0.25e0 * pf * wf

  for k in range(nz-2):
    for j in range(ny-2):
      for i in range(nx-2):
        # Pressure Dilatation
        PDl[k,j,i] = (-pu[k,j,i] + pu[k,j,i+1]) / dx \
                   + (-pv[k,j,i] + pv[k,j+1,i]) / (-y[j] + y[j+1]) \
                   + (-pw[k,j,i] + pw[k+1,j,i]) / dz
        # Pressure Diffusion
        pf = p[k+1,j+1,i+1] - P[k+1,j+1,i+1]
        PDf[k,j,i] = pf * ((-(u[k+1,j+1,i] + U[k+1,j+1,i]) + (u[k+1,j+1,i+2] - U[k+1,j+1,i+2])) / (2.e0 * dx) \
                         + (-(v[k+1,j,i+1] + V[k+1,j,i+1]) + (v[k+1,j+2,i+1] - V[k+1,j+2,i+1])) / (-y[j] + y[j+2]) \
                         + (-(w[k,j+1,i+1] + W[k,j+1,i+1]) + (w[k+2,j+1,i+1] - W[k+2,j+1,i+1])) / (2.e0 * dz))
  return -PDl + PDf


@njit(cache=True, fastmath=True, nogil=True)
def MolecularDiffusion(nx, ny, nz, dx, y, dz, Rgas, rho, u, v, w, p, U, V, W):
  # rho, u, v, w, p : instantaneus
  # U, V, W         : Favre average
  txx, tyx, tzx, txy, tyy, tzy, txz, tyz, tzz = tau(nx, ny, nz, dx, y, dz, Rgas, rho, u, v, w, p)
  Txx, Tyx, Tzx, Txy, Tyy, Tzy, Txz, Tyz, Tzz = tau(nx, ny, nz, dx, y, dz, Rgas, rho, U, V, W, p)
  txx = txx - Txx
  tyx = tyx - Tyx
  Tzx = tzx - Tzx
  txy = txy - Txy
  tyy = tyy - Tyy
  tzy = tzy - Tzy
  txz = txz - Txz
  tyz = tyz - Tyz
  tzz = tzz - Tzz
  Ev  = np.zeros((nz-2,ny-2,nx-1), dtype=np.float32)
  Fv  = np.zeros((nz-2,ny-1,nx-2), dtype=np.float32)
  Gv  = np.zeros((nz-1,ny-2,nx-2), dtype=np.float32)
  D   = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  # Favere decompose
  uF = u - U
  vF = v - V
  wF = w - W
  for k in range(1,nz-1):
    for j in range(1,ny-1):
      for i in range(nx-1):#centered
        # Molecular Diffusion
        utxx = 0.5e0 * (uF[k,j,i] + uF[k,j,i+1]) * txx[k-1,j-1,i]
        vtyx = 0.5e0 * (vF[k,j,i] + vF[k,j,i+1]) * tyx[k-1,j-1,i]
        wtzx = 0.5e0 * (wF[k,j,i] + wF[k,j,i+1]) * tzx[k-1,j-1,i]
        Ev[k-1,j-1,i] = utxx + vtyx + wtzx

  for k in range(1,nz-1):
    for j in range(ny-1):#centered
      for i in range(1,nx-1):
        # Molecular Diffusion
        utxy = 0.5e0 * (uF[k,j,i] + uF[k,j+1,i]) * txy[k-1,j,i-1]
        vtyy = 0.5e0 * (vF[k,j,i] + vF[k,j+1,i]) * tyy[k-1,j,i-1]
        wtzy = 0.5e0 * (wF[k,j,i] + wF[k,j+1,i]) * tzy[k-1,j,i-1]
        Fv[k-1,j,i-1] = utxy + vtyy + wtzy
  
  for k in range(nz-1):#centered
    for j in range(1,ny-1):
      for i in range(1,nx-1):
        # Molecular Diffusion
        utxz = 0.5e0 * (uF[k,j,i] + uF[k+1,j,i]) * txz[k,j-1,i-1]
        vtyz = 0.5e0 * (vF[k,j,i] + vF[k+1,j,i]) * tyz[k,j-1,i-1]
        wtzz = 0.5e0 * (wF[k,j,i] + wF[k+1,j,i]) * tzz[k,j-1,i-1]
        Gv[k,j-1,i-1] = utxz + vtyz + wtzz
  
  for k in range(nz-2):
    for j in range(ny-2):
      for i in range(nx-2):
        # Molecular Diffusion
        D[k,j,i] = (-Ev[k,j,i] + Ev[k,j,i+1]) / dx \
                 + (-Fv[k,j,i] + Fv[k,j+1,i]) / (-y[j] + y[j+1]) \
                 + (-Gv[k,j,i] + Gv[k+1,j,i]) / dz
  return D


@njit(cache=True, fastmath=True, nogil=True)
def Dissipation(nx, ny, nz, dx, y, dz, Rgas, rho, u, v, w, p, Rho, U, V, W, P, UF, VF, WF):
  # rho, u, v, w, p : instantaneus
  # Rho, U, V, W, P : Reynolds average
  # UF, VF, WF      : Favre average
  txx, tyx, tzx, txy, tyy, tzy, txz, tyz, tzz = tau(nx, ny, nz, dx, y, dz, Rgas, rho, u, v, w, p)
  Txx, Tyx, Tzx, Txy, Tyy, Tzy, Txz, Tyz, Tzz = tau(nx, ny, nz, dx, y, dz, Rgas, Rho, U, V, W, P)
  txx = txx - Txx
  tyx = tyx - Tyx
  Tzx = tzx - Tzx
  txy = txy - Txy
  tyy = tyy - Tyy
  tzy = tzy - Tzy
  txz = txz - Txz
  tyz = tyz - Tyz
  tzz = tzz - Tzz
  eps = np.zeros((nz-2,ny-2,nx-2), dtype=np.float32)
  # Favere decompose
  uF = u - UF
  vF = v - VF
  wF = w - WF
  for k in range(nz-2):
    for j in range(ny-2):
      for i in range(nx-2):
        # Dissipation
        txxc = 0.5e0 * (txx[k,j,i] + txx[k,j,i+1])
        tyxc = 0.5e0 * (tyx[k,j,i] + tyx[k,j,i+1])
        tzxc = 0.5e0 * (tzx[k,j,i] + tzx[k,j,i+1])
        txyc = 0.5e0 * (txy[k,j,i] + txy[k,j+1,i])
        tyyc = 0.5e0 * (tyy[k,j,i] + tyy[k,j+1,i])
        tzyc = 0.5e0 * (tzy[k,j,i] + tzy[k,j+1,i])
        txzc = 0.5e0 * (txz[k,j,i] + txz[k+1,j,i])
        tyzc = 0.5e0 * (tyz[k,j,i] + tyz[k+1,j,i])
        tzzc = 0.5e0 * (tzz[k,j,i] + tzz[k+1,j,i])
        eps[k,j,i] = txxc * (-uF[k+1,j+1,i] + uF[k+1,j+1,i+2]) / (2.e0 * dx) \
                   + tyxc * (-vF[k+1,j+1,i] + vF[k+1,j+1,i+2]) / (2.e0 * dx) \
                   + tzxc * (-wF[k+1,j+1,i] + wF[k+1,j+1,i+2]) / (2.e0 * dx) \
                   + txyc * (-uF[k+1,j,i+1] + uF[k+1,j+2,i+1]) / (-y[j] + y[j+2]) \
                   + tyyc * (-vF[k+1,j,i+1] + vF[k+1,j+2,i+1]) / (-y[j] + y[j+2]) \
                   + tzyc * (-wF[k+1,j,i+1] + wF[k+1,j+2,i+1]) / (-y[j] + y[j+2]) \
                   + txzc * (-uF[k,j+1,i+1] + uF[k+2,j+1,i+1]) / (2.e0 * dz) \
                   + tyzc * (-vF[k,j+1,i+1] + vF[k+2,j+1,i+1]) / (2.e0 * dz) \
                   + tzzc * (-wF[k,j+1,i+1] + wF[k+2,j+1,i+1]) / (2.e0 * dz)
  return eps


def print_TKE(x, y, z, P, T, PI, D, eps, directory, name):
  P1d   = np.float32(P.flatten())
  T1d   = np.float32(T.flatten())
  PI1d  = np.float32(PI.flatten())
  D1d   = np.float32(D.flatten())
  eps1d = np.float32(eps.flatten())

  os.makedirs(directory, exist_ok=True)
  filename  = name + ".vtr" 
  filepath  = os.path.join(directory, filename)

  x_coords = vtk.vtkFloatArray()
  y_coords = vtk.vtkFloatArray()
  z_coords = vtk.vtkFloatArray()
  x_coords.SetName("X-Axis")
  y_coords.SetName("Y-Axis")
  z_coords.SetName("Z-Axis")

  nx = len(x)
  ny = len(y)
  nz = len(z)
  
  for i in range(nx):
    x_coords.InsertNextValue(x[i])
  for j in range(ny):
    y_coords.InsertNextValue(y[j])
  for k in range(nz):
    z_coords.InsertNextValue(z[k])
  
  grid = vtk.vtkRectilinearGrid()
  grid.SetDimensions(nx, ny, nz)
  grid.SetXCoordinates(x_coords)
  grid.SetYCoordinates(y_coords)
  grid.SetZCoordinates(z_coords)

  P = vtk.vtkFloatArray()
  P.SetName("P")
  for i in range(nx * ny * nz):
    P.InsertNextValue(P1d[i])
  grid.GetPointData().AddArray(P)

  T = vtk.vtkFloatArray()
  T.SetName("T")
  for i in range(nx * ny * nz):
    T.InsertNextValue(T1d[i])
  grid.GetPointData().AddArray(T)

  PI = vtk.vtkFloatArray()
  PI.SetName("PI")
  for i in range(nx * ny * nz):
    PI.InsertNextValue(PI1d[i])
  grid.GetPointData().AddArray(PI)

  D = vtk.vtkFloatArray()
  D.SetName("D")
  for i in range(nx * ny * nz):
    D.InsertNextValue(D1d[i])
  grid.GetPointData().AddArray(D)

  eps = vtk.vtkFloatArray()
  eps.SetName("eps")
  for i in range(nx * ny * nz):
    eps.InsertNextValue(eps1d[i])
  grid.GetPointData().AddArray(eps)

  writer = vtk.vtkXMLRectilinearGridWriter()
  writer.SetFileName(filepath)
  writer.SetInputData(grid)
  writer.SetDataModeToAppended()
  writer.EncodeAppendedDataOff()
  writer.Write()

