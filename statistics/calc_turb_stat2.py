import numpy as np
import os
import re
from tqdm import tqdm
from numba import jit
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl, Sutherland

Q_directory = "../3D_solver/TBL/data"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)
dx = -x[0] + x[1]
dz = -z[0] + z[1]

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

gamma = 1.4e0
Rgas = 287.03

@jit(nopython=True, cache=True, fastmath=True)
def tau(dx, y, dz, Rgas, rho, u, v, w, p):
  mu  = 1.716e-5 * ((273.2e0 + 111.e0) / (p / (rho * Rgas) + 111.e0)) * ((p / (rho * Rgas)) / 273.2e0)**1.5
  mux = 0.5e0 * (mu[1:-1,1:-1,:-1] + mu[1:-1,1:-1,1:])
  muy = 0.5e0 * (mu[1:-1,:-1,1:-1] + mu[1:-1,1:,1:-1])
  muz = 0.5e0 * (mu[:-1,1:-1,1:-1] + mu[1:,1:-1,1:-1])
  txx = np.zeros((Nz-2,Ny-2,Nx-1), dtype=np.float32)
  tyx = np.zeros((Nz-2,Ny-2,Nx-1), dtype=np.float32)
  tzx = np.zeros((Nz-2,Ny-2,Nx-1), dtype=np.float32)
  txy = np.zeros((Nz-2,Ny-1,Nx-2), dtype=np.float32)
  tyy = np.zeros((Nz-2,Ny-1,Nx-2), dtype=np.float32)
  tzy = np.zeros((Nz-2,Ny-1,Nx-2), dtype=np.float32)
  txz = np.zeros((Nz-1,Ny-2,Nx-2), dtype=np.float32)
  tyz = np.zeros((Nz-1,Ny-2,Nx-2), dtype=np.float32)
  tzz = np.zeros((Nz-1,Ny-2,Nx-2), dtype=np.float32)
  for k in range(1,Nz-1):
    for j in range(1,Ny-1):
      for i in range(Nx-1):
        uy = 0.50e0 * (-u[k,j-1,i] + u[k,j+1,i] - u[k,j-1,i+1] + u[k,j+1,i+1]) / (-y[j-1] + y[j+1])
        vy = 0.50e0 * (-v[k,j-1,i] + v[k,j+1,i] - v[k,j-1,i+1] + v[k,j+1,i+1]) / (-y[j-1] + y[j+1])
        wz = 0.25e0 * (-w[k-1,j,i] + w[k+1,j,i] - w[k-1,j,i+1] + w[k+1,j,i+1]) / dz
        uz = 0.25e0 * (-u[k-1,j,i] + u[k+1,j,i] - u[k-1,j,i+1] + u[k+1,j,i+1]) / dz
        txx[k-1,j-1,i] = 2.e0 * mux[k-1,j-1,i] * (2.e0 * (-u[k,j,i] + u[k,j,i+1]) / dx - vy - wz) / 3.e0
        tyx[k-1,j-1,i] = mux[k-1,j-1,i] * ((-v[k,j,i] + v[k,j,i+1]) / dx + uy)
        tzx[k-1,j-1,i] = mux[k-1,j-1,i] * ((-w[k,j,i] + w[k,j,i+1]) / dx + uz)
  for k in range(1,Nz-1):
    for j in range(Ny-1):
      for i in range(1,Nx-1):
        vz = 0.25e0 * (-v[k-1,j,i] + v[k+1,j,i] - v[k-1,j+1,i] + v[k+1,j+1,i]) / dz
        wz = 0.25e0 * (-w[k-1,j,i] + w[k+1,j,i] - w[k-1,j+1,i] + w[k+1,j+1,i]) / dz
        ux = 0.25e0 * (-u[k,j,i-1] + u[k,j,i+1] - u[k,j+1,i-1] + u[k,j+1,i+1]) / dx
        vx = 0.25e0 * (-v[k,j,i-1] + v[k,j,i+1] - v[k,j+1,i-1] + v[k,j+1,i+1]) / dx
        txy[k-1,j,i-1] = muy[k-1,j,i-1] * ((-u[k,j,i] + u[k,j+1,i]) / (-y[j] + y[j+1]) + vx)
        tyy[k-1,j,i-1] = 2.e0 * muy[k-1,j,i-1] * (2.e0 * (-v[k,j,i] + v[k,j+1,i]) / (-y[j] + y[j+1]) - wz - ux) / 3.e0
        tzy[k-1,j,i-1] = muy[k-1,j,i-1] * ((-w[k,j,i] + w[k,j+1,i]) / (-y[j] + y[j+1]) + vz)
  for k in range(Nz-1):
    for j in range(1,Ny-1):
      for i in range(1,Nx-1):
        wx = 0.25e0 * (-w[k,j,i-1] + w[k,j,i+1] - w[k+1,j,i-1] + w[k+1,j,i+1]) / dx
        ux = 0.25e0 * (-u[k,j,i-1] + u[k,j,i+1] - u[k+1,j,i-1] + u[k+1,j,i+1]) / dx
        vy = 0.50e0 * (-v[k,j-1,i] + v[k,j+1,i] - v[k+1,j-1,i] + v[k+1,j+1,i]) / (-y[j-1] + y[j+1])
        wy = 0.50e0 * (-w[k,j-1,i] + w[k,j+1,i] - w[k+1,j-1,i] + w[k+1,j+1,i]) / (-y[j-1] + y[j+1])
        txz[k,j-1,i-1] = muz[k,j-1,i-1] * ((-u[k,j,i] + u[k+1,j,i]) / dz + wx)
        tyz[k,j-1,i-1] = muz[k,j-1,i-1] * ((-v[k,j,i] + v[k+1,j,i]) / dz + wy)
        tzz[k,j-1,i-1] = 2.e0 * muz[k,j-1,i-1] * (2.e0 * (-w[k,j,i] + w[k+1,j,i]) / dz - ux - vy) / 3.e0
  return txx, tyx, tzx, txy, tyy, tzy, txz, tyz, tzz

def main():
  mean_path  = os.path.join(Q_directory, "Qmean.vtr")
  rhom       = getScalar(mean_path, Nx, Ny, Nz, 'rho')
  um, vm, wm = getVector(mean_path, Nx, Ny, Nz, 'velocity')
  pm         = getScalar(mean_path, Nx, Ny, Nz, 'p')
  
  umF_path = os.path.join(Q_directory, "uF.npy")
  vmF_path = os.path.join(Q_directory, "vF.npy")
  wmF_path = os.path.join(Q_directory, "wF.npy")
  TmF_path = os.path.join(Q_directory, "TF.npy")
  TtF_path = os.path.join(Q_directory, "TtF.npy")
  umF = np.load(umF_path)
  vmF = np.load(vmF_path)
  wmF = np.load(wmF_path)
  TmF = np.load(TmF_path)
  TtF = np.load(TtF_path)

  uT  = np.zeros(Ny-2, dtype=np.float32)
  vT  = np.zeros(Ny-2, dtype=np.float32)
  vTt = np.zeros(Ny-2, dtype=np.float32)
  uu  = np.zeros(Ny-2, dtype=np.float32)
  vv  = np.zeros(Ny-2, dtype=np.float32)
  T2  = np.zeros(Ny-2, dtype=np.float32)

  K  = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  uv = np.zeros(Ny-2, dtype=np.float32)
  uw = np.zeros(Ny-2, dtype=np.float32)
  vu = np.zeros(Ny-2, dtype=np.float32)
  vw = np.zeros(Ny-2, dtype=np.float32)
  wu = np.zeros(Ny-2, dtype=np.float32)
  wv = np.zeros(Ny-2, dtype=np.float32)
  ww = np.zeros(Ny-2, dtype=np.float32)
  uK = np.zeros((Ny-2,Nx), dtype=np.float32)
  vK = np.zeros(Ny,        dtype=np.float32)
  wK = np.zeros((Nz,Ny-2), dtype=np.float32)
  up = np.zeros((Ny-2,Nx), dtype=np.float32)
  vp = np.zeros(Ny,        dtype=np.float32)
  wp = np.zeros((Nz,Ny-2), dtype=np.float32)
  pdiv = np.zeros(Ny-2, dtype=np.float32) 
  UF = np.zeros(Ny-2, dtype=np.float32)
  VF = np.zeros(Ny-2, dtype=np.float32)
  WF = np.zeros(Ny-2, dtype=np.float32)
  utxx = np.zeros((Ny-2,Nx-1), dtype=np.float32)
  vtyx = np.zeros((Ny-2,Nx-1), dtype=np.float32)
  wtzx = np.zeros((Ny-2,Nx-1), dtype=np.float32)
  utxy = np.zeros((Ny-1),      dtype=np.float32)
  vtyy = np.zeros((Ny-1),      dtype=np.float32)
  wtzy = np.zeros((Ny-1),      dtype=np.float32)
  utxz = np.zeros((Nz-1,Ny-2), dtype=np.float32)
  vtyz = np.zeros((Nz-1,Ny-2), dtype=np.float32)
  wtzz = np.zeros((Nz-1,Ny-2), dtype=np.float32)
  PI   = np.zeros(Ny-2, dtype=np.float32)
  eps  = np.zeros(Ny-2, dtype=np.float32)

  Q = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
  Q[0,:,:,:], Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:], Q[4,:,:,:] = rhom, um, vm, wm, pm
  yp, _, _, _, ut, _ = non_dim_tbl(Q, x, y, z)
  del Q
  zt   = Sutherland(pm[:,0,:] / (Rgas * rhom[:,0,:])) / (np.mean(rhom[:,0,:]) * ut)
  Norm = np.mean(rhom[:,0,:]) * ut**3 / np.mean(zt)

  # turbulent Pr
  dTtdT = np.mean(-TtF[:,:-2,:] + TtF[:,2:,:], axis=(0,2)) / np.mean(-TmF[:,:-2,:] + TmF[:,2:,:], axis=(0,2))

  Y = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  for k in range(Nz):
    for j in range(Ny):
      for i in range(Nx):
        Y[k,j,i] = y[j]

  txx, tyx, tzx, txy, tyy, tzy, txz, tyz, tzz = tau(dx, y, dz, Rgas, rhom, um, vm, wm, pm)
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_directory, Q_file)
    if file_path == os.path.join(Q_directory, "Qmean.vtr"):
      continue
    rho     = getScalar(file_path, Nx, Ny, Nz, 'rho')
    u, v, w = getVector(file_path, Nx, Ny, Nz, 'velocity')
    p       = getScalar(file_path, Nx, Ny, Nz, 'p')
    uf  = u - um
    vf  = v - vm
    wf  = w - wm
    pf  = p - pm
    uF  = u - umF
    vF  = v - vmF
    wF  = w - wmF
    TF  = p / (rho * Rgas) - TmF
    Ttotal = p / (rho * Rgas) + 0.5e0 * (gamma - 1.e0) * (u**2 + v**2 + w**2) / gamma
    TtF = Ttotal - TtF
    # thermal statistics
    uT  += np.mean(rho[:,1:-1,:] * uF[:,1:-1,:] *  TF[:,1:-1,:], axis=(0,2))
    vT  += np.mean(rho[:,1:-1,:] * vF[:,1:-1,:] *  TF[:,1:-1,:], axis=(0,2))
    vTt += np.mean(rho[:,1:-1,:] * vF[:,1:-1,:] * TtF[:,1:-1,:], axis=(0,2))
    uu  += np.mean(rho[:,1:-1,:] * uF[:,1:-1,:]**2, axis=(0,2))
    vv  += np.mean(rho[:,1:-1,:] * vF[:,1:-1,:]**2, axis=(0,2))
    T2  += np.mean(rho[:,1:-1,:] * TF[:,1:-1,:]**2, axis=(0,2))
    # C
    K  += 0.5e0 * rho * (uF**2 + vF**2 + wF**2)
    # P
    uv += np.mean(rho[:,1:-1,:] * uF[:,1:-1,:] * vF[:,1:-1,:], axis=(0,2))
    uw += np.mean(rho[:,1:-1,:] * uF[:,1:-1,:] * wF[:,1:-1,:], axis=(0,2))
    vw += np.mean(rho[:,1:-1,:] * vF[:,1:-1,:] * wF[:,1:-1,:], axis=(0,2))
    ww += np.mean(rho[:,1:-1,:] * wF[:,1:-1,:]**2, axis=(0,2))
    # T
    uK += np.mean(0.5e0 * rho[:,1:-1,:] * uF[:,1:-1,:] * (uF[:,1:-1,:]**2 + vF[:,1:-1,:]**2 + wF[:,1:-1,:]**2), axis=0)
    vK += np.mean(0.5e0 * rho *           vF           * (uF**2           + vF**2           + wF**2),       axis=(0,2))
    wK += np.mean(0.5e0 * rho[:,1:-1,:] * wF[:,1:-1,:] * (uF[:,1:-1,:]**2 + vF[:,1:-1,:]**2 + wF[:,1:-1,:]**2), axis=2)
    # PI
    up   += np.mean(uF[:,1:-1,:] * pf[:,1:-1,:], axis=0)
    vp   += np.mean(vF           * pf,       axis=(0,2))
    wp   += np.mean(wF[:,1:-1,:] * pf[:,1:-1,:], axis=2)
    pdiv += np.mean(pf[1:-1,1:-1,1:-1] * ((-uF[1:-1,1:-1,:-2] + uF[1:-1,1:-1,2:]) / (2.e0 * dx) \
                                        + (-vF[1:-1,:-2,1:-1] + vF[1:-1,2:,1:-1]) / (-Y[1:-1,:-2,1:-1] + Y[1:-1,2:,1:-1]) \
                                        + (-wF[:-2,1:-1,1:-1] + wF[2:,1:-1,1:-1]) / (2.e0 * dz)), axis=(0,2))
    # M
    UF += np.mean(uF[:,1:-1,:], axis=(0,2))
    VF += np.mean(vF[:,1:-1,:], axis=(0,2))
    WF += np.mean(wF[:,1:-1,:], axis=(0,2))
    txxf, tyxf, tzxf, txyf, tyyf, tzyf, txzf, tyzf, tzzf = tau(dx, y, dz, Rgas, rho, u, v, w, p)
    # D
    utxx += 0.5e0 * np.mean((uf[1:-1,1:-1,:-1] + uf[1:-1,1:-1,1:]) * (txxf - txx), axis=0)
    vtyx += 0.5e0 * np.mean((vf[1:-1,1:-1,:-1] + vf[1:-1,1:-1,1:]) * (tyxf - tyx), axis=0)
    wtzx += 0.5e0 * np.mean((wf[1:-1,1:-1,:-1] + wf[1:-1,1:-1,1:]) * (tzxf - tzx), axis=0)
    utxy += 0.5e0 * np.mean((uf[1:-1,:-1,1:-1] + uf[1:-1,1:,1:-1]) * (txyf - txy), axis=(0,2))
    vtyy += 0.5e0 * np.mean((vf[1:-1,:-1,1:-1] + vf[1:-1,1:,1:-1]) * (tyyf - tyy), axis=(0,2))
    wtzy += 0.5e0 * np.mean((wf[1:-1,:-1,1:-1] + wf[1:-1,1:,1:-1]) * (tzyf - tzy), axis=(0,2))
    utxz += 0.5e0 * np.mean((uf[:-1,1:-1,1:-1] + uf[1:,1:-1,1:-1]) * (txzf - txz), axis=2)
    vtyz += 0.5e0 * np.mean((vf[:-1,1:-1,1:-1] + vf[1:,1:-1,1:-1]) * (tyzf - tyz), axis=2)
    wtzz += 0.5e0 * np.mean((wf[:-1,1:-1,1:-1] + wf[1:,1:-1,1:-1]) * (tzzf - tzz), axis=2)
    # eps
    eps += np.mean((txxf - txx) * 0.5e0 * (-uf[1:-1,1:-1,:-1] + uf[1:-1,1:-1,1:]), axis=(0,2)) / dx \
         + np.mean((tyxf - tyx) * 0.5e0 * (-vf[1:-1,1:-1,:-1] + vf[1:-1,1:-1,1:]), axis=(0,2)) / dx \
         + np.mean((tzxf - tzx) * 0.5e0 * (-wf[1:-1,1:-1,:-1] + wf[1:-1,1:-1,1:]), axis=(0,2)) / dx \
         + np.mean(0.5e0 * (txyf[:,:-1,:] + txyf[:,:-1,:] - txy[:,1:,:] - txy[:,1:,:]) \
         * (-uf[1:-1,:-2,1:-1] + uf[1:-1,2:,1:-1]), axis=(0,2)) / (-y[:-2] + y[2:]) \
         + np.mean(0.5e0 * (tyyf[:,:-1,:] + tyyf[:,:-1,:] - tyy[:,1:,:] - tyy[:,1:,:]) \
         * (-vf[1:-1,:-2,1:-1] + vf[1:-1,2:,1:-1]), axis=(0,2)) / (-y[:-2] + y[2:]) \
         + np.mean(0.5e0 * (tzyf[:,:-1,:] + tzyf[:,:-1,:] - tzy[:,1:,:] - tzy[:,1:,:]) \
         * (-wf[1:-1,:-2,1:-1] + wf[1:-1,2:,1:-1]), axis=(0,2)) / (-y[:-2] + y[2:]) \
         + np.mean((txzf - txz) * 0.5e0 * (-uf[:-1,1:-1,1:-1] + uf[1:,1:-1,1:-1]), axis=(0,2)) / dz \
         + np.mean((tyzf - tyz) * 0.5e0 * (-vf[:-1,1:-1,1:-1] + vf[1:,1:-1,1:-1]), axis=(0,2)) / dz \
         + np.mean((tzzf - tzz) * 0.5e0 * (-wf[:-1,1:-1,1:-1] + wf[1:,1:-1,1:-1]), axis=(0,2)) / dz

  uT   /= (np.float32(num_files) * np.mean(rho[:,1:-1,:], axis=(0,2)))
  vT   /= (np.float32(num_files) * np.mean(rho[:,1:-1,:], axis=(0,2)))
  vTt  /= (np.float32(num_files) * np.mean(rho[:,1:-1,:], axis=(0,2)))
  uu   /= (np.float32(num_files) * np.mean(rho[:,1:-1,:], axis=(0,2)))
  vv   /= (np.float32(num_files) * np.mean(rho[:,1:-1,:], axis=(0,2)))
  T2   /= (np.float32(num_files) * np.mean(rho[:,1:-1,:], axis=(0,2)))
  K    /= (np.float32(num_files) * rho)
  uv   /= (np.float32(num_files) * np.mean(rho[:,1:-1,:], axis=(0,2)))
  uw   /= (np.float32(num_files) * np.mean(rho[:,1:-1,:], axis=(0,2)))
  vw   /= (np.float32(num_files) * np.mean(rho[:,1:-1,:], axis=(0,2)))
  ww   /= (np.float32(num_files) * np.mean(rho[:,1:-1,:], axis=(0,2)))
  uK   /= (np.float32(num_files))
  vK   /= (np.float32(num_files))
  wK   /= (np.float32(num_files))
  up   /= np.float32(num_files)
  vp   /= np.float32(num_files)
  wp   /= np.float32(num_files)
  pdiv /= np.float32(num_files)
  UF   /= np.float32(num_files)
  VF   /= np.float32(num_files)
  WF   /= np.float32(num_files)
  utxx /= np.float32(num_files)
  vtyx /= np.float32(num_files)
  wtzx /= np.float32(num_files)
  utxy /= np.float32(num_files)
  vtyy /= np.float32(num_files)
  wtzy /= np.float32(num_files)
  utxz /= np.float32(num_files)
  vtyz /= np.float32(num_files)
  wtzz /= np.float32(num_files)
  PI   /= np.float32(num_files)
  eps  /= np.float32(num_files)

  # thermal statistics
  RuT = uT / (np.sqrt(uu) * np.sqrt(T2))
  RvT = vT / (np.sqrt(vv) * np.sqrt(T2))
  Prt = (1.e0 - vTt / vT) / (1.e0 - dTtdT)
  # TKE budget
  C   = np.mean((-rho[:,1:-1,:-2] * umF[:,1:-1,:-2] * K[:,1:-1,:-2] \
                + rho[:,1:-1,2:]  * umF[:,1:-1,2:]  * K[:,1:-1,2:]) / (2.e0 * dx),  axis=(0,2)) \
      + np.mean((-rho[:,:-2,:]    * vmF[:,:-2,:]    * K[:,:-2,:]    \
                + rho[:,2:,:] * vmF[:,2:,:] * K[:,2:,:]), axis=(0,2)) / (-y[:-2] + y[2:]) \
      + np.mean((-rho[:-2,1:-1,:] * wmF[:-2,1:-1,:] * K[:-2,1:-1,:] \
                + rho[2:,1:-1,:]  * wmF[2:,1:-1,:]  * K[2:,1:-1,:]) / (2.e0 * dz),  axis=(0,2))
  P   = -np.mean(rho[:,1:-1,:], axis=(0,2)) * (\
         + uu * np.mean(-umF[:,1:-1,:-2] + umF[:,1:-1,2:], axis=(0,2)) / (2.e0 * dx) \
         + uv * np.mean(-vmF[:,1:-1,:-2] + vmF[:,1:-1,2:], axis=(0,2)) / (2.e0 * dx) \
         + wu * np.mean(-wmF[:,1:-1,:-2] + wmF[:,1:-1,2:], axis=(0,2)) / (2.e0 * dx) \
         + uv * np.mean(-umF[:,:-2,:] + umF[:,2:,:], axis=(0,2)) / (-y[:-2] + y[2:]) \
         + vv * np.mean(-vmF[:,:-2,:] + vmF[:,2:,:], axis=(0,2)) / (-y[:-2] + y[2:]) \
         + vw * np.mean(-wmF[:,:-2,:] + wmF[:,2:,:], axis=(0,2)) / (-y[:-2] + y[2:]) \
         + wu * np.mean(-umF[:-2,1:-1,:] + umF[2:,1:-1,:], axis=(0,2)) / (2.e0 * dz) \
         + vw * np.mean(-vmF[:-2,1:-1,:] + vmF[2:,1:-1,:], axis=(0,2)) / (2.e0 * dz) \
         + ww * np.mean(-wmF[:-2,1:-1,:] + wmF[2:,1:-1,:], axis=(0,2)) / (2.e0 * dz) )
  T   = -np.mean((- uK[:,:-2] + uK[:,2:]) / (2.e0 * dx), axis=1)\
        -       ((- vK[:-2]   + vK[2:])   / (-y[:-2] + y[2:]))\
        -np.mean((- wK[:-2,:] + wK[2:,:]) / (2.e0 * dz), axis=0)
  PI  = -np.mean((-up[:,:-2] + up[:,2:]) / (2.e0 * dx), axis=1) \
        -        (-vp[:-2]   + vp[2:]) / (-y[:-2] + y[2:]) \
        -np.mean((-wp[:-2,:] + wp[2:,:]) / (2.e0 * dz), axis=0) + pdiv
  M   = np.mean(rhom[:,1:-1,:], axis=(0,2)) * \
        (UF * (np.mean(-txx[:,:,:-1] + txx[:,:,1:], axis=(0,2)) / dx \
             + np.mean(-txy[:,:-1,:] + txy[:,1:,:], axis=(0,2)) / (-y[:-2] + y[2:]) \
             + np.mean(-txz[:-1,:,:] + txz[1:,:,:], axis=(0,2)) / dz \
             - np.mean(-pm[:,1:-1,:-2] + pm[:,1:-1,2:], axis=(0,2)) / dx) \
       + VF * (np.mean(-tyx[:,:,:-1] + tyx[:,:,1:], axis=(0,2)) / dx \
             + np.mean(-tyy[:,:-1,:] + tyy[:,1:,:], axis=(0,2)) / (-y[:-2] + y[2:]) \
             + np.mean(-tyz[:-1,:,:] + tyz[1:,:,:], axis=(0,2)) / dz \
             - np.mean(-pm[:,:-2,:] + pm[:,2:,:], axis=(0,2)) / (-y[:-2] + y[2:])) \
       + WF * (np.mean(-tzx[:,:,:-1] + tzx[:,:,1:], axis=(0,2)) / dx \
             + np.mean(-tzy[:,:-1,:] + tzy[:,1:,:], axis=(0,2)) / (-y[:-2] + y[2:]) \
             + np.mean(-tzz[:-1,:,:] + tzz[1:,:,:], axis=(0,2)) / dz \
             - np.mean(-pm[:-2,1:-1,:] + pm[2:,1:-1,:], axis=(0,2)) / dz))
  D   = np.mean(-utxx[:,:-1] + utxx[:,1:] - vtyx[:,:-1] + vtyx[:,1:] - wtzx[:,:-1] + wtzx[:,1:], axis=1) / dx \
      + (-utxy[:-1] + utxy[1:] - vtyy[:-1] + vtyy[1:] - wtzy[:-1] + wtzy[1:]) / (-y[:-2] + y[2:]) \
      + np.mean(-utxz[:-1,:] + utxz[1:,:] - vtyz[:-1,:] + vtyz[1:,:] - wtzz[:-1,:] + wtzz[1:,:], axis=0) / dz

  # normalize
  C   = C   / Norm
  P   = P   / Norm
  T   = T   / Norm
  PI  = PI  / Norm
  M   = M   / Norm
  D   = D   / Norm
  eps = eps / Norm

  delta = 2.e-3
  for j in range(Ny):
    if np.mean(um[:,j,:]) >= 0.99e0 * np.mean(um[:,-1,:]):
      delta = y[j]
      break

  save_path = os.path.join(Q_directory, "turb_stat2.d")
  with open(save_path, "w", encoding="UTF-8") as f:
    print("# yp      y         Rut        RvT        Prt        C         P          T         PI         M          D      eps", file=f)
    for j in range(1,Ny-1):
      print(f'{yp[j]:.3e}', f'{y[j-1]/delta:.3e}', f'{RuT[j-1]:.3e}', f'{RvT[j-1]:.3e}', f'{Prt[j-1]:.3e}', \
      f'{C[j-1]:.3e}', f'{P[j-1]:.3e}', f'{T[j-1]:.3e}', f'{PI[j-1]:.3e}', f'{M[j-1]:.3e}', f'{D[j-1]:.3e}', f'{eps[j-1]:.3e}', file=f)

main()

