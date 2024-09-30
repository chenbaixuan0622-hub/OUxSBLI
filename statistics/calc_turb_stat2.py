import numpy as np
import os
import re
from tqdm import tqdm
from numba import jit
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl, Sutherland

Q_directory = "../../../../../mnt/data1/TBL_5thHRSLAU2_20240922"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]
num_files = len(Q_files)

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

mean_path  = os.path.join(Q_directory, "Qmean.vtr")
rhom       = getScalar(mean_path, Nx, Ny, Nz, 'rho')
um, vm, wm = getVector(mean_path, Nx, Ny, Nz, 'velocity')
pm         = getScalar(mean_path, Nx, Ny, Nz, 'p')

Q  = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)

umF = np.load()
vmF = np.load()
wmF = np.load()
TmF = np.load()
TtF = np.load()

uT = np.zeros(Ny-2, dtype=np.float32)
vT = np.zeros(Ny-2, dtype=np.float32)
uu = np.zeros(Ny-2, dtype=np.float32)
vv = np.zeros(Ny-2, dtype=np.float32)
T2 = np.zeros(Ny-2, dtype=np.float32)

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
UF = np.zeros(Ny-2, dtype=np.float32)
VF = np.zeros(Ny-2, dtype=np.float32)
WF = np.zeros(Ny-2, dtype=np.float32)

yp, _, _, tw, ut, _ = non_dim_tbl(Q, x, y, z)

# turbulent Pr
dTtdT = np.mean((-TtF[:,:-2,:] + TtF[:,2:,:]) / (-TmF[:,:-2,:] + TmF[:,2:,:]), axis=(0,2))


@jit(nopython=True, cached=True, fastmath=True)
def tau(dx, y, dz, Rgas, rho, u, v, w, p):
  mu  = 1.716e-5 * ((273.2e0 + 111.e0) / (p / (rho * Rgas) + 111.e0)) * ((p / (rho * Rgas)) / 273.2e0)**1.5
  mux = 0.5e0 * (mu[1:-1,1:-1,:-1] + mu[1:-1,1:-1,1:])
  muy = 0.5e0 * (mu[1:-1,:-1,1:-1] + mu[1:-1,1:,1:-1])
  muz = 0.5e0 * (mu[:-1,1:-1,1:-1] + mu[1:,1:-1,1:-1])
  uy = 0.50e0 * (-u[1:-1,:-2,:-1] + u[1:-1,2:,:-1] - u[1:-1,:-2,1:] + u[1:-1,2:,1:]) / (-y[:-2] + y[2:])
  vy = 0.50e0 * (-v[1:-1,:-2,:-1] + v[1:-1,2:,:-1] - v[1:-1,:-2,1:] + v[1:-1,2:,1:]) / (-y[:-2] + y[2:])
  wz = 0.25e0 * (-w[:-2,1:-1,:-1] + w[2:,1:-1,:-1] - w[:-2,1:-1,1:] + w[2:,1:-1,1:]) / dz
  uz = 0.25e0 * (-u[:-2,1:-1,:-1] + u[2:,1:-1,:-1] - u[:-2,1:-1,1:] + u[2:,1:-1,1:]) / dz
  txx = 2.e0 * mux * (2.e0 * (-u[1:-1,1:-1,:-1] + u[1:-1,1:-1,1:]) / dx- vy - wz) / 3.e0
  tyx = mux * ((-v[1:-1,1:-1,:-1] + v[1:-1,1:-1,1:]) / dx + uy)
  tzx = mux * ((-w[1:-1,1:-1,:-1] + w[1:-1,1:-1,1:]) / dx + uz)
  vz = 0.25e0 * (-v[:-2,:-1,1:-1] + v[2:,:-1,1:-1] - v[:-2,1:,1:-1] + v[2:,1:,1:-1]) / dz
  wz = 0.25e0 * (-w[:-2,:-1,1:-1] + w[2:,:-1,1:-1] - w[:-2,1:,1:-1] + w[2:,1:,1:-1]) / dz
  ux = 0.25e0 * (-u[1:-1,:-1,:-2] + u[1:-1,:-1,2:] - u[1:-1,1:,:-2] + u[1:-1,1:,2:]) / dx
  vx = 0.25e0 * (-v[1:-1,:-1,:-2] + v[1:-1,:-1,2:] - v[1:-1,1:,:-2] + v[1:-1,1:,2:]) / dx
  txy = muy * ((-u[1:-1,:-1,1:-1] + u[1:-1,1:,1:-1]) / (-y[:-1] + y[1:])+ vx)
  tyy = 2.e0 * muy * (2.e0 * (-v[1:-1,:-1,1:-1] + v[1:-1,1:,1:-1]) / (-y[:-1] + y[1:])- wz - ux) / 3.e0
  tzy = muy * ((-w[1:-1,:-1,1:-1] + w[1:-1,1:,1:-1]) / (-y[:-1] + y[1:])+ vz)
  wx = 0.25e0 * (-w[:-1,1:-1,:-2] + w[:-1,1:-1,2:] - w[1:,1:-1,:-2] + w[1:,1:-1,2:]) / dx
  ux = 0.25e0 * (-u[:-1,1:-1,:-2] + u[:-1,1:-1,2:] - u[1:,1:-1,:-2] + u[1:,1:-1,2:]) / dx
  vy = 0.50e0 * (-v[:-1,:-2,1:-1] + v[:-1,2:,1:-1] - v[1:,:-2,1:-1] + v[1:,2:,1:-1]) / (-y[:-2] + y[2:])
  wy = 0.50e0 * (-w[:-1,:-2,1:-1] + w[:-1,2:,1:-1] - w[1:,:-2,1:-1] + w[1:,2:,1:-1]) / (-y[:-2] + y[2:])
  txz = muz * ((-u[:-1,1:-1,1:-1] + u[1:,1:-1,1:-1]) / dz + wx)
  tyz = muz * ((-v[:-1,1:-1,1:-1] + v[1:,1:-1,1:-1]) / dz + wy)
  tzz = 2.e0 * muz * (2.e0 * (-w[:-1,1:-1,1:-1] + w[1:,1:-1,1:-1]) / dz - ux - vy) / 3.e0
  return txx, tyx, tzx, txy, tyy, tzy, txz, tyz, tzz


Rgas = 287.03
txx, tyx, tzx, txy, tyy, tzy, txz, tyz, tzz = tau(-x[0] + x[1], y, -z[0] + z[1], Rgas, rhom, um, vm, wm, pm)
for Q_file in tqdm(Q_files):
  file_path = os.path.join(Q_directory, Q_file)
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
  # turbulent kinetic energy budget
  K  += 0.5e0 * rho * (uF**2 + vF**2 + wF**2)
  uv += np.mean(rho[:,1:-1,:] * uF[:,1:-1,:] * vF[:,1:-1,:], axis=(0,2))
  uw += np.mean(rho[:,1:-1,:] * uF[:,1:-1,:] * wF[:,1:-1,:], axis=(0,2))
  vw += np.mean(rho[:,1:-1,:] * vF[:,1:-1,:] * wF[:,1:-1,:], axis=(0,2))
  ww += np.mean(rho[:,1:-1,:] * wF[:,1:-1,:]**2, axis=(0,2))
  uK += np.mean(u[:,1:-1,:] * K[:,1:-1,:], axis=0)
  vK += np.mean(v * K,           axis=(0,2))
  wK += np.mean(w[:,1:-1,:] * K[:,1:-1,:], axis=-1)
  up += np.mean(uf[:,1:-1,:] * pf[:,1:-1,:], axis=0)
  vp += np.mean(vf           * pf[:,1:-1,:], axis=(0,2))
  wp += np.mean(wf[:,1:-1,:] * pf[:,1:-1,:], axis=2)
  UF += np.mean(uF[:,1:-1,:], axis=(0,2))
  VF += np.mean(vF[:,1:-1,:], axis=(0,2))
  WF += np.mean(wF[:,1:-1,:], axis=(0,2))
  txxf, tyxf, tzxf, txyf, tyyf, tzyf, txzf, tyzf, tzzf = tau(-x[0] + x[1], y, -z[0] + z[1], Rgas, rho, u, v, w, p)
  # D
  utxx += 0.5e0 * np.mean(uf[:,1:-1,:-1] + uf[:,1:-1,1:], axis=0)  * (txxf - txx)
  vtyx += 0.5e0 * np.mean(vf[:,1:-1,:-1] + vf[:,1:-1,1:], axis=0)  * (tyxf - tyx)
  wtzx += 0.5e0 * np.mean(wf[:,1:-1,:-1] + wf[:,1:-1,1:], axis=0)  * (tzxf - tzx)
  utxy += 0.5e0 * np.mean(uf[:,:-1,:]    + uf[:,1:,:], axis=(0,2)) * (txyf - txy)
  vtyy += 0.5e0 * np.mean(vf[:,:-1,:]    + vf[:,1:,:], axis=(0,2)) * (tyyf - tyy)
  wtzy += 0.5e0 * np.mean(wf[:,:-1,:]    + wf[:,1:,:], axis=(0,2)) * (tzyf - tzy)
  utxz += 0.5e0 * np.mean(uf[:-1,1:-1,:] + uf[1:,1:-1,:], axis=2)  * (txzf - txz)
  vtyz += 0.5e0 * np.mean(vf[:-1,1:-1,:] + vf[1:,1:-1,:], axis=2)  * (tyzf - tyz)
  wtzz += 0.5e0 * np.mean(wf[:-1,1:-1,:] + wf[1:,1:-1,:], axis=2)  * (tzzf - tzz)
  # eps
  eps += np.mean(txxf - txx, axis=(0,2)) * 0.5e0 * np.mean(-uf[:,1:-1,:-1] + uf[:,1:-1,1:], axis=(0,2)) / dx \
       + np.mean(tyxf - tyx, axis=(0,2)) * 0.5e0 * np.mean(-vf[:,1:-1,:-1] + vf[:,1:-1,1:], axis=(0,2)) / dx \
       + np.mean(tzxf - tzx, axis=(0,2)) * 0.5e0 * np.mean(-wf[:,1:-1,:-1] + wf[:,1:-1,1:], axis=(0,2)) / dx \
       + np.mean(txyf - txy, axis=(0,2)) * np.mean(-uf[:,:-1,:] + uf[:,1:,:], axis=(0,2)) / (-y[:-2] + y[2:]) \
       + np.mean(tyyf - tyy, axis=(0,2)) * np.mean(-vf[:,:-1,:] + vf[:,1:,:], axis=(0,2)) / (-y[:-2] + y[2:]) \
       + np.mean(tzyf - tzy, axis=(0,2)) * np.mean(-wf[:,:-1,:] + wf[:,1:,:], axis=(0,2)) / (-y[:-2] + y[2:]) \
       + np.mean(txzf - txz, axis=(0,2)) * 0.5e0 * np.mean(-uf[:-1,1:-1,:] + uf[1:,1:-1,:], axis=(0,2)) / dz \
       + np.mean(tyzf - tyz, axis=(0,2)) * 0.5e0 * np.mean(-vf[:-1,1:-1,:] + vf[1:,1:-1,:], axis=(0,2)) / dz \
       + np.mean(tzzf - tzz, axis=(0,2)) * 0.5e0 * np.mean(-wf[:-1,1:-1,:] + wf[1:,1:-1,:], axis=(0,2)) / dz

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
uK   /= (np.float32(num_files) * np.mean(rho[:,1:-1,:], axis=0)))
vK   /= (np.float32(num_files) * np.mean(rho,           axis=(0,2)))
wK   /= (np.float32(num_files) * np.mean(rho[:,1:-1,:], axis=-1)))
pu   /= np.float32(num_files)
pv   /= np.float32(num_files)
pw   /= np.float32(num_files)
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
eps  /= np.float32(num_files)

# thermal statistics
RuT = uT / (np.sqrt(uu) * np.sqrt(T2))
RvT = vT / (np.sqrt(vv) * np.sqrt(T2))
Prt = (1.e0 - vTt / vT) / (1.e0 - dTtdT)
# TKE budget
C   = np.mean((-rho[:,1:-1,:-2] * umF[:,1:-1,:-2] * K[:,1:-1,:-2] \
              + rho[:,1:-1,2:]  * umF[:,1:-1,2:]  * K[:,1:-1,2:]) / (2.e0 * dx),  axis=(0,2)) \
    + np.mean((-rho[:,:-2,:]    * vmF[:,:-2,:]    * K[:,:-2,:]    \
              + rho[:,2:,:] * vmF[:,2:,:] * K[:,2:,:]) / (-y[:-2] + y[2:]), axis=(0,2)) \
    + np.mean((-rho[:-2,1:-1,:] * wmF[:-2,1:-1,:] * K[:-2,1:-1,:] \
              + rho[2:,1:-1,:]  * wmF[2:,1:-1,:]  * K[2:,1:-1,:]) / (2.e0 * dz),  axis=(0,2))
P   = -np.mean(rho[:,1:-1,:], axis=(0,2)) * (uu + uv + uw) * np.mean(-umF[:,1:-1,:-2] + umF[:,1:-1,2:], axis=(0,2)) / (2.e0 * dx) \
      -np.mean(rho[:,1:-1,:], axis=(0,2)) * (uv + vv + vw) * np.mean(-vmF[:,:-2,:]    + vmF[:,2:,:],    axis=(0,2)) / (-y[:-2] + y[2:]) \
      -np.mean(rho[:,1:-1,:], axis=(0,2)) * (uw + vw + ww) * np.mean(-wmF[:-2,1:-1,:] + wmF[2:,1:-1,:], axis=(0,2)) / (2.e0 * dz)
T   = -np.mean((- np.mean(rho[:,1:-1,:-2],  axis=0) * uK[:,:-2] \
                + np.mean(rho[:,1:-1,2:],   axis=0) * uK[:,2:] - up[:,:-2] + up[:,2:]) / (2.e0 * dx), axis)\
      -       ((- np.mean(rho[:,:-2,:], axis=(0,2)) * vK[:-2] \
                + np.mean(rho[:,2:,:],  axis=(0,2)) * vK[2:]   - vp[:-2]   + vp[2:])   / (-y[:-2] + y[2:]))\
      -np.mean((- np.mean(rho[:-2,1:-1,:],  axis=2) * wK[:-2,:] \
                + np.mean(rho[2:,1:-1,:],   axis=2) * wK[2:,:] - wp[:-2,:] + wp[2:,:]) / (2.e0 * dz), axis=)
PI  = np.mean(pf[:,1:-1,:], axis=(0,2)) * (np.mean(-uf[:,1:-1,:-2] + uf[:,1:-1,2:], axis=(0,2)) / (2.e0 * dx) \
                                         + np.mean(-vf[:,:-2,:]    + vf[:,2:,:],    axis=(0,2)) / (-y[:-2] + y[2:]) \
                                         + np.mean(-wf[:-2,1:-1,:] + wf[2:,1:-1,:], axis=(0,2)) / (2.e0 * dz))
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
           - np.mean(-pm[:-2,1:-1,:] + pm[2:,1:-1,:], axis=(0,2)) / dz)
D   = (utxx + vtyx + wtzx) / dx \
      (utxy + vtyy + wtzy) / (-y[:-2] + y[2:]) \
      (utxz + vtyz + wtzz) / dz

save_path = os.path.join(Q_directory, "turb_stat2.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# yp     Rut      RvT      Prt      C       P       T       PI       M       D      eps", file=f)
  for j in range(1,Ny-1):
    print(f'{yp[j]:.3e}', f'{Rut[j]:.3e}', f'{RvT[j]:.3e}', f'{Prt[j]:.3e}', \
    f'{C[j]:.3e}', f'{P[j]:.3e}', f'{T[j]:.3e}', f'{PI[j]:.3e}', f'{M[j]:.3e}', f'{D[j]:.3e}', f'{eps[j]:.3e}', file=f)

