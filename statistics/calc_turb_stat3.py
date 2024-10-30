import numpy as np
import os
import re
from tqdm import tqdm
from numba import jit
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl, Sutherland, tau

Q_directory = "../../a100"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]

gamma = 1.4e0
Rgas  = 287.03

# target yplus
yp1 = 10.e0
yp2 = 80.e0

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, X, y, z = getGrid(first_path)
nx0 = int(0.1*Nx)
nx  = int(0.2*Nx)
x = X[nx0:nx0+nx]

def calc_corr(Nx, Ny, Nz, nx0, nx, x, y, z, yp1, yp2, rhom, um, vm, wm, pm, umF, vmF, wmF, TmF, TtmF):
  # Reynolds average
  # rhom, um, vm, wm, pm
  
  # Favre average
  # umF, vmF, wmF, TmF, TtF
  
  # x direction correlation
  uux = np.zeros((Nz,2,nx), dtype=np.float32)
  vvx = np.zeros((Nz,2,nx), dtype=np.float32)
  wwx = np.zeros((Nz,2,nx), dtype=np.float32)
  # z direction correlation
  uuz = np.zeros((Nz,2,nx), dtype=np.float32)
  vvz = np.zeros((Nz,2,nx), dtype=np.float32)
  wwz = np.zeros((Nz,2,nx), dtype=np.float32)

  # Strong Reynolds analogy
  uu   = np.zeros_like(um, dtype=np.float32)
  vv   = np.zeros_like(um, dtype=np.float32)
  uv   = np.zeros_like(um, dtype=np.float32)
  TT   = np.zeros_like(um, dtype=np.float32)
  TtTt = np.zeros_like(um, dtype=np.float32)
  uT   = np.zeros_like(um, dtype=np.float32)
  vT   = np.zeros_like(um, dtype=np.float32)
  vTt  = np.zeros_like(um, dtype=np.float32)

  # calc yp
  yp_path = os.path.join(Q_directory, 'yp.npy')
  yp = np.load(yp_path)
  for j in range(Ny):
    if yp[j] >= yp1:
      ny1 = j
      break

  for j in range(Ny):
    if yp[j] >= yp2:
      ny2 = j
      break


  # calc boundary layer thickness
  u0 = np.mean(um[:,-1,:])
  for j in range(Ny):
    if np.mean(um[:,j,:]) >= 0.99e0 * u0:
      delta = y[j] - (-y[j-1] + y[j]) * (np.mean(um[:,j,:]) - 0.99e0 * u0) / (-np.mean(um[:,j-1,:]) + np.mean(um[:,j,:]))
      break

  itr = 0.e0 
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_directory, Q_file)
    if file_path == os.path.join(Q_directory, "TKE.vtr") or \
       file_path == os.path.join(Q_directory, "Qmean.vtr"):
      continue
    Rho     = getScalar(file_path, Nx, Ny, Nz, 'rho')
    U, V, W = getVector(file_path, Nx, Ny, Nz, 'velocity')
    P       = getScalar(file_path, Nx, Ny, Nz, 'p')
    rho = Rho[:,:,nx0:nx0+nx]
    u   =   U[:,:,nx0:nx0+nx]
    v   =   V[:,:,nx0:nx0+nx]
    w   =   W[:,:,nx0:nx0+nx]
    p   =   P[:,:,nx0:nx0+nx]
    del Rho, U, V, W, P
    uf  = u - um
    vf  = v - vm
    wf  = w - wm
    pf  = p - pm
    uF  = u - umF
    vF  = v - vmF
    wF  = w - wmF
    TF  = p / (rho * Rgas) - TmF
    Tt  = p / (rho * Rgas) + 0.5e0 * (gamma - 1.e0) * (u**2 + v**2 + w**2) / gamma
    TtF = Tt - TtmF

    # x corr
    for k in range(Nz):
      for i in range(nx):
        uux[k,0,i] += uF[k,ny1,0] * uF[k,ny1,i]
        uux[k,1,i] += uF[k,ny2,0] * uF[k,ny2,i]
        vvx[k,0,i] += vF[k,ny1,0] * vF[k,ny1,i]
        vvx[k,1,i] += vF[k,ny2,0] * vF[k,ny2,i]
        wwx[k,0,i] += wF[k,ny1,0] * wF[k,ny1,i]
        wwx[k,1,i] += wF[k,ny2,0] * wF[k,ny2,i]

    # z corr
    for k in range(Nz):
      for i in range(nx):
        uuz[k,0,i] += uF[0,ny1,i] * uF[k,ny1,i]
        uuz[k,1,i] += uF[0,ny2,i] * uF[k,ny2,i]
        vvz[k,0,i] += vF[0,ny1,i] * vF[k,ny1,i]
        vvz[k,1,i] += vF[0,ny2,i] * vF[k,ny2,i]
        wwz[k,0,i] += wF[0,ny1,i] * wF[k,ny1,i]
        wwz[k,1,i] += wF[0,ny2,i] * wF[k,ny2,i]


    # Strong Reynolds analogy
    uu   += uF  * uF
    vv   += vF  * vF
    uv   += uF  * vF
    TT   += TF  * TF
    TtTt += TtF * TtF
    uT   += uF  * TF
    vT   += vF  * TF
    vTt  += vF  * TtF

    itr += 1.e0

  # corr
  UUx = np.mean(uux, axis=0)  / itr
  VVx = np.mean(vvx, axis=0)  / itr
  WWx = np.mean(wwx, axis=0)  / itr
  UUz = np.mean(uuz, axis=-1) / itr
  VVz = np.mean(vvz, axis=-1) / itr
  WWz = np.mean(wwz, axis=-1) / itr

  # Strong Reynolds analogy
  uu   /= itr
  vv   /= itr
  uv   /= itr
  TT   /= itr
  TtTt /= itr
  uT   /= itr
  vT   /= itr
  vTt  /= itr

  Ruv = uv / (np.sqrt(uu) * np.sqrt(vv))
  RvT = vT / (np.sqrt(vv) * np.sqrt(TT))
  RuT = uT / (np.sqrt(uu) * np.sqrt(TT))
  # theoretical relation
  RUV = -RvT * (1.e0 - vTt / vT)
  RUT = -1.e0 + 0.5e0 * TtTt / TT
  

  save_path = os.path.join(Q_directory, "corr_x.d")
  with open(save_path, "w", encoding="UTF-8") as f:
    print("# x         u1         u2          v1          v2          w1          w2", file=f)
    for i in range(nx):
      print(f'{x[i]:.3e}', f'{UUx[0,i]:.3e}', f'{UUx[1,i]:.3e}', f'{VVx[0,i]:.3e}', \
            f'{VVx[1,i]:.3e}', f'{WWx[0,i]:.3e}', f'{WWx[1,i]:.3e}', file=f)


  save_path = os.path.join(Q_directory, "corr_z.d")
  with open(save_path, "w", encoding="UTF-8") as f:
    print("# z         u1         u2          v1          v2          w1          w2", file=f)
    for k in range(Nz):
      print(f'{z[k]:.3e}', f'{UUz[k,0]:.3e}', f'{UUz[k,1]:.3e}', f'{VVz[k,0]:.3e}', \
            f'{VVz[k,1]:.3e}', f'{WWz[k,0]:.3e}', f'{WWz[k,1]:.3e}', file=f)


  save_path = os.path.join(Q_directory, "SRA.d")
  with open(save_path, "w", encoding="UTF-8") as f:
    print("# y         Ruv         RvT          RuT        RUV       RUT", file=f)
    for j in range(Ny):
      print(f'{y[j]/delta:.3e}', f'{Ruv[j]:.3e}', f'{RvT[j]:.3e}', f'{RuT[j]:.3e}', \
      f'{RUV[j]:.3e}', f'{RUT[j]:.3e}', file=f)


def main(Nx, Ny, Nz, nx0, nx, yp1, yp2, x, y, z):
  rho_path = os.path.join(Q_directory, "rho.npy")
  u_path   = os.path.join(Q_directory, "u.npy")
  v_path   = os.path.join(Q_directory, "v.npy")
  w_path   = os.path.join(Q_directory, "w.npy")
  p_path   = os.path.join(Q_directory, "p.npy")
  Rho = np.load(rho_path)
  U   = np.load(u_path)
  V   = np.load(v_path)
  W   = np.load(w_path)
  P   = np.load(p_path)


  uF_path  = os.path.join(Q_directory, "uF.npy")
  vF_path  = os.path.join(Q_directory, "vF.npy")
  wF_path  = os.path.join(Q_directory, "wF.npy")
  TF_path  = os.path.join(Q_directory, "TF.npy")
  TtF_path = os.path.join(Q_directory, "TtF.npy")
  UF  = np.load(uF_path)
  VF  = np.load(vF_path)
  WF  = np.load(wF_path)
  Tf  = np.load(wF_path)
  Ttf = np.load(wF_path)


  rho = Rho[:,:,nx0:nx0+nx]
  u   =   U[:,:,nx0:nx0+nx]
  v   =   V[:,:,nx0:nx0+nx]
  w   =   W[:,:,nx0:nx0+nx]
  p   =   P[:,:,nx0:nx0+nx]
  uF  =  UF[:,:,nx0:nx0+nx]
  vF  =  VF[:,:,nx0:nx0+nx]
  wF  =  WF[:,:,nx0:nx0+nx]
  TF  =  Tf[:,:,nx0:nx0+nx]
  TtF = Ttf[:,:,nx0:nx0+nx]

  del Rho, U, V, W, P, UF, VF, WF, Tf, Ttf

  calc_corr(Nx, Ny, Nz, nx0, nx, x, y, z, yp1, yp2, rho, u, v, w, p, uF, vF, wF, TF, TtF)

main(Nx, Ny, Nz, nx0, nx, yp1, yp2, x, y, z)

