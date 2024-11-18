import numpy as np
import os
import re
from tqdm import tqdm
from numba import jit
from statsmodels.tsa.stattools import acf
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import non_dim_tbl, Sutherland, tau
from mod.mod_SRA import SRA, corr_x, corr_z

#Q_directory = "../../SBLI_1delta/stat03ms_05ms_SLAU"
Q_directory = "../../../../../../media/user/HD-EDS-E/hatayama/TBL/SBLI_05delta/stat03ms_09ms_SLAU"

Q_files   = [f for f in os.listdir(Q_directory) if f.endswith(".vtr")]

def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')

Q_files.sort(key=extract_number)

gamma = 1.4e0
Rgas  = 287.03

# target yplus
yp1 = 5.e0
yp2 = 30.e0

endT = 0.6e-3

first_path          = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, X, y, z = getGrid(first_path)
delta = 25.e0

nx0 = int(6.e0 / delta * Nx)
nx2 = int(8.e0 / delta * Nx)

nx  = nx2 - nx0
x = X[nx0:nx0+nx]

def calc_corr(Nx, Ny, Nz, nx0, nx, x, y, z, yp1, yp2, endT, rhom, um, vm, wm, pm, umF, vmF, wmF, TmF, TtmF):
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
  uu   = np.zeros((Nz,Ny-1,nx), dtype=np.float32)
  vv   = np.zeros((Nz,Ny-1,nx), dtype=np.float32)
  uv   = np.zeros((Nz,Ny-1,nx), dtype=np.float32)
  TT   = np.zeros((Nz,Ny-1,nx), dtype=np.float32)
  TtTt = np.zeros((Nz,Ny-1,nx), dtype=np.float32)
  uT   = np.zeros((Nz,Ny-1,nx), dtype=np.float32)
  vT   = np.zeros((Nz,Ny-1,nx), dtype=np.float32)
  vTt  = np.zeros((Nz,Ny-1,nx), dtype=np.float32)

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

  # calc auto correlation
  t      = np.linspace(0.e0, endT, len(Q_files), dtype=np.float32)
  ua     = np.zeros((Nz,2,nx,len(Q_files)), dtype=np.float32)
  va     = np.zeros((Nz,2,nx,len(Q_files)), dtype=np.float32)
  wa     = np.zeros((Nz,2,nx,len(Q_files)), dtype=np.float32)

  itr = 0.e0 
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_directory, Q_file)
    if file_path == os.path.join(Q_directory, "TKE.vtr") or \
       file_path == os.path.join(Q_directory, "ReynoldsStress.vtr") or \
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

    corr_x(nx, ny1, ny2, Nz, uf, vf, wf, uux, vvx, wwx)
    corr_z(nx, ny1, ny2, Nz, uf, vf, wf, uuz, vvz, wwz)
    SRA(nx, Ny, Nz, uF, vF, TF, TtF, uu, vv, uv, TT, TtTt, uT, vT, vTt)
    
    # auto correlation
    ua[:,0,:,int(itr)] = u[:,ny1,:]
    va[:,0,:,int(itr)] = v[:,ny1,:]
    wa[:,0,:,int(itr)] = w[:,ny1,:]
    ua[:,1,:,int(itr)] = u[:,ny2,:]
    va[:,1,:,int(itr)] = v[:,ny2,:]
    wa[:,1,:,int(itr)] = w[:,ny2,:]
  
    itr += 1.e0

  # corr
  UUx = np.mean(uux, axis=0)  / itr
  VVx = np.mean(vvx, axis=0)  / itr
  WWx = np.mean(wwx, axis=0)  / itr
  UUz = np.mean(uuz, axis=-1) / itr
  VVz = np.mean(vvz, axis=-1) / itr
  WWz = np.mean(wwz, axis=-1) / itr

  # auto corr
  save_path = os.path.join(Q_directory, "ua_upstream.npy")
  np.save(save_path, ua[:,:,:,:int(itr)])
  save_path = os.path.join(Q_directory, "va_upstream.npy")
  np.save(save_path, va[:,:,:,:int(itr)])
  save_path = os.path.join(Q_directory, "wa_upstream.npy")
  np.save(save_path, wa[:,:,:,:int(itr)])

  u_acf1 = np.zeros(int(itr), dtype=np.float32)
  v_acf1 = np.zeros(int(itr), dtype=np.float32)
  w_acf1 = np.zeros(int(itr), dtype=np.float32)
  u_acf2 = np.zeros(int(itr), dtype=np.float32)
  v_acf2 = np.zeros(int(itr), dtype=np.float32)
  w_acf2 = np.zeros(int(itr), dtype=np.float32)
  for k in range(Nz):
    for i in range(nx):
      u_acf1 += acf(ua[k,0,i,:int(itr)], nlags=int(itr))
      v_acf1 += acf(va[k,0,i,:int(itr)], nlags=int(itr))
      w_acf1 += acf(wa[k,0,i,:int(itr)], nlags=int(itr))
      u_acf2 += acf(ua[k,1,i,:int(itr)], nlags=int(itr))
      v_acf2 += acf(va[k,1,i,:int(itr)], nlags=int(itr))
      w_acf2 += acf(wa[k,1,i,:int(itr)], nlags=int(itr))

  u_acf1 /= float(nx*Nz)
  v_acf1 /= float(nx*Nz)
  w_acf1 /= float(nx*Nz)
  u_acf2 /= float(nx*Nz)
  v_acf2 /= float(nx*Nz)
  w_acf2 /= float(nx*Nz)

  # Strong Reynolds analogy
  uu   /= itr
  vv   /= itr
  uv   /= itr
  TT   /= itr
  TtTt /= itr
  uT   /= itr
  vT   /= itr
  vTt  /= itr

  Ruv = np.mean(uv / (np.sqrt(uu) * np.sqrt(vv)), axis=(0,-1))
  RvT = np.mean(vT / (np.sqrt(vv) * np.sqrt(TT)), axis=(0,-1))
  RuT = np.mean(uT / (np.sqrt(uu) * np.sqrt(TT)), axis=(0,-1))
  # theoretical relation
  RUV = -RvT * np.mean((1.e0 - vTt) / vT, axis=(0,-1))
  RUT = -1.e0 + 0.5e0 * np.mean(TtTt / TT, axis=(0,-1))
  
  save_path = os.path.join(Q_directory, "corr_x_upstream.d")
  with open(save_path, "w", encoding="UTF-8") as f:
    print("# x       u1        u2         v1         v2         w1         w2", file=f)
    for i in range(nx):
      print(f'{x[i]:.3e}', f'{UUx[0,i]/UUx[0,0]:.3e}', f'{UUx[1,i]/UUx[1,0]:.3e}', f'{VVx[0,i]/VVx[0,0]:.3e}', \
            f'{VVx[1,i]/VVx[1,0]:.3e}', f'{WWx[0,i]/WWx[0,0]:.3e}', f'{WWx[1,i]/WWx[1,0]:.3e}', file=f)


  save_path = os.path.join(Q_directory, "corr_z_upstream.d")
  with open(save_path, "w", encoding="UTF-8") as f:
    print("# z       u1        u2         v1         v2         w1         w2", file=f)
    for k in range(Nz):
      print(f'{z[k]:.3e}', f'{UUz[k,0]/UUz[0,0]:.3e}', f'{UUz[k,1]/UUz[0,1]:.3e}', f'{VVz[k,0]/VVz[0,0]:.3e}', \
            f'{VVz[k,1]/VVz[0,1]:.3e}', f'{WWz[k,0]/WWz[0,0]:.3e}', f'{WWz[k,1]/WWz[0,1]:.3e}', file=f)


  save_path = os.path.join(Q_directory, "SRA_upstream.d")
  with open(save_path, "w", encoding="UTF-8") as f:
    print("# y       Ruv         RvT          RuT        RUV       RUT", file=f)
    for j in range(Ny-1):
      print(f'{y[j]/delta:.3e}', f'{Ruv[j]:.3e}', f'{RvT[j]:.3e}', f'{RuT[j]:.3e}', \
      f'{RUV[j]:.3e}', f'{RUT[j]:.3e}', file=f)


  save_path = os.path.join(Q_directory, "auto_corr_upstream.d")
  with open(save_path, "w", encoding="UTF-8") as f:
    print("# t       u_acf       u_acf       v_acf        v_acf       w_acf       w_acf", file=f)
    for i in range(int(itr)):
      print(f'{t[i]:.3e}', f'{u_acf1[i]:.3e}', f'{u_acf2[i]:.3e}', f'{v_acf1[i]:.3e}', f'{v_acf2[i]:.3e}', \
            f'{w_acf1[i]:.3e}', f'{w_acf2[i]:.3e}', file=f)


def main(Nx, Ny, Nz, nx0, nx, yp1, yp2, endT, x, y, z):
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

  calc_corr(Nx, Ny, Nz, nx0, nx, x, y, z, yp1, yp2, endT, rho, u, v, w, p, uF, vF, wF, TF, TtF)

main(Nx, Ny, Nz, nx0, nx, yp1, yp2, endT, x, y, z)

