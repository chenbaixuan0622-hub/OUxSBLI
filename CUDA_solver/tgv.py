import numpy as np
import cupy as cp


def set_grid(nx, ny, nz, Lx, Ly, Lz):
  dx = Lx / np.float64(nx-2) 
  dy = Ly / np.float64(ny-2) 
  dz = Lz / np.float64(nz-2) 
  x  = np.zeros(nx+1, dtype=np.float64)
  y  = np.zeros(ny+1, dtype=np.float64)
  z  = np.zeros(nz+1, dtype=np.float64)
  xc = np.zeros(nx, dtype=np.float64)
  yc = np.zeros(ny, dtype=np.float64)
  zc = np.zeros(nz, dtype=np.float64)
  for i in range(1, nx):
    x[i] = dx * np.float64(i-1)
  x[0]  = x[1]  - dx
  x[-1] = x[-2] + dx
  y  = x
  z  = y
  xc = 0.5e0 * (x[:-1] + x[1:])
  yc = 0.5e0 * (y[:-1] + y[1:])
  zc = 0.5e0 * (z[:-1] + z[1:])
  return xc, yc, zc, dx, dy, dz


def set_init(nx, ny, nz, x, y, z, gamma, Rgas, M0, rho0):
  Q = np.zeros((nz, ny, nx, 5), dtype=np.float64)
  for k in range(nz):
    for j in range(ny):
      for i in range(nx):
        Q[k,j,i,0] =  rho0
        Q[k,j,i,1] =  rho0 * M0 * np.sin(x[i]) * np.cos(y[j]) * np.cos(z[k])
        Q[k,j,i,2] = -rho0 * M0 * np.cos(x[i]) * np.sin(y[j]) * np.cos(z[k])
        Q[k,j,i,3] = 0.e0
        Q[k,j,i,4] = (1.e0 / gamma + 0.0625e0 * rho0 * (M0**2) * \
        (np.cos(2.e0 * x[i]) + np.cos(2.e0 * y[j])) * (np.cos(2.e0 * z[k]) + 2.e0) \
        ) / (gamma - 1.e0) + 0.5e0 * (Q[k,j,i,1]**2 + Q[k,j,i,2]**2) / Q[k,j,i,0]
  # x direction
  Q[1:-1,1:-1,0,:]  = Q[1:-1,1:-1,-2,:]
  Q[1:-1,1:-1,-1,:] = Q[1:-1,1:-1,1,:]
  # y direction
  Q[1:-1,0,1:-1,:]  = Q[1:-1,-2,1:-1,:]
  Q[1:-1,-1,1:-1,:] = Q[1:-1,1,1:-1,:]
  # corner
  Q[1:-1,0,0,:]   = Q[1:-1,-2,-2,:]
  Q[1:-1,0,-1,:]  = Q[1:-1,-2,1,:]
  Q[1:-1,-1,0,:]  = Q[1:-1,1,-2,:]
  Q[1:-1,-1,-1,:] = Q[1:-1,1,1,:]
  # z direction
  Q[0,:,:,:]  = Q[-2,:,:,:]
  Q[-1,:,:,:] = Q[1,:,:,:]
  return Q


def set_bc(Q, nx, ny, nz):
  Q = cp.reshape(Q, (nz,ny,nx,5))
  # x direction
  Q[1:-1,1:-1,0,:]  = Q[1:-1,1:-1,-2,:]
  Q[1:-1,1:-1,-1,:] = Q[1:-1,1:-1,1,:]
  # y direction
  Q[1:-1,0,1:-1,:]  = Q[1:-1,-2,1:-1,:]
  Q[1:-1,-1,1:-1,:] = Q[1:-1,1,1:-1,:]
  # corner
  Q[1:-1,0,0,:]   = Q[1:-1,-2,-2,:]
  Q[1:-1,0,-1,:]  = Q[1:-1,-2,1,:]
  Q[1:-1,-1,0,:]  = Q[1:-1,1,-2,:]
  Q[1:-1,-1,-1,:] = Q[1:-1,1,1,:]
  # z direction
  Q[0,:,:,:]  = Q[-2,:,:,:]
  Q[-1,:,:,:] = Q[1,:,:,:]
  return cp.ravel(Q)

