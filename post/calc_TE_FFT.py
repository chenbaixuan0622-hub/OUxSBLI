import numpy as np
import os
from tqdm import tqdm
from mod.mod_FFT import getV, reconstruct_u

########################################
Q_dir = "../../../../../../media/user/HD-EDS-E/hatayama/TBL/SBLI_05delta/stat03ms_09ms_SLAU/post"

u_path = os.path.join(Q_dir, "ua_upstream.npy")
v_path = os.path.join(Q_dir, "va_upstream.npy")
w_path = os.path.join(Q_dir, "wa_upstream.npy")

# non-dimensionalized length
Lx = 2.e0
Lz = 0.5e0

u0 = 506.67e0

########################################

# U[Nz,Ny,Nx,Nt], V[Nz,Ny,Nx,Nt], W[Nz,Ny,Nx,Nt]
U = np.load(u_path)
V = np.load(v_path)
W = np.load(w_path)

# non-dimensionalize
U /= u0
V /= u0
W /= u0

Nx = len(U[0,0,:,0])
Nz = len(U[:,0,0,0])

x = np.linspace(0.e0, Lx, Nx)
z = np.linspace(0.e0, Lz, Nz)

nx0 = 0
nx1 = int(Nx/4)
nz0 = 0
nz1 = int(Nz/4)

reconstruct_u(x[nx0:nx1], z[nz0:nz1], U[nz0:nz1,0,nx0:nx1,0], V[nz0:nz1,0,nx0:nx1,0], W[nz0:nz1,0,nx0:nx1,0])

'''
ke     = np.sum(0.5e0 * (U[nz0:nz1,0,nx0:nx1,0]**2 + V[nz0:nz1,0,nx0:nx1,0]**2 + W[nz0:nz1,0,nx0:nx1,0]**2))
ke_fft = reconstruct_ke(U[nz0:nz1,0,nx0:nx1,0], V[nz0:nz1,0,nx0:nx1,0], W[nz0:nz1,0,nx0:nx1,0])

print("total kinetic energy = ", ke)
print("reconstructed kinetic energy = ", ke_fft)
print(ke_fft / ke * 1e2, " % kinetic energy was reconstructed")
'''
