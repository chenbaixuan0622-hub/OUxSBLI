import numpy as np
from scipy.sparse import diags
import jax
import jax.numpy as jnp
import jax.scipy.linalg as jsl
from functools import partial
from scipy.linalg import lu_factor, lu_solve
import os
import matplotlib.pyplot as plt
from tqdm import tqdm
from mod.mod_plot import set_Params
from mod.mod_read import getGrid, getVector, extract_number, Vector_Data
from mod.mod_POD import make_data, make_grid
from mod.mod_nabla import Vector, Scalar
#from mod.mod_matrix import Poisson_coef, create_matrix, wall_Neumann_bc_A, wall_Neumann_bc_B
from mod.mod_matrix_jax import Poisson_coef, create_matrix, wall_Neumann_bc_A, wall_Neumann_bc_B
from mod.mod_Helmholtz import plot_Helmholtz_decomposition, save_Helmholtz_decomposition


set_Params()

# parameter
Q_dir   = "../3D_solver/TBL/data"
Lx1     = 24.e-3#4.e-3
Lx2     = 40.e-3#20.e-3
Ly1     = 0.e-3
Ly2     = 10.e-3
stridex = 4
stridez = 4
ny      = 100
endT    = 0.1e-3


def HD(x, y, div, x0, x1, y1, mu=0.e0):
  nx, ny = len(x), len(y)

  a, A = Poisson_coef(x, y, mu)
  A    = wall_Neumann_bc_A(A, nx, ny)
  factor, piv = lu_factor(A)
  f = np.zeros_like(div)
  f[1:-1,1:-1] = div[1:-1,1:-1] / a

  B   = wall_Neumann_bc_B(f, nx, ny, x0, x1, y1)
  phi = lu_solve((factor, piv), B)
  return phi.reshape((ny,nx))


@partial(jax.jit, static_argnums=(0, 1))
def HD_jax(nx, ny, x, y, div, x0, x1, y1, mu=None):
  if mu is None:
    a, A = Poisson_coef(x, y)
    A    = wall_Neumann_bc_A(A, nx, ny)
    factor, piv = jsl.lu_factor(A)
    f = jnp.zeros_like(div)
    f = f.at[1:-1,1:-1].set(div[1:-1,1:-1] / a)

    B   = wall_Neumann_bc_B(f, nx, ny, x0, x1, y1)
  
    phi = jsl.lu_solve((factor, piv), B)
    return phi.reshape((ny,nx))
  else:
    a, A = Poisson_coef(x, y, mu)
    A    = wall_Neumann_bc_A(A, nx, ny)
    factor, piv = jsl.lu_factor(A)
    f = jnp.zeros_like(div)
    f = f.at[1:-1,1:-1].set(div[1:-1,1:-1] / a)

    B   = wall_Neumann_bc_B(f, nx, ny, x0, x1, y1)
  
    factor = jax.lax.complex(factor, jnp.zeros_like(factor))
    phi = jsl.lu_solve((factor, piv), B)
    return phi.reshape((ny,nx))


def HD3D(x, y, z, div, x0, x1, y1):
  nx, ny, nz = len(x), len(y), len(z)
  dz  = -z[0] + z[1]
  div = np.fft.fft(div, axis=0)
  x0  = np.fft.fft(x0,  axis=0)
  x1  = np.fft.fft(x1,  axis=0)
  y1  = np.fft.fft(y1,  axis=0)
  p   = np.zeros_like(div)

  for k in tqdm(range(nz)):
    mu = -2.e0 * (1.e0 - np.cos(2.e0 * np.pi * k / nz)) / dz**2
    p[k,:,:] = HD(x, y, div[k,:,:], x0, x1, y1, mu)
  return np.real(np.fft.ifft(p, axis=0))


@partial(jax.jit, static_argnums=(0, 1, 2))
def HD3D_jax(nx, ny, nz, x, y, z, div, x0, x1, y1):
  dz   = -z[0] + z[1]
  div  = jnp.fft.fft(div, axis=0)
  x0   = jnp.fft.fft(x0,  axis=0)
  x1   = jnp.fft.fft(x1,  axis=0)
  y1   = jnp.fft.fft(y1,  axis=0)
  p    = jnp.zeros_like(div)
  div0 = jnp.real(div[0,:,:])
  x00  = jnp.real(x0)
  x10  = jnp.real(x1)
  y10  = jnp.real(y1)
  p    = p.at[0,:,:].set(HD_jax(nx, ny, x, y, div0, x00, x10, y10))

  def body(k, p):
    mu = -2.e0 * (1.e0 - jnp.cos(2.e0 * jnp.pi * k / nz)) / dz**2
    p  = p.at[k,:,:].set(HD_jax(nx, ny, x, y, div[k,:,:], x0, x1, y1, mu))
    return p

  p = jax.lax.fori_loop(1, nz, body, p)
  return jnp.real(jnp.fft.ifft(p, axis=0))


def main(Q_dir, Lx1, Lx2, Ly1, Ly2, stridex, stridez, ny, endT):
  Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)
  Data = Vector_Data(Q_dir)
  umean, vmean, wmean, _, _, _ = Data.getMeanVector_interp('velocity', Lx1, Lx2, Ly1, Ly2, stridex, stridez, ny)
  
  itr  = 0
  for Q_file in tqdm(Q_files):
    file_path = os.path.join(Q_dir, Q_file)
    u, v, w, x, y, z = Data.getVector_interp(file_path, 'velocity', Lx1, Lx2, Ly1, Ly2, stridex, stridez, ny)
    U = np.array([u, v, w])
  
    x0 = np.zeros(len(y), dtype=np.float32)
    x1 = np.zeros(len(y), dtype=np.float32)
    y1 = np.zeros(len(x), dtype=np.float32)
  
    '''
    # 2D
    U   = U[0:2,5,:,:]

    div = Vector(U, x, y).divergence()

    #phi = HD(x, y, div, x0, x1, y1)
    div = jnp.array(div)
    phi = HD_jax(len(x), len(y), x, y, div, x0, x1, y1)
    phi = np.array(phi)
    plot_Helmholtz_decomposition(U, phi, x, y)
    '''
    # 3D
    div = Vector(U, x, y, z).divergence()
    #phi = HD3D(x, y, z, div, x0, x1, y1)

    div = jnp.array(div)
    phi = HD3D_jax(len(x), len(y), len(z), x, y, z, div, x0, x1, y1)
    phi = np.array(phi)
  
    name = "HD" + str(itr).zfill(5)
    HD_dir = os.path.join(Q_dir, "HD")
    os.makedirs(HD_dir, exist_ok=True)
    save_Helmholtz_decomposition(U, phi, x, y, z, HD_dir, name)
    itr += 1

main(Q_dir, Lx1, Lx2, Ly1, Ly2, stridex, stridez, ny, endT)

