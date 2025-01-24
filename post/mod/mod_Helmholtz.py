import numpy as np
import jax
import jax.numpy as jnp
import jax.scipy.linalg as jsl
from functools import partial
from scipy.linalg import lu_factor, lu_solve
import os
import vtk
import matplotlib.pyplot as plt
from mod.mod_nabla import Scalar, Vector, bc
#from mod.mod_matrix import Poisson_coef, create_matrix, wall_Neumann_bc_A, wall_Neumann_bc_B
from mod.mod_matrix_jax import Poisson_coef, create_matrix, wall_Neumann_bc_A, wall_Neumann_bc_B


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


def plot_Helmholtz_decomposition(U, phi, x, y, z=None):
  X, Y = np.meshgrid(x*1e3, y*1e3)

  if z is None:
    gradPhi    = Scalar(phi, x, y).gradient()
    rotA       = U - gradPhi
    rotgradPhi = Vector(gradPhi, x, y).rotation()
    divgradPhi = Vector(gradPhi, x, y).divergence()
    divrotA    = Vector(rotA, x, y).divergence()
    rotrotA    = Vector(rotA, x, y).rotation()

    plt.contourf(X, Y, phi, levels=100, cmap='jet')
    plt.colorbar()
    plt.show()
    plt.close()

    fig, ax  = plt.subplots(2, 3, figsize=(12, 6))
    contour1 = ax[0,0].contourf(X, Y, gradPhi[0,:,:], levels=100, cmap='jet')
    contour2 = ax[0,1].contourf(X, Y, rotgradPhi,     levels=100, cmap='jet')
    contour3 = ax[0,2].contourf(X, Y, divgradPhi,     levels=100, cmap='jet')
    contour4 = ax[1,0].contourf(X, Y, rotA[0,:,:],    levels=100, cmap='jet')
    contour5 = ax[1,1].contourf(X, Y, divrotA,        levels=100, cmap='jet')
    contour6 = ax[1,2].contourf(X, Y, rotrotA,        levels=100, cmap='jet')
  else:
    nz         = len(z) // 2
    gradPhi    = Scalar(phi, x, y, z).gradient()
    rotA       = U - gradPhi
    rotgradPhi = Vector(gradPhi, x, y, z).rotation()
    divgradPhi = Vector(gradPhi, x, y, z).divergence()
    divrotA    = Vector(rotA, x, y, z).divergence()
    rotrotA    = Vector(rotA, x, y, z).rotation()

    plt.contourf(X, Y, phi[nz,:,:], levels=100, cmap='jet')
    plt.colorbar()
    plt.show()
    plt.close()

    fig, ax  = plt.subplots(2, 3, figsize=(12, 6))
    contour1 = ax[0,0].contourf(X, Y, gradPhi[0,nz,:,:],    levels=100, cmap='jet')
    contour2 = ax[0,1].contourf(X, Y, rotgradPhi[0,nz,:,:], levels=100, cmap='jet')
    contour3 = ax[0,2].contourf(X, Y, divgradPhi[nz,:,:],   levels=100, cmap='jet')
    contour4 = ax[1,0].contourf(X, Y, rotA[0,nz,:,:],       levels=100, cmap='jet')
    contour5 = ax[1,1].contourf(X, Y, divrotA[nz,:,:],      levels=100, cmap='jet')
    contour6 = ax[1,2].contourf(X, Y, rotrotA[0,nz,:,:],    levels=100, cmap='jet')

  cbar1    = fig.colorbar(contour1, ax=ax[0,0], \
                          orientation='horizontal', pad=0.1, fraction=0.046, location='top')
  cbar2    = fig.colorbar(contour2, ax=ax[0,1], \
                          orientation='horizontal', pad=0.1, fraction=0.046, location='top')
  cbar3    = fig.colorbar(contour3, ax=ax[0,2], \
                          orientation='horizontal', pad=0.1, fraction=0.046, location='top')
  cbar4    = fig.colorbar(contour4, ax=ax[1,0], \
                          orientation='horizontal', pad=0.1, fraction=0.046, location='top')
  cbar5    = fig.colorbar(contour5, ax=ax[1,1], \
                          orientation='horizontal', pad=0.1, fraction=0.046, location='top')
  cbar6    = fig.colorbar(contour6, ax=ax[1,2], \
                          orientation='horizontal', pad=0.1, fraction=0.046, location='top')
  plt.show()


def save_Helmholtz_decomposition(U, phi, X, Y, Z, dir, name):
  gradPhi = Scalar(phi, X, Y, Z).gradient()
  for i in range(3):
    gradPhi[i,:,:,:] = bc(gradPhi[i,:,:,:]).periodic(z=True)
    gradPhi[i,:,:,:] = bc(gradPhi[i,:,:,:]).Neumann(x1=True, x2=True, y2=True)
    gradPhi[i,:,:,:] = bc(gradPhi[i,:,:,:]).Dirichlet(y1=0.e0)
  rotA  = U - gradPhi

  phi1d = np.float32(phi[1:-1,1:-1,1:-1].flatten())
  
  '''
  rotA = Vector(A, X, Y, Z).rotation()
  gradPhi = U - rotA
  '''

  gradPhix = np.float32(gradPhi[0,1:-1,1:-1,1:-1].flatten())
  gradPhiy = np.float32(gradPhi[1,1:-1,1:-1,1:-1].flatten())
  gradPhiz = np.float32(gradPhi[2,1:-1,1:-1,1:-1].flatten())
  rotAx    = np.float32(rotA[0,1:-1,1:-1,1:-1].flatten())
  rotAy    = np.float32(rotA[1,1:-1,1:-1,1:-1].flatten())
  rotAz    = np.float32(rotA[2,1:-1,1:-1,1:-1].flatten())

  x = X[1:-1]
  y = Y[1:-1]
  z = Z[1:-1]

  os.makedirs(dir, exist_ok=True)
  filename  = name + ".vtr" 
  filepath  = os.path.join(dir, filename)

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

  grid = vtk.vtkRectilinearGrid()
  grid.SetDimensions(nx, ny, nz)
  grid.SetXCoordinates(x_coords)
  grid.SetYCoordinates(y_coords)
  grid.SetZCoordinates(z_coords)

  phi = vtk.vtkFloatArray()
  phi.SetName("phi")
  for i in range(nx*ny*nz):
    phi.InsertNextValue(phi1d[i])
  grid.GetPointData().AddArray(phi)

  grad = vtk.vtkFloatArray()
  grad.SetName("gradPhi")
  grad.SetNumberOfComponents(3)
  for i in range(nx*ny*nz):
    grad.InsertNextTuple3(gradPhix[i], gradPhiy[i], gradPhiz[i])
  grid.GetPointData().AddArray(grad)

  rot = vtk.vtkFloatArray()
  rot.SetName("rotA")
  rot.SetNumberOfComponents(3)
  for i in range(nx*ny*nz):
    rot.InsertNextTuple3(rotAx[i], rotAy[i], rotAz[i])
  grid.GetPointData().AddArray(rot)

  writer = vtk.vtkXMLRectilinearGridWriter()
  writer.SetFileName(filepath)
  writer.SetInputData(grid)
  writer.Write()

