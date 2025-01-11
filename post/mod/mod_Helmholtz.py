import numpy as np
import os
import vtk
import matplotlib.pyplot as plt
from mod.mod_nabla import Scalar, Vector, bc


def plot_Helmholtz_decomposition(U, phi, x, y, z=None):
  X, Y = np.meshgrid(x*1e3, y*1e3)

  if z is None:
    gradPhi    = Scalar(phi, x, y).gradient()
    rotA       = U + gradPhi
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
  rotA    = U - gradPhi

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

