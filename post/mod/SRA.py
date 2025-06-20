import numpy as np
import vtk
import os
from numba import njit


@njit(cache=True, fastmath=True, nogil=True)
def SRA(uf, vf, Tf, Ttf, uT, vT, uv, vT1, vTt, T2, Tt2):
  uT  += uf * Tf / np.sqrt(uf**2 * Tf**2)
  vT  += vf * Tf / np.sqrt(vf**2 * Tf**2)
  uv  += uf * vf / np.sqrt(uf**2 * vf**2)
  vT1 += vf * Tf
  vTt += vf * Ttf
  T2  += Tf**2
  Tt2 += Ttf**2


def print_SRA(x, y, z, uT, vT, uv, directory, name):
  #uveq   = -vT * factor1
  #uTeq   = -1.e0 + factor2
  uT1d   = np.float32(uT.flatten())
  vT1d   = np.float32(vT.flatten())
  uv1d   = np.float32(uv.flatten())
  #uveq1d = np.float32(uveq.flatten())
  #uTeq1d = np.float32(uTeq.flatten())


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

  uT = vtk.vtkFloatArray()
  uT.SetName("uT")
  for i in range(nx * ny * nz):
    uT.InsertNextValue(uT1d[i])
  grid.GetPointData().AddArray(uT)

  vT = vtk.vtkFloatArray()
  vT.SetName("vT")
  for i in range(nx * ny * nz):
    vT.InsertNextValue(vT1d[i])
  grid.GetPointData().AddArray(vT)

  uv = vtk.vtkFloatArray()
  uv.SetName("uv")
  for i in range(nx * ny * nz):
    uv.InsertNextValue(uv1d[i])
  grid.GetPointData().AddArray(uv)

  '''
  uveq = vtk.vtkFloatArray()
  uveq.SetName("uv_eq")
  for i in range(nx * ny * nz):
    uv.InsertNextValue(uveq1d[i])
  grid.GetPointData().AddArray(uveq)

  uTeq = vtk.vtkFloatArray()
  uTeq.SetName("uT_eq")
  for i in range(nx * ny * nz):
    uv.InsertNextValue(uTeq1d[i])
  grid.GetPointData().AddArray(uTeq)
  '''

  writer = vtk.vtkXMLRectilinearGridWriter()
  writer.SetFileName(filepath)
  writer.SetInputData(grid)
  writer.Write()

