import numpy as np
import vtk
import os
from numba import njit


@njit(cache=True, fastmath=True, nogil=True)
def ReynoldsStress(Nx, Ny, Nz, Rho, uf, vf, wf, rhow, ut, ruu, rvv, rww, ruv):
  ruu += Rho * uf**2 / (rhow * ut**2)
  rvv += Rho * vf**2 / (rhow * ut**2)
  rww += Rho * wf**2 / (rhow * ut**2)
  ruv += Rho * uf * vf / (rhow * ut**2)

def print_ReynoldsStress(x, y, z, ruu, rvv, rww, ruv, directory, name):
  ruu1d = np.float32(ruu.flatten())
  rvv1d = np.float32(rvv.flatten())
  rww1d = np.float32(rww.flatten())
  ruv1d = np.float32(ruv.flatten())

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

  ruu = vtk.vtkFloatArray()
  ruu.SetName("ruu")
  for i in range(nx * ny * nz):
    ruu.InsertNextValue(ruu1d[i])
  grid.GetPointData().AddArray(ruu)

  rvv = vtk.vtkFloatArray()
  rvv.SetName("rvv")
  for i in range(nx * ny * nz):
    rvv.InsertNextValue(rvv1d[i])
  grid.GetPointData().AddArray(rvv)

  rww = vtk.vtkFloatArray()
  rww.SetName("rww")
  for i in range(nx * ny * nz):
    rww.InsertNextValue(rww1d[i])
  grid.GetPointData().AddArray(rww)

  ruv = vtk.vtkFloatArray()
  ruv.SetName("ruv")
  for i in range(nx * ny * nz):
    ruv.InsertNextValue(ruv1d[i])
  grid.GetPointData().AddArray(ruv)

  writer = vtk.vtkXMLRectilinearGridWriter()
  writer.SetFileName(filepath)
  writer.SetInputData(grid)
  writer.SetDataModeToAppended()
  writer.EncodeAppendedDataOff()
  writer.Write()

