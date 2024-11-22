import numpy as np
import vtk
import os
from numba import njit
from mod.mod_turb_stat import Sutherland


@njit(cache=True, fastmath=True, nogil=True)
def tau_wall(nx, nz, y, rho, u, p):
  Rgas = 287.03e0
  tw   = np.zeros((nz,nx), dtype=np.float32)
  for k in range(nz):
    for i in range(nx):
      mu      = Sutherland(p[k,0,i] / (Rgas * rho[k,0,i]))
      tw[k,i] = mu * (-u[k,0,i] + u[k,1,i]) / (-y[0] + y[1])
  return tw


def print_tau(x, y, z, tau, directory, name):
  tau1d = np.float32(tau.flatten())

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
  nz = len(z)
  
  for i in range(nx):
    x_coords.InsertNextValue(x[i])
  for j in range(1):
    y_coords.InsertNextValue(0.e0)
  for k in range(nz):
    z_coords.InsertNextValue(z[k])
  
  grid = vtk.vtkRectilinearGrid()
  grid.SetDimensions(nx, 1, nz)
  grid.SetXCoordinates(x_coords)
  grid.SetYCoordinates(y_coords)
  grid.SetZCoordinates(z_coords)

  tau = vtk.vtkFloatArray()
  tau.SetName("tau")
  for i in range(nx * nz):
    tau.InsertNextValue(tau1d[i])
  grid.GetPointData().AddArray(tau)

  writer = vtk.vtkXMLRectilinearGridWriter()
  writer.SetFileName(filepath)
  writer.SetInputData(grid)
  writer.Write()

