import os
import numpy as np
import vtk
from vtk.util import numpy_support

def getGrid(file_path):
  # make VTK Structured Grid Reader
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(file_path)
  reader.Update()
  
  # get grid
  grid = reader.GetOutput()
  x = numpy_support.vtk_to_numpy(grid.GetXCoordinates())
  y = numpy_support.vtk_to_numpy(grid.GetYCoordinates())
  z = numpy_support.vtk_to_numpy(grid.GetZCoordinates())
  return len(x), len(y), len(z), x, y, z

def getVelocity(file_path,Nx,Ny,Nz):
  # make VTK Structured Grid Reader
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(file_path)
  reader.Update()

  # get dataset
  Q = reader.GetOutput()
  V = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray("velocity"))

  u = np.reshape(V[:,0], [Nz,Ny,Nx])
  v = np.reshape(V[:,1], [Nz,Ny,Nx])
  w = np.reshape(V[:,2], [Nz,Ny,Nx])
  return u, v, w

def getMean(directory_path,vtk_files,Nx,Ny,Nz):
  # decleare
  umean = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  vmean = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  wmean = np.zeros((Nz,Ny,Nx), dtype=np.float32)

  # calc mean velocity
  n = 1.e0
  for vtk_file in vtk_files:
    # get path
    file_path = os.path.join(directory_path, vtk_file)
    u, v, w = getVelocity(file_path,Nx,Ny,Nz)
    umean = ((n - 1.e0) * umean + u) / n
    vmean = ((n - 1.e0) * vmean + v) / n
    wmean = ((n - 1.e0) * wmean + w) / n
    n += 1.e0
  return umean, vmean, wmean

