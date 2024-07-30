import os
import numpy as np
import vtk
from vtk.util import numpy_support

def gridInfo(data_directory):
  file_path = os.path.join(data_directory, "x.npy")
  x         = np.load(file_path)
  file_path = os.path.join(data_directory, "y.npy")
  y         = np.load(file_path)
  file_path = os.path.join(data_directory, "z.npy")
  z         = np.load(file_path)
  Nx        = len(x)
  Ny        = len(y)
  Nz        = len(z)
  return np.float32(x), np.float32(y), np.float32(z), Nx, Ny, Nz

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
  return len(x), len(y), len(z), np.float32(x), np.float32(y), np.float32(z)

def getVector(file_path,Nx,Ny,Nz,name):
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(file_path)
  reader.Update()

  # get dataset
  Q = reader.GetOutput()
  V = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray(name))

  u = np.reshape(V[:,0], [Nz,Ny,Nx])
  v = np.reshape(V[:,1], [Nz,Ny,Nx])
  w = np.reshape(V[:,2], [Nz,Ny,Nx])
  return np.float32(u), np.float32(v), np.float32(w)

def getScalar(file_path,Nx,Ny,Nz,name):
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(file_path)
  reader.Update()

  # get dataset
  Q = reader.GetOutput()
  a = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray(name))
  a = np.reshape(a, [Nz,Ny,Nx])
  return np.float32(a)

def getMeanVector(directory_path,vtk_files,Nx,Ny,Nz,name):
  # decleare
  umean = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  vmean = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  wmean = np.zeros((Nz,Ny,Nx), dtype=np.float32)

  # calc mean velocity
  n = 1.e0
  for vtk_file in vtk_files:
    # get path
    file_path = os.path.join(directory_path, vtk_file)
    u, v, w = getVector(file_path,Nx,Ny,Nz,name)
    umean = ((n - 1.e0) * umean + u) / n
    vmean = ((n - 1.e0) * vmean + v) / n
    wmean = ((n - 1.e0) * wmean + w) / n
    n += 1.e0
  return umean, vmean, wmean

def getMeanScalar(directory_path,vtk_files,Nx,Ny,Nz,name):
  # decleare
  amean = np.zeros((Nz,Ny,Nx), dtype=np.float32)

  # calc mean velocity
  n = 1.e0
  for vtk_file in vtk_files:
    # get path
    file_path = os.path.join(directory_path, vtk_file)
    a = getScalar(file_path,Nx,Ny,Nz,name)
    amean = ((n - 1.e0) * amean + a) / n
    n += 1.e0
  return amean

