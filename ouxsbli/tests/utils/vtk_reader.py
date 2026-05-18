import numpy as np
import os
import re
import vtk
from vtk.util import numpy_support


def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')


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
  return len(x), len(y), len(z), \
         x.astype(np.float32, copy=False), \
         y.astype(np.float32, copy=False), \
         z.astype(np.float32, copy=False)


def getVector(file_path, Nx, Ny, Nz, name):
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(file_path)
  reader.GetPointDataArraySelection().DisableAllArrays()
  reader.GetPointDataArraySelection().EnableArray(name)
  reader.Update()
  # get dataset
  Q = reader.GetOutput()
  V = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray(name))
  V = V.reshape((Nz,Ny,Nx,3))
  u = V[:,:,:,0]
  v = V[:,:,:,1]
  w = V[:,:,:,2]
  return u.astype(np.float32, copy=False), \
         v.astype(np.float32, copy=False), \
         w.astype(np.float32, copy=False)


def getScalar(file_path, Nx, Ny, Nz, name):
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(file_path)
  reader.GetPointDataArraySelection().DisableAllArrays()
  reader.GetPointDataArraySelection().EnableArray(name)
  reader.Update()
  # get dataset
  Q = reader.GetOutput()
  a = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray(name))
  a = a.reshape((Nz,Ny,Nx))
  return a.astype(np.float32, copy=False)


def getQ(file_path, Nx, Ny, Nz, reader=None):
  if reader is None:
    reader = vtk.vtkXMLRectilinearGridReader()
    reader.GetPointDataArraySelection().DisableAllArrays()
    reader.GetPointDataArraySelection().EnableArray("rho")
    reader.GetPointDataArraySelection().EnableArray("velocity")
    reader.GetPointDataArraySelection().EnableArray("p")
  reader.SetFileName(file_path)
  reader.Update()
  # get dataset
  Q   = reader.GetOutput()
  rho = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray("rho"))
  rho = rho.reshape((Nz,Ny,Nx))
  V   = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray("velocity"))
  V   = V.reshape((Nz,Ny,Nx,3))
  u   = V[:,:,:,0]
  v   = V[:,:,:,1]
  w   = V[:,:,:,2]
  p   = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray("p"))
  p   = p.reshape((Nz,Ny,Nx))
  return rho.astype(np.float32, copy=False), \
         u.astype(np.float32, copy=False), \
         v.astype(np.float32, copy=False), \
         w.astype(np.float32, copy=False), \
         p.astype(np.float32, copy=False)


def initial_vtr(data_dir):
  import glob
  files = glob.glob(os.path.join(str(data_dir), "Q*.vtr"))
  if not files:
    raise FileNotFoundError(f"No Q*.vtr files found in {data_dir}")
  return min(files, key=lambda f: extract_number(os.path.basename(f)))


def latest_vtr(data_dir):
  import glob
  files = glob.glob(os.path.join(str(data_dir), "Q*.vtr"))
  if not files:
    raise FileNotFoundError(f"No Q*.vtr files found in {data_dir}")
  return max(files, key=lambda f: extract_number(os.path.basename(f)))

