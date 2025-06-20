import numpy as np
import os
import re
import vtk
from vtk.util import numpy_support
from tqdm import tqdm


def extract_number(filename):
  match = re.search(r'Q(\d+)\.vtr$', filename)
  if match:
    return int(match.group(1))
  return float('inf')


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
  reader.GetPointDataArraySelection().DisableAllArrays()
  reader.GetPointDataArraySelection().EnableArray(name)
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
  reader.GetPointDataArraySelection().DisableAllArrays()
  reader.GetPointDataArraySelection().EnableArray(name)
  reader.Update()
  # get dataset
  Q = reader.GetOutput()
  a = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray(name))
  a = np.reshape(a, [Nz,Ny,Nx])
  return np.float32(a)


def getMeanVector(directory_path,vtk_files,Nx,Ny,Nz,name):
  um = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  vm = np.zeros((Nz,Ny,Nx), dtype=np.float32)
  wm = np.zeros((Nz,Ny,Nx), dtype=np.float32)

  for vtk_file in tqdm(vtk_files):
    file_path = os.path.join(directory_path, vtk_file)
    u, v, w   = getVector(file_path,Nx,Ny,Nz,name)
    um += u
    vm += v
    wm += w
  n  = float(len(vtk_files))
  um = um / n
  vm = vm / n
  wm = wm / n
  return um, vm, wm


def getMeanScalar(directory_path,vtk_files,Nx,Ny,Nz,name):
  am = np.zeros((Nz,Ny,Nx), dtype=np.float32)

  for vtk_file in tqdm(vtk_files):
    file_path = os.path.join(directory_path, vtk_file)
    a  = getScalar(file_path,Nx,Ny,Nz,name)
    am += a
  n  = float(len(vtk_files))
  am = am / n
  return am


class Vector_Data:
  def __init__(self, dir):
    self.dir = dir
    files = [f for f in os.listdir(self.dir) if f.endswith(".vtr")]
    files.sort(key=extract_number)
    self.files = files
    file = os.path.join(self.dir, files[0])
    self.Nx, self.Ny, self.Nz, self.X, self.Y, self.Z = getGrid(file)

  def getVector(self, file_path, name):
    reader = vtk.vtkXMLRectilinearGridReader()
    reader.SetFileName(file_path)
    reader.Update()

    # get dataset
    Q = reader.GetOutput()
    V = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray(name))

    u = np.reshape(V[:,0], [self.Nz,self.Ny,self.Nx])
    v = np.reshape(V[:,1], [self.Nz,self.Ny,self.Nx])
    w = np.reshape(V[:,2], [self.Nz,self.Ny,self.Nx])
    return np.float32(u), np.float32(v), np.float32(w)
  
  def interp_xz(self, Lx1, Lx2, Ly1, Ly2, stridex, stridez, ny):
    nx1 = int(Lx1 / self.X[-1] * self.Nx)
    nx2 = int(Lx2 / self.X[-1] * self.Nx)
    ix  = np.arange(nx1,  nx2, stridex)
    iy  = np.arange(0, len(self.Y), 1)
    iz  = np.arange(0, len(self.Z), stridez)
    x   = self.X[ix]
    z   = self.Z[iz]
    y   = np.linspace(Ly1, Ly2, ny)
    ix, iy, iz = np.meshgrid(ix, iy, iz)
    return ix, iy, iz, x, y, z

  def getVector_interp(self, file_path, name, Lx1, Lx2, Ly1, Ly2, stridex, stridez, ny):
    U, V, W = self.getVector(file_path, name)
    ix, iy, iz, x, y, z = self.interp_xz(Lx1, Lx2, Ly1, Ly2, stridex, stridez, ny)
    u = U[iz,iy,ix].transpose(2,0,1)
    v = V[iz,iy,ix].transpose(2,0,1)
    w = W[iz,iy,ix].transpose(2,0,1)
    nx, nz = len(x), len(z)
    ui = np.zeros((nz,ny,nx), dtype=np.float32)
    vi = np.zeros((nz,ny,nx), dtype=np.float32)
    wi = np.zeros((nz,ny,nx), dtype=np.float32)
    for k in range(nz):
      for i in range(nx):
        ui[k,:,i] = np.interp(y, self.Y, u[k,:,i])
        vi[k,:,i] = np.interp(y, self.Y, v[k,:,i])
        wi[k,:,i] = np.interp(y, self.Y, w[k,:,i])
    return np.float32(ui), np.float32(vi), np.float32(wi), x, y, z
  
  def getMeanVector_interp(self, name, Lx1, Lx2, Ly1, Ly2, stridex, stridez, ny, path=None):
    if path is not None:
      umean, vmean, wmean, x, y, z = self.getVector_interp(path, 'velocity', Lx1, Lx2, Ly1, Ly2, stridex, stridez, ny)
    else:
      ix, iy, iz, x, y, z = self.interp_xz(Lx1, Lx2, Ly1, Ly2, stridex, stridez, ny)
      nx, ny, nz = len(x), len(y), len(z)
      umean = np.zeros((nz,ny,nx), dtype=np.float32)
      vmean = np.zeros((nz,ny,nx), dtype=np.float32)
      wmean = np.zeros((nz,ny,nx), dtype=np.float32)
      for file in tqdm(self.files):
        file_path = os.path.join(self.dir, file)
        u, v, w, x, y, z = self.getVector_interp(file_path, name, Lx1, Lx2, Ly1, Ly2, stridex, stridez, ny)
        umean += u
        vmean += v
        wmean += w
      umean /= len(self.files)
      vmean /= len(self.files)
      wmean /= len(self.files)
    return umean, vmean, wmean, x, y, z


class Data:
  def __init__(self, dir, endT, Lx1=None, Lx2=None, Ly1=None, Ly2=None, Lz1=None, Lz2=None, \
               stridex=1, stridey=1, stridez=1):
    self.dir  = dir
    self.endT = endT
    
    if Lx1 is not None:
      self.Lx1 = Lx1
    else:
      self.Lx1 = None
    if Lx2 is not None:
      self.Lx2 = Lx2
    else:
      self.Lx2 = None
    if Ly1 is not None:
      self.Ly1 = Ly1
    else:
      self.Ly1 = None
    if Ly2 is not None:
      self.Ly2 = Ly2
    else:
      self.Ly2 = None
    if Lz1 is not None:
      self.Lz1 = Lz1
    else:
      self.Lz1 = None
    if Lz2 is not None:
      self.Lz2 = Lz2
    else:
      self.Lz2 = None

    self.stridex = stridex
    self.stridey = stridey
    self.stridez = stridez


  @staticmethod
  def indices(X, stride, L1=None, L2=None):
    def index(x, L):
      if L >= x[-1]:
        n = len(x)
      elif L <= x[0]:
        n = 0
      else:
        for i in range(len(x)):
          if L < x[i]:
            n = i
            break
      return n
    
    if L1 is None:
      n1 = 0
    else:
      n1 = index(X, L1)
    if L2 is None:
      n2 = len(X)
    else:
      n2 = index(X, L2)
    return np.arange(n1, n2, stride)
 

  def make_grid(self):
    files = [f for f in os.listdir(self.dir) if f.endswith(".vtr")]
    files.sort(key=extract_number)
    Nx, Ny, Nz, X, Y, Z = getGrid(os.path.join(self.dir, files[0]))
    t = np.linspace(0.e0, self.endT, len(files))
      
    if self.Lx1 is None:
      if self.Lx2 is None:
        indicesx = np.arange(0, len(X), self.stridex)
      else:
        indicesx = self.indices(X, self.stridex, L1=None, L2=self.Lx2)
    else:
      if self.Lx2 is None:
        indicesx = self.indices(X, self.stridex, L1=self.Lx1, L2=None)
      else:
        indicesx = self.indices(X, self.stridex, self.Lx1, self.Lx2)
    x = X[indicesx]
    
    if self.Ly1 is None:
      if self.Ly2 is None:
        indicesy = np.arange(0, len(Y), self.stridey)
      else:
        indicesy = self.indices(Y, self.stridey, L1=None, L2=self.Ly2)
    else:
      if self.Ly2 is None:
        indicesy = self.indices(Y, self.stridey, L1=self.Ly1, L2=None)
      else:
        indicesy = self.indices(Y, self.stridey, self.Ly1, self.Ly2)
    y = Y[indicesy]
    
    if self.Lz1 is None:
      if self.Lz2 is None:
        indicesz = np.arange(0, len(Z), self.stridez)
      else:
        indicesz = self.indices(Z, self.stridez, L1=None, L2=self.Lz2)
    else:
      if self.Lz2 is None:
        indicesz = self.indices(Z, self.stridez, L1=self.Lz1, L2=None)
      else:
        indicesz = self.indices(Z, self.stridez, self.Lz1, self.Lz2)
    z = Z[indicesz]
    indicesx, indicesy, indicesz = np.meshgrid(indicesx, indicesy, indicesz)
    return Nx, Ny, Nz, indicesx, indicesy, indicesz, x, y, z, t

