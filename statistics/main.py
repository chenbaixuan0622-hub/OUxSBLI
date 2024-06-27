import os
import time
import pickle
import numpy as np
from scipy.stats import gaussian_kde
from scipy.fft import fft
import matplotlib.pyplot as plt

# mylibs
import readVTK
import turbulent_statistics as ts
import myplot as myplt


start_time = time.time()

# search vtk files
directory_path = "./data/"
vtk_files = [f for f in os.listdir(directory_path) if f.endswith(".vtr")]
num_files = len(vtk_files)

# get grid information
first_path = os.path.join(directory_path, vtk_files[0])
Nx, Ny, Nz, x, y, z = readVTK.getGrid(first_path)

# check save data
means = ["./data/umean.npy", "./data/vmean.npy", "./data/wmean.npy"]
if all(os.path.exists(mean) for mean in means):
  # read pickle
  umean = np.load("./data/umean.npy")
  vmean = np.load("./data/vmean.npy")
  wmean = np.load("./data/wmean.npy")
else:
  # get mean velocity
  umean, vmean, wmean = readVTK.getMean(directory_path,vtk_files,Nx,Ny,Nz)
  # save mean velocity
  np.save("./data/umean.npy", umean)
  np.save("./data/vmean.npy", vmean)
  np.save("./data/wmean.npy", wmean)

# plot mean velocity
myplt.plot_velocity(x,y,umean[int(0.5*Nz),:,:],"./data/umean.png")
myplt.plot_velocity(x,y,vmean[int(0.5*Nz),:,:],"./data/vmean.png")
myplt.plot_velocity(x,y,wmean[int(0.5*Nz),:,:],"./data/wmean.png")

# calc integral length scale
# reference point
Nx1 = 100
Ny1 = 20

length = 200
span = Nz-1
# velocity time series data
ut = np.zeros((num_files,span,length))
vt = np.zeros((num_files,span,length))
wt = np.zeros((num_files,span,length))
# fluctuating velocity time series data
uf = np.zeros((num_files,span,length))
vf = np.zeros((num_files,span,length))
wf = np.zeros((num_files,span,length))
# velocity gradient time series data
du = np.zeros((num_files,span,length))
dv = np.zeros((num_files,span,length))
dw = np.zeros((num_files,span,length))

t = 0
for vtk_file in vtk_files:
  # get path
  file_path = os.path.join(directory_path, vtk_file)
  u, v, w = readVTK.getVelocity(file_path,Nx,Ny,Nz)

  # time series data
  for j in range(span):
    for i in range(length):
      ut[t,j,i] = u[j,Ny1,i+Nx1]
      ut[t,j,i] = v[j,Ny1,i+Nx1]
      wt[t,j,i] = w[j,Ny1,i+Nx1]
      uf[t,j,i] = ut[t,j,i] - umean[j,Ny1,i+Nx1]
      vf[t,j,i] = vt[t,j,i] - vmean[j,Ny1,i+Nx1]
      wf[t,j,i] = wt[t,j,i] - wmean[j,Ny1,i+Nx1]
  t += 1
  print("read No.", vtk_file, "file")

t0 = time.time()
# numpy
R11, R22 = ts.longitudinal_and_lateral_corr(span,length,Nx1,Ny1,x,y,uf,vf)
L11, L22 = ts.integral_scale(x[Nx1:Nx1+length-1],R11,R22)

t1 = time.time()
print("correlation and scale calculation time:", t1 - t0)

print("longitudinal integral scale is",L11)
print("lateral integral scale is",L22)
myplt.plot_corr(x,Nx1,length,R11,R22)

# calculate probability density function
# PDF for 1 point velocity profile
u_range, u_pdf = ts.PDF(ut[:,int(0.5*Nz),Nx1])
myplt.plot_pdf(u_range,u_pdf,"u")

v_range, v_pdf = ts.PDF(vt[:,int(0.5*Nz),Nx1])
myplt.plot_pdf(v_range,v_pdf,"v")

w_range, w_pdf = ts.PDF(wt[:,int(0.5*Nz),Nx1])
myplt.plot_pdf(w_range,w_pdf,"w")

# PDF for 1 point fluctuating velocity profile
u_range, u_pdf = ts.PDF(uf[:,int(0.5*Nz),Nx1])
myplt.plot_pdf(u_range,u_pdf,"uf")

v_range, v_pdf = ts.PDF(vf[:,int(0.5*Nz),Nx1])
myplt.plot_pdf(v_range,v_pdf,"vf")

w_range, w_pdf = ts.PDF(wf[:,int(0.5*Nz),Nx1])
myplt.plot_pdf(w_range,w_pdf,"wf")

end_time = time.time()

# write output file
with open("./data/output.txt", "w", encoding="UTF-8") as f:
  print("grid information",file=f)
  print("Nx =",Nx,file=f)
  print("Ny =",Ny,file=f)
  print("Nz =",Nz,file=f)
  print("longitudinal integral scale is",L11,file=f)
  print("lateral integral scale is",L22,file=f)
  print("elapsed time:", end_time - start_time,file=f)

