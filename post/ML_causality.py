import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import optim
from pytorch_memlab import profile
from sklearn.model_selection import train_test_split
from tqdm import tqdm
import os
import sys
import cv2
import matplotlib.pyplot as plt
# AI modules
from mod_AI.plot  import plot_Dataset, plot_loss, plot_result
from mod_AI.utils import Dataset, divide_into_batch, trainNN 
from mod_AI.block import noBlock, BasicBlock, DenseBlock, SpectralTransform, FourierBlock
from mod_AI.cnn   import ResNet, DenseNet, UNet
# vtk modules
from mod.mod_read import getGrid, getVector, getScalar
from mod.mod_turb_stat import tau_2d


Q_directory = "../../../../../media/user/HD-EDS-E/hatayama/TBL/TBL20240819_KEEP6thVisc4th"
Q_files     = [f for f in os.listdir(Q_directory) if f.endswith("vtr")]
num_files   = len(Q_files)

AI_directory = os.path.join(Q_directory, "AI")
os.makedirs(AI_directory, exist_ok = True)

# read grid information
first_path = os.path.join(Q_directory, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)
yp_path = os.path.join(Q_directory, "yp.npy")
yp      = np.load(yp_path)

# teaching_data (viscous layer)
yp1 = 5
# test_data (log layer)
yp2 = 7

# calc yp
for j in range(Ny):
  if yp[j] > yp1:
    Ny_target1 = j
    break

for j in range(Ny):
  if yp[j] > yp2:
    Ny_target2 = j
    break

print("yp =", yp1, " when j = ", Ny_target1)
print("yp =", yp2, " when j = ", Ny_target2)

nx = 32
nz = 32

# set GPU and NN property
device       = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
#net_instance = UNet(FourierBlock, [64,32,16,8]).to(device)
net_instance = ResNet(BasicBlock, [64,32,16,8]).to(device)
net          = torch.jit.script(net_instance).to(device)
criterion    = nn.L1Loss()#nn.MSELoss()
optimizer    = optim.Adam(net.parameters())
epoch        = 20
print(net)
print("Device: {}".format(device))

num_split = 100

um_path = os.path.join(Q_directory, "um.npy"  )
vm_path = os.path.join(Q_directory, "vm.npy"  )
wm_path = os.path.join(Q_directory, "wm.npy"  )
um = np.load(um_path)
vm = np.load(vm_path)
wm = np.load(wm_path)
n1 = nz * np.floor(Nz / nz).astype(int)
n2 = nx * np.floor(Nx / nx).astype(int)
kp = np.zeros((num_files // num_split,2,n1,n2), dtype=np.float32)

for i in range(num_split):
  filename  = "model" + str(i-1) + ".pth"
  file_path = os.path.join(AI_directory, filename)

  if os.path.isfile(file_path):
    net.load_state_dict(torch.load(file_path))
    print("load previous parameters")

  start = int(num_files / num_split) * i
  end   = int(num_files / num_split) * (i + 1)
  print("read from ", start, " to ", end)

  itr = 0
  print("reading VTK files")
  Q = np.zeros((5,Nz,Ny,Nx), dtype=np.float32)
  for Q_file in tqdm(Q_files[start:end]):
    file_path = os.path.join(Q_directory, Q_file)
    Q[0,:,:,:]                         = getScalar(file_path, Nx, Ny, Nz, "rho")
    Q[1,:,:,:], Q[2,:,:,:], Q[3,:,:,:] = getVector(file_path, Nx, Ny, Nz, "velocity")
    Q[4,:,:,:]                         = getScalar(file_path, Nx, Ny, Nz, "p")
    tau, _ = tau_2d(Q, x, y, z)
    ut = np.sqrt(np.mean(tau) / np.mean(Q[0,:,0,:]))

    kp[itr,0,:,:] = 0.5e0 * ((Q[1,:n1,Ny_target1,:n2] - um[:n1,Ny_target1,:n2])**2 \
                           + (Q[2,:n1,Ny_target1,:n2] - vm[:n1,Ny_target1,:n2])**2 \
                           + (Q[3,:n1,Ny_target1,:n2] - wm[:n1,Ny_target1,:n2])**2) / ut**2
    kp[itr,1,:,:] = 0.5e0 * ((Q[1,:n1,Ny_target2,:n2] - um[:n1,Ny_target2,:n2])**2 \
                           + (Q[2,:n1,Ny_target2,:n2] - vm[:n1,Ny_target2,:n2])**2 \
                           + (Q[3,:n1,Ny_target2,:n2] - wm[:n1,Ny_target2,:n2])**2) / ut**2

    itr += 1

  teaching_data = kp[:,0,:,:].reshape([-1,1,nz,nx])
  test_data     = kp[:,1,:,:].reshape([-1,1,nz,nx])
  
  print(len(teaching_data))
  print(len(test_data))

  # normalize data
  teaching_mean = np.mean(teaching_data)
  test_mean     = np.mean(test_data)
  teaching_std  = np.std(teaching_data)
  test_std      = np.std(test_data)
  teaching_data = (teaching_data - teaching_mean) / teaching_std
  test_data     = (test_data     - test_mean)     / test_std
  
  del Q

  batchsize = min(len(teaching_data[:,0,0,0]), 200)

  # dataset is torch tensor
  dataset                     = Dataset(teaching_data, test_data)
  train_dataset, test_dataset = train_test_split(dataset, test_size=0.3, shuffle=False)

  # make dataloader
  #plot_Dataset(1e3*x[:nx],1e3*z[:nz],train_dataset[10],AI_directory)
  train_batch, test_batch = divide_into_batch(train_dataset,test_dataset,batchsize)

  net, train_loss_list, test_loss_list = trainNN(net,device,optimizer,criterion,train_batch,test_batch,epoch)
  plot_loss(epoch,train_loss_list,test_loss_list,AI_directory,i)
  filename = "model" + str(i) + ".pth"
  torch.save(net.state_dict(), os.path.join(AI_directory, filename))

  plot_result(1e3*x[:nx],1e3*z[:nz],net,device,test_batch,AI_directory,'kp_log','kp_vis',i,\
  teaching_mean,test_mean,teaching_std,test_std)

