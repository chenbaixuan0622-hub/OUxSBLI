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
#my libs
import readVTK
from plot  import plot_Dataset, plot_loss, plot_result
from utils import Dataset, divide_into_batch, trainNN 
from block import noBlock, BasicBlock, DenseBlock, SpectralTransform, FourierBlock
from cnn import ResNet, DenseNet, UNet


# read grid information
dir_path   = "./data/0"
vtk_files  = [f for f in os.listdir(dir_path) if f.endswith("vtr")]
first_path = os.path.join(dir_path, vtk_files[0])
Nx, Ny, Nz, xtemp, ytemp, ztemp = readVTK.getGrid(first_path)
# dataset size
nx = 200
ny = Ny
nz = Nz
# reference point
nxr = 100
nyr = 0
nzr = 0
x = np.zeros((nx), dtype = np.float32)
y = np.zeros((ny), dtype = np.float32)
z = np.zeros((nz), dtype = np.float32)
x = xtemp[nxr:nxr+nx]
y = ytemp[nyr:nyr+ny]
z = ztemp[nzr:nzr+nz]


# set GPU and NN property
device       = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
net_instance = UNet(FourierBlock, [64,32,16,8]).to(device)
net          = torch.jit.script(net_instance).to(device)
criterion    = nn.MSELoss()
optimizer    = optim.Adam(net.parameters())
epoch        = 50
batchsize    = 100
print(net)
print("Device: {}".format(device))

num_dir = 10

for i in range(num_dir):
  # load NN parameters
  pdir_path = os.path.join("./data/", str(i-1))
  if os.path.isfile(os.path.join(pdir_path, "model.pth")):
    net.load_state_dict(torch.load(os.path.join(pdir_path, "model.pth")))
    print("load previous parameters")

  # read VTK files
  dir_path  = os.path.join("./data/", str(i))
  vtk_files = [f for f in os.listdir(dir_path) if f.endswith("vtr")]
  num_files = len(vtk_files)  

  Q = np.zeros((num_files,nz,ny,nx), dtype = np.float32)
  p = np.zeros((num_files,nz,ny,nx), dtype = np.float32)

  itr = 0
  Qtemp = np.zeros((Nz,Ny,Nx), dtype = np.float32) 
  ptemp = np.zeros((Nz,Ny,Nx), dtype = np.float32)
  print("reading VTK files")
  for vtk_file in tqdm(vtk_files):
    file_path = os.path.join(dir_path, vtk_file)
    Qtemp        = readVTK.getPointData(file_path,Nx,Ny,Nz,"Qcriterion")
    Q[itr,:,:,:] = Qtemp[nzr:nzr+nz,nyr:nyr+ny,nxr:nxr+nx]
    ptemp        = readVTK.getPointData(file_path,Nx,Ny,Nz,"p")
    p[itr,:,:,:] = ptemp[nzr:nzr+nz,nyr:nyr+ny,nxr:nxr+nx]
    itr += 1
  del Qtemp
  del ptemp

  teaching_data = np.empty((num_files*nz,1,ny,nx))
  test_data     = np.empty_like(teaching_data)

  # set teaching and test data
  itr = 0
  Qmean = np.mean(Q)
  pmean = np.mean(p)
  Qstd  = np.std(Q)
  pstd  = np.std(p)
  for j in range(num_files):
    for i in range(nz):
      teaching_data[itr,0,:,:] = (Q[j,i,:,:] - Qmean) / Qstd
      test_data[itr,0,:,:]     = (p[j,i,:,:] - pmean) / pstd
      itr += 1

  del Q
  del p

  # dataset is torch tensor
  dataset                     = Dataset(teaching_data, test_data)
  train_dataset, test_dataset = train_test_split(dataset, test_size=0.1, shuffle=False)

  # make dataloader
  plot_Dataset(x,y,train_dataset[50],dir_path)
  train_batch, test_batch = divide_into_batch(train_dataset,test_dataset,batchsize)

  net, train_loss_list, test_loss_list = trainNN(net,device,optimizer,criterion,train_batch,test_batch,epoch)
  plot_loss(epoch,train_loss_list,test_loss_list,dir_path)
  torch.save(net.state_dict(), os.path.join(dir_path, 'model.pth'))

  plot_result(x,y,net,device,test_batch,dir_path,'p','Q')

