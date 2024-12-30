import torch
import torch.nn as nn
import torch.nn.functional as F
#from group.p4m import RotatedConv1, RotatedConv, GPReLU, GSequential, Gdown, GBatchNorm2d, GMaxPool2d

# no block
class noBlock(nn.Module):
  def __init__(self, in_channels):
    super().__init__()

  def forward(self, x):
    return x


# BasicBlock for ResNet
class BasicBlock(nn.Module):
  def __init__(self, in_channels, out_channels, stride=1):
    super().__init__()
    self.conv1    = nn.Conv2d(in_channels,  out_channels, kernel_size=5, stride=1, padding=2)
    self.relu     = nn.PReLU()
    self.conv2    = nn.Conv2d(out_channels, out_channels, kernel_size=5, stride=1, padding=2)
    self.down     = nn.Conv2d(in_channels,  out_channels, kernel_size=1, stride=1, padding=0)
    self.shortcut = nn.Sequential()

  def forward(self, x):
    out = self.conv1(x)
    out = self.relu(out)
    out = self.conv2(out)
    out = self.shortcut(self.down(x),out)
    return out


# DenseBlock for DenseNet
class DenseBlock(nn.Module):
  def __init__(self, in_channels, out_channels, stride=1):
    super().__init__()
    self.down       = nn.Conv2d(in_channels,  out_channels, kernel_size=1, stride=1, padding=0)
    self.conv1      = nn.Conv2d(in_channels,  out_channels, kernel_size=5, stride=1, padding=2)
    self.relu1      = nn.PReLU()
    self.shortcut1  = nn.Sequential()
    self.conv2      = nn.Conv2d(out_channels, out_channels, kernel_size=5, stride=1, padding=2)
    self.relu2      = nn.PReLU()
    self.shortcut2  = nn.Sequential()
    self.conv3      = nn.Conv2d(out_channels, out_channels, kernel_size=5, stride=1, padding=2)
    self.relu3      = nn.PReLU()
    self.shortcut3  = nn.Sequential()
    self.conv4      = nn.Conv2d(out_channels, out_channels, kernel_size=5, stride=1, padding=2)
    self.relu4      = nn.PReLU()
    self.shortcut4  = nn.Sequential()

  def forward(self, x):
    # 1st layer
    out1 = self.conv1(x)
    out1 = self.relu1(out1)
    out1 = self.shortcut1(self.down(x),out1)
    # 2nd layer
    out2 = self.conv2(out1)
    out2 = self.relu2(out2)
    out2 = self.shortcut2(self.down(x),out2)
    out2 = self.shortcut2(out1,out2)
    # 3rd layer
    out3 = self.conv3(out2)
    out3 = self.relu3(out3)
    out3 = self.shortcut3(self.down(x),out3)
    out3 = self.shortcut3(out1,out3)
    out3 = self.shortcut3(out2,out3)
    # 4 th layer
    out4 = self.conv4(out3)
    out4 = self.relu4(out4)
    out4 = self.shortcut4(self.down(x),out4)
    out4 = self.shortcut4(out1,out4)
    out4 = self.shortcut4(out2,out4)
    out4 = self.shortcut4(out3,out4)
    return out4


# SpectralTransform
class SpectralTransform(nn.Module):
  def __init__(self, in_channels, out_channels):
    super().__init__()

    channels   = 2*in_channels
    self.conv1 = nn.Conv2d(in_channels=2*in_channels,out_channels=channels, kernel_size=5, stride=1, padding=0)
    self.relu1 = nn.PReLU()
    self.conv2 = nn.ConvTranspose2d(in_channels=channels,out_channels=2*out_channels, kernel_size=5, stride=1, padding=0)
    self.relu2 = nn.PReLU()

  def forward(self,x):
    # x[batch_size, channel_size, height, width]
    x_size = x[0,0,:,:].size()
    # FFT
    x = torch.fft.rfftn(x, dim=(-2,-1), norm="ortho")
    # stack real part and imag part
    x = torch.cat([x.real, x.imag], dim=1)

    # Fourier Convolution
    x = self.conv1(x)
    x = self.relu1(x)
    x = self.conv2(x)
    x = self.relu2(x)

    # split into 2 parts, real part and imag part
    x_real, x_imag = torch.chunk(x, chunks=2, dim=1)
    # IFFT
    x = torch.fft.irfftn(torch.complex(x_real, x_imag), dim=(-2,-1), norm="ortho")
    x = F.interpolate(x, size=x_size, mode='bicubic', align_corners=False)
    return x


# Fourier Block
class FourierBlock(nn.Module):
  def __init__(self, channels):
    super().__init__()

    mid_channels = channels
    # input
    self.maxpool = nn.MaxPool2d(kernel_size=1, stride=1, padding=0)
    # global
    self.ffc     = SpectralTransform(in_channels=channels, out_channels=channels)
    self.conv1   = nn.Sequential(nn.Conv2d(in_channels=channels, out_channels=mid_channels, kernel_size=5, stride=1, padding=0),nn.PReLU())
    self.conv2   = nn.Sequential(nn.ConvTranspose2d(in_channels=mid_channels, out_channels=channels, kernel_size=5, stride=1, padding=0),nn.PReLU())
    # local
    self.conv3   = nn.Sequential(nn.Conv2d(in_channels=channels, out_channels=mid_channels, kernel_size=3, stride=1, padding=0),nn.PReLU())
    self.conv4   = nn.Sequential(nn.ConvTranspose2d(in_channels=mid_channels, out_channels=channels, kernel_size=3, stride=1, padding=0),nn.PReLU())
    self.conv5   = nn.Sequential(nn.Conv2d(in_channels=channels, out_channels=mid_channels, kernel_size=3, stride=1, padding=0),nn.PReLU())
    self.conv6   = nn.Sequential(nn.ConvTranspose2d(in_channels=mid_channels, out_channels=channels, kernel_size=3, stride=1, padding=0),nn.PReLU())
    # merge
    self.connect = nn.Sequential()
    self.relu7   = nn.PReLU()
    # output
    self.conv    = nn.Conv2d(in_channels=channels, out_channels=channels, kernel_size=1, stride=1, padding=0)

  def forward(self, X):
    #print(X.shape) [group_size(p4m), batch_size, channel_size, height, width]
    x = self.maxpool(X)
    
    # global
    x1 = self.ffc(x)
    x2 = self.conv1(x)
    x2 = self.conv2(x2)

    # local
    x3 = self.conv3(x)
    x3 = self.conv4(x3)
    x4 = self.conv5(x)
    x4 = self.conv6(x4)

    # merge
    # global
    x_global = self.connect(x1)+self.connect(x3)
    x_global = self.relu7(x_global)
    # local
    x_local = self.connect(x2)+self.connect(x4)
    x_local = self.relu7(x_local)
    x = self.connect(x_global)+self.connect(x_local)

    # output
    X = self.conv(x)

    return X

