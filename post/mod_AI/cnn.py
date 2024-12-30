import torch
import torch.nn as nn
import torch.nn.functional as F
#from group.p4m import RotatedConv1, RotatedConv, RotatedTConv, GPReLU, GSequential, Gdown, GBatchNorm2d, GMaxPool2d


class ResNet(nn.Module):
  def __init__(self, block, layers, num_classes=1000):
    super().__init__()

    self.in_channels = 64
    self.conv1       = nn.Conv2d(2, 64, kernel_size=7, stride=2, padding=3)
    self.layer1      = self._make_layer(block, 64, layers[0], stride=1)
    self.layer2      = self._make_layer(block, 32, layers[1], stride=1)
    self.layer3      = self._make_layer(block, 32, layers[2], stride=1)
    self.layer4      = self._make_layer(block, 2,  layers[3], stride=1)
    self.maxpool     = nn.MaxPool2d(kernel_size=2, stride=1, padding=1)
    self.dropout     = nn.Dropout(0.05)

  def _make_layer(self, block, channels, blocks, stride):
    layers = []

    # 1st Residual Block
    layers.append(block(self.in_channels, channels, stride))

    # the other Residual Blocks
    self.in_channels = channels
    for _ in range(1, blocks):
      layers.append(block(channels, channels))

    return nn.Sequential(*layers)

  def forward(self, x):
    input_size = x[0,0,:,:].size()
    X = self.conv1(x)
    X = self.layer1(X)
    X = self.layer2(X)
    X = self.layer3(X)
    X = self.layer4(X)
    x = self.maxpool(X)
    x = self.dropout(x)
    x = F.interpolate(x,size=input_size, mode='bicubic',align_corners=False)
    return x


class DenseNet(nn.Module):
  def __init__(self, block, layers, num_classes=1000):
    super().__init__()

    self.in_channels = 64
    self.conv1       = nn.Conv2d(2, 64, kernel_size=7, stride=2, padding=3)
    self.layer1      = self._make_layer(block, 64, layers[0], stride=1)
    self.bn1         = nn.BatchNorm2d(64)
    self.layer2      = self._make_layer(block, 32, layers[1], stride=1)
    self.bn2         = nn.BatchNorm2d(32)
    self.layer3      = self._make_layer(block, 32, layers[2], stride=1)
    self.bn3         = nn.BatchNorm2d(32)
    self.layer4      = self._make_layer(block, 2,  layers[3], stride=1)
    self.bn4         = nn.BatchNorm2d(2)
    self.maxpool     = nn.MaxPool2d(kernel_size=2, stride=1, padding=1)
    self.dropout     = nn.Dropout(0.1)

  def _make_layer(self, block, channels, blocks, stride):
    layers = []
    # 1st Residual Block
    layers.append(block(self.in_channels, channels, stride))
    # the other Residual Blocks
    self.in_channels = channels
    for _ in range(1, blocks):
      layers.append(block(channels, channels))
    return nn.Sequential(*layers)

  def forward(self, x):
    input_size = x[0,0,:,:].size()
    X = self.conv1(x)
    X = self.layer1(X)
    X = self.bn1(X)
    X = self.layer2(X)
    X = self.bn2(X)
    X = self.layer3(X)
    X = self.bn3(X)
    X = self.layer4(X)
    X = self.bn4(X)
    x = self.maxpool(X)
    x = self.dropout(x)
    x = F.interpolate(x,size=input_size, mode='bicubic',align_corners=False)
    return x


class UNet(nn.Module):
  def __init__(self, block, num_of_ch):
    super().__init__()
    channels1, channels2, channels3, channels4 = num_of_ch
    # encoder
    # layer1
    self.conv1 = nn.Conv2d(1, channels1, kernel_size=3, stride=1, padding=0)
    self.relu1 = nn.PReLU()
    # layer2
    self.conv2 = nn.Conv2d(channels1, channels2, kernel_size=3, stride=1, padding=0)
    self.relu2 = nn.PReLU()
    # layer3
    self.conv3 = nn.Conv2d(channels2, channels3, kernel_size=3, stride=1, padding=0)
    self.relu3 = nn.PReLU()
    # layer4
    self.conv4 = nn.Conv2d(channels3, channels4, kernel_size=3, stride=1, padding=0)
    self.relu4 = nn.PReLU()

    # Block
    self.layer = self._make_layer(block, channels4)
    
    # decoder
    # layer5
    self.tconv1   = nn.ConvTranspose2d(channels4, channels3, kernel_size=3, stride=1, padding=0)
    self.relu5    = nn.PReLU()
    # layer6
    self.tconv2   = nn.ConvTranspose2d(channels3, channels2, kernel_size=3, stride=1, padding=0)
    self.relu6    = nn.PReLU()
    # layer7
    self.tconv3   = nn.ConvTranspose2d(channels2, channels1, kernel_size=3, stride=1, padding=0)
    self.relu7    = nn.PReLU()
    # layer8
    self.tconv4   = nn.ConvTranspose2d(channels1, 1, kernel_size=3, stride=1, padding=0)
    self.relu8    = nn.PReLU()

    self.connect  = nn.Sequential()
    self.dropout  = nn.Dropout(0.1)

  def _make_layer(self, block, in_channels):
    layers = []
    layers.append(block(in_channels))
    return nn.Sequential(*layers)
  
  def forward(self, x):
    # encoder
    # layer1
    x1 = self.conv1(x)
    x1 = self.relu1(x1)
    # layer2
    x2 = self.conv2(x1)
    x2 = self.relu2(x2)
    # layer3
    x3 = self.conv3(x2)
    x3 = self.relu3(x3)
    # layer4
    x4 = self.conv4(x3)
    x4 = self.relu4(x4)

    # block
    x4f = self.layer(x4)
    x4  = self.layer(x4f) + x4f

    # decoder
    # layer5
    x5 = self.tconv1(x4)
    x5 = self.relu5(x5) + x3
    # layer6
    x6 = self.tconv2(x5)
    x6 = self.relu6(x6) + x2
    # layer7
    x7 = self.tconv3(x6)
    x7 = self.relu7(x7) + x1
    # layer8
    x8 = self.tconv4(x7)
    x  = self.relu8(x8)

    # dropout
    x = self.dropout(x)
    return x


class Encoder(nn.Module):
  def __init__(self, block, num_of_ch, size_flatten):
    super().__init__()
    channels1, channels2, channels3, channels4 = num_of_ch
    # layer1
    self.conv1 = nn.Conv2d(1, channels1, kernel_size=3, stride=1, padding=0)
    self.relu1 = nn.PReLU()
    # layer2
    self.conv2 = nn.Conv2d(channels1, channels2, kernel_size=3, stride=1, padding=0)
    self.relu2 = nn.PReLU()
    # layer3
    self.conv3 = nn.Conv2d(channels2, channels3, kernel_size=3, stride=1, padding=0)
    self.relu3 = nn.PReLU()
    # layer4
    self.conv4 = nn.Conv2d(channels3, channels4, kernel_size=3, stride=1, padding=0)
    self.relu4 = nn.PReLU()
    # Block
    self.layer = self._make_layer(block, channels4)
    # dropout
    self.dropout = nn.Dropout(0.1)
    # linear
    self.fc = nn.Linear(size_flatten, 10)

  def _make_layer(self, block, in_channels):
    layers = []
    layers.append(block(in_channels))
    return nn.Sequential(*layers)
  
  def forward(self, x):
    # layer1
    x = self.conv1(x)
    x = self.relu1(x)
    # layer2
    x = self.conv2(x)
    x = self.relu2(x)
    # layer3
    x = self.conv3(x)
    x = self.relu3(x)
    # layer4
    x = self.conv4(x)
    x = self.relu4(x)
    # block
    x = self.layer(x)

    # dropout
    x = self.dropout(x)
    x = x.view(x.size(0), -1)
    #print(x.size())
    x = self.fc(x)
    return x

