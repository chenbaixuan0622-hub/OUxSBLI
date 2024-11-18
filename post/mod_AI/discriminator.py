import torch
import torch.nn as nn

class SAMDiscriminator(nn.Module):
  def __init__(self, nfeatures, dnh, hlayers):
    super(SAMDiscriminator, self).__init__()

    self.nfeatures = nfeatures

    layers = []
    layers.append(nn.Linear(nfeatures, dnh))
    layers.append(nn.BatchNorm1d(dnh))
    layers.append(nn.LeakyReLU(.2))

    for i in range(hlayers-1):
      layers.append(nn.Linear(dnh, dnh))
      layers.append(nn.BatchNorm1d(dnh))
      layers.append(nn.LeakyReLU(.2))

    layers.append(nn.Linear(dnh, 1))

    self.layers = nn.Sequential(*layers)

    mask = torch.eye(nfeatures, nfeatures)
    self.register_buffer("mask", mask.unsqueeze(0))

  def forward(self, input, obs_data=None):
    if obs_data is not None:
      return [self.layers(i) for i in torch.unbind(
        obs_data.unsqueeze(1) * (1 - self.mask)
          + input.unsqueeze(1) * self.mask, 1)]
    else:
      return self.layers(input)
  
  def reset_parameters(self):
    for layer in self.layers:
      if hasattr(layer, 'reset_parameters'):
        layer.reset_parameters()

