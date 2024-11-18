import torch
import torch.nn as nn
from cdt.utils.torch import ChannelBatchNorm1d, MatrixSampler, Linear3D

class SAMGenerator(nn.Module):
  def __init__(self, data_shape, nh):
    super(SAMGenerator, self).__init__()
    
    nb_vars = data_shape[1]
    skeleton = 1 - torch.eye(nb_vars + 1, nb_vars)

    self.register_buffer('skeleton', skeleton)

    self.input_layer = Linear3D(
      (nb_vars, nb_vars + 1, nh))

    layers = []

    layers.append(ChannelBatchNorm1d(nb_vars, nh))
    layers.append(nn.Tanh())
    self.layers = nn.Sequential(*layers)

    self.output_layer = Linear3D((nb_vars, nh, 1))

  def forward(self, data, noise, adj_matrix, drawn_neurons=None):
    x = self.input_layer(data, noise, adj_matrix * self.skeleton)

    x = self.layers(x)

    output = self.output_layer(x, noise=None, adj_matrix=drawn_neurons)

    return output.squeeze(2)

  def reset_parameters(self):
    self.input_layer.reset_parameters()
    self.output_layer.reset_parameters()

    for layer in self.layers:
      if hasattr(layer, 'reset_parameters'):
        layer.reset_parameters()

