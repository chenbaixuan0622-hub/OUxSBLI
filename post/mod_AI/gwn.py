import torch
import torch.nn as nn
import torch.nn.functional as F


class gwnet(nn.Module):
  def __init__(self, device, in_dim, out_dim, residual_channels, dilation_channels, skip_channels, end_channels, dropout, nx, nt, Ainit):
    super(gwnet, self).__init__()
    # Ainit : initial matrix
    # E1init, E2init : initial node embedding

    m, p, n = torch.svd(Ainit)
    self.E1 = torch.mm(m[:,:10], torch.diag(p[:10]**0.5))
    self.E2 = torch.mm(torch.diag(p[:10]**0.5), n[:,:10].t())

    self.start_conv = nn.Conv2d(in_channels=in_dim,
                                out_channels=residual_channels,
                                kernel_size=(1,1))

    self.filter_convs = nn.ModuleList()
    self.gate_convs   = nn.ModuleList()

    kernel_size = 2
    self.layers = nt - 1
    for i in range(self.layers):
      # dilated causal convolution
      self.filter_convs.append(nn.Conv2d(in_channels=residual_channels,
                                         out_channels=dilation_channels,
                                         kernel_size=(1,kernel_size),
                                         stride=1, padding=0, dilation=1))
      self.gate_convs.append(nn.Conv2d(in_channels=residual_channels,
                                       out_channels=dilation_channels,
                                       kernel_size=(1,kernel_size),
                                       stride=1, padding=0, dilation=1))

    self.residual_conv = nn.Conv2d(in_channels=dilation_channels,
                                   out_channels=skip_channels,
                                   kernel_size=(1,1))
    self.skip_conv = nn.Conv1d(in_channels=residual_channels,
                               out_channels=skip_channels,
                               kernel_size=(1))

    # graph convolution
    self.weight  = nn.Parameter(torch.randn(nx,1).to(device),  requires_grad=True).to(device)
    self.A       = nn.Parameter(torch.randn(nx,nx).to(device), requires_grad=True).to(device)
    self.dropout = nn.Dropout(dropout)

    self.end_conv_1 = nn.Conv2d(in_channels=skip_channels,
                                out_channels=end_channels,
                                kernel_size=(1,1),
                                bias=True)
    self.end_conv_2 = nn.Conv2d(in_channels=end_channels,
                                out_channels=out_dim,
                                kernel_size=(1,1),
                                bias=True)

  def forward(self, x):
    # x[B,F,S,T]
    x = self.start_conv(x)

    Ainit = F.softmax(F.relu(torch.mm(self.E1, self.E2)), dim=1)
    
    skip = x[:,:,:,-1]
    for i in range(self.layers):
      residual = x

      # dilated causal convolution
      filter = self.filter_convs[i](residual)
      filter = torch.tanh(filter)

      gate   = self.gate_convs[i](residual)
      gate   = torch.sigmoid(gate)
      x      = filter * gate

      x = x + residual[:,:,:,-x.size(3):]

    # graph convolution
    support = self.weight.unsqueeze(0).unsqueeze(0)
    support = x * support
    A = Ainit + self.A
    x = torch.matmul(A, support)
    x = self.dropout(x)

    x = self.residual_conv(x)

    skip = self.skip_conv(skip)
    x = x + torch.reshape(skip, x.shape)

    x = F.relu(x)
    x = F.relu(self.end_conv_1(x))
    x = self.end_conv_2(x)
    return x

