import torch
import torch.nn as nn

def notears_constr(adj_m, max_pow=None):
  m_exp = [adj_m]
  if max_pow is None:
    max_pow = adj_m.shape[1]
  while(m_exp[-1].sum() > 0 and len(m_exp) < max_pow):
    m_exp.append(m_exp[-1] @ adj_m/len(m_exp))

  return sum([i.diag().sum() for idx, i in enumerate(m_exp)])

