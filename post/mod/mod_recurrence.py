import numpy as np
from sklearn.metrics.pairwise import pairwise_distances
import matplotlib.pyplot as plt

def recurrence_plot(s, eps=None, steps=None):
  if eps == None: eps=0.1
  if steps == None: steps=10
  d = pairwise_distances(s)
  d = np.floor(d / eps)
  d[d > steps] = steps
  return d

