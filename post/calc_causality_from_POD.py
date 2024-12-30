import numpy as np
import pandas as pd
import os
from scipy.signal import correlate, correlation_lags
import matplotlib.pyplot as plt
from lingam import DirectLiNGAM, ICALiNGAM
from statsmodels.tsa.stattools import grangercausalitytests
from cdt.causality.graph import SAM
import networkx as nx
from mod.mod_info import TE, Hxy


# for plot
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['font.size'] = 12

# parameter
endT  = 0.1e-3
Q_dir = "../../data"
name1 = "POD_time_coef_SW.npy"
name2 = "POD_time_coef_BL.npy"
mode  = 4

# data [time, mode]
path1 = os.path.join(Q_dir, name1)
data1 = np.load(path1)
path2 = os.path.join(Q_dir, name2)
data2 = np.load(path2)
data  = np.concatenate([data1[:,:mode], data2[:,:mode]], axis=-1)


# Transfer Entropy
Ate = np.zeros((2*mode, 2*mode), dtype=np.float32)

for j in range(mode*2):
  for i in range(mode*2):
    if i > j:
      xt = data[:,i]
      yt = data[:,j]
      Ate[i,j], Ate[j,i] = TE(xt, yt, 10, 3)
    else:
      Ate[i,i] = np.nan


# Conditional Entropy
Ah = np.zeros((2*mode, 2*mode), dtype=np.float32)

for j in range(mode*2):
  for i in range(mode*2):
    if i > j:
      xt = data[:,i]
      yt = data[:,j]
      Ah[i,j], Ah[j,i] = Hxy(xt, yt, 5)
    else:
      Ah[i,i] = np.nan


# LiNGAM
model = DirectLiNGAM()
model.fit(data)
Alingam = model.adjacency_matrix_
np.fill_diagonal(Alingam, np.nan)


# correlation
Ac = np.zeros((2*mode, 2*mode), dtype=np.float32)

for j in range(mode*2):
  for i in range(mode*2):
    if i > j:
      xt   = data[:,i] - data[:,i].mean()
      yt   = data[:,j] - data[:,j].mean()
      Cxy  = correlate(xt, yt, mode='full')
      Cxy /= (np.linalg.norm(xt, ord=2) * np.linalg.norm(yt, ord=2))
      center = len(xt) - 1
      stride = 10
      Ac[i,j] = np.max(Cxy[center:center+stride])
      Ac[j,i] = np.max(Cxy[center-stride:center])
    else:
      Ac[i,i] = np.nan


'''
# ICALiNGAM
model = ICALiNGAM()
model.fit(data)
Aicalingam = model.adjacency_matrix_
np.fill_diagonal(Aicalingam, np.nan)


# Granger Causality
Agc = np.zeros((2*mode, 2*mode), dtype=np.float32)
max_lag = 10

for i in range(len(Agc[0,:])):
  for j in range(len(Agc[0,:])):
    if i != j:
      test_result = grangercausalitytests(data[:, [i,j]], max_lag, verbose=False)
      p_values = [test[0]['ssr_ftest'][1] for test in test_result.values()]
      Agc[i,j] = -np.log10(min(p_values)) if min(p_values) > 0 else np.nan
    else:
      Agc[i,j] = np.nan


# SAM
model = SAM()
data_df = pd.DataFrame(data, columns=[f'Features_{i}' for i in range(data.shape[1])])
sam_graph = model.predict(data_df)
Asam = nx.to_numpy_array(sam_graph)
np.fill_diagonal(Asam, np.nan)
'''

# plot result
plt.subplot(131)
plt.imshow(Ate, cmap='jet', extent=None)
plt.xlabel("Effect", fontsize=16)
plt.ylabel("Cause", fontsize=16)
plt.title("TE")

'''
plt.subplot(142)
plt.imshow(Ah, cmap='jet', extent=None)
plt.xlabel("Effect", fontsize=16)
plt.ylabel("Cause", fontsize=16)
plt.title("Conditional Entropy")

plt.subplot(143)
plt.imshow(Alingam, cmap='jet', extent=None)
plt.xlabel("Effect", fontsize=16)
plt.ylabel("Cause", fontsize=16)
plt.title("LiNGAM")
'''

plt.subplot(132)
plt.imshow(Ac, cmap='jet', extent=None)
plt.xlabel("Effect", fontsize=16)
plt.ylabel("Cause", fontsize=16)
plt.title("Corr")

plt.subplot(133)
plt.imshow(Ac*Ate, cmap='jet', extent=None)
plt.xlabel("Effect", fontsize=16)
plt.ylabel("Cause", fontsize=16)
plt.title("TE * Corr")

'''
plt.subplot(154)
plt.imshow(Agc, cmap='jet', extent=None)
plt.xlabel("Effect", fontsize=16)
plt.ylabel("Cause", fontsize=16)
plt.title("Granger Causality")

plt.subplot(155)
plt.imshow(Asam, cmap='jet', extent=None)
plt.xlabel("Effect", fontsize=16)
plt.ylabel("Cause", fontsize=16)
plt.title("SAM")
'''
plt.show()

