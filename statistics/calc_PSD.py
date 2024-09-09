import numpy as np
import scipy as sp
import os
import matplotlib.pyplot as plt
# read module
from mod.mod_read import getGrid

Qk6_directory = "../../../../../media/user/HD-EDS-E/hatayama/TBL/TBL20240819_KEEP6thVisc4th"
Qs_directory  = "../../../../../media/user/HD-EDS-E/hatayama/TBL/TBL20240822_SLAU"

Qk_files   = [f for f in os.listdir(Qk6_directory) if f.endswith(".vtr")]
num_files = len(Qk_files)
Qs_files   = [f for f in os.listdir(Qs_directory) if f.endswith(".vtr")]
num_files = len(Qs_files)

endT = 0.385384e-3
t    = np.linspace(0.e0, endT, num_files)

def PSD(uf, t, num_files):
  fft  = np.fft.fft(uf[0,:])
  freq = np.fft.fftfreq(num_files, d=-t[0]+t[1])
  Amp  = abs(fft / (0.5e0 * float(num_files)))
  ps   = Amp**2
  psd  = ps / float(num_files)
  return freq, psd

uf_path = os.path.join(Qk6_directory, "uf.npy")
# uf[space, time]
uf = np.load(uf_path)
freq, psd = PSD(uf, t, num_files)

save_path = os.path.join(Qk6_directory, "PSD_u.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# k       PSD", file=f)
  for i in range(1,int(0.5*num_files)):
    print(f'{freq[i]:.3e}', f'{psd[i]:.3e}', file=f)

uf_path = os.path.join(Qs_directory, "uf.npy")
# uf[space, time]
uf = np.load(uf_path)
freq, psd = PSD(uf, t, num_files)

save_path = os.path.join(Qs_directory, "PSD_u.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# k       PSD", file=f)
  for i in range(1,int(0.5*num_files)):
    print(f'{freq[i]:.3e}', f'{psd[i]:.3e}', file=f)



