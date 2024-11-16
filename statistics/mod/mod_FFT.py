import numpy as np

def PSD(data, t):
  n    = len(t)
  fft  = np.fft.fft(data)
  freq = np.fft.fftfreq(n, d=-t[0]+t[1])
  Amp  = abs(fft / (0.5e0 * float(n)))
  ps   = Amp**2
  psd  = ps / float(n)
  return freq, psd

