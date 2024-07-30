import numpy as np
from pyinform.dist import Dist
from pyinform.shannon import entropy

N = 1000

obs = np.zeros(N, dtype=np.float32)
for i in range(N):
  obs[i] = np.float32(np.random.randint(1, 7))

px, x = np.histogram(obs, bins=6)

d = Dist(px)

print(entropy(d))

