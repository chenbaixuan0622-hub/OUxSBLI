import numpy as np
from scipy.optimize import bisect
from oblique_shock import beta, downstream_mach


M1 = 2.0
theta = 8.0

beta_i = beta(M1, theta)
M2 = downstream_mach(M1, beta_i, theta)
beta_r = beta(M2, theta)

print(f"beta_i = {beta_i:.3f} deg")
print(f"beta_r = {beta_r:.3f} deg")
