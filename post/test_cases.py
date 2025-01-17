import numpy as np
import matplotlib.pyplot as plt
from mod.mod_plot import set_Params
from mod.mod_info import EE
from mod.mod_ds import case1, case2, case3, case4

set_Params()


nt  = 1000
w   = 0.1e0 

# case 1
x, y = case1(nt)
EE1xy = EE(x, y, p=5)
EE1yx = EE(y, x, p=5)

print("case 1 x->y:", EE1xy, "y->x", EE1yx)

# case 2
x, y = case2(nt, w)
EE2xy = EE(x, y, p=5)
EE2yx = EE(y, x, p=5)

print("case 2 x->y:", EE2xy, "y->x", EE2yx)

# case 3
x, y = case3(nt, w)
EE3xy = EE(x, y, p=5)
EE3yx = EE(y, x, p=5)

print("case 3 x->y:", EE3xy, "y->x", EE3yx)

# case 4
x, y, z = case4(nt, w)
EE4xy = EE(x, y, p=7)
EE4yx = EE(y, x, p=7)
EE4yz = EE(y, z, p=7)
EE4zy = EE(z, y, p=7)
EE4zx = EE(z, x, p=7)
EE4xz = EE(x, z, p=7)

print("case 4 x->y:", EE4xy, "y->x", EE4yx)
print("case 4 y->z:", EE4yz, "z->y", EE4zy)
print("case 4 z->x:", EE4zx, "x->z", EE4xz)

