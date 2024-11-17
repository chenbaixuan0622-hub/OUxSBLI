import numpy as np
import jax
import jax.numpy as jnp
import os
import time
from src.calc_time_dev import Runge_Kutta
from src.print import print_vtk


def main(nx, ny, nz, gamma, R, dt, rho0, M0):
  Nt = 1
  Np = 50 
  dt = 0.01e0


  Q = set_init(nx, ny, nz, gamma, rho0, M0, x, y, z, rho0, M0)
  print_vtk(x, y, z, gamma, Q, 0)

  start = time.time()
  for i in range(Np):
    xs = nx, ny, nz, dxdy, dydz, dzdx, J, Q
    xs = jax.lax.fori_loop(0, Nt, Runge_Kutta, xs)
    print_vtk(np.float32(x), np.float32(y), np.float32(z), gamma, Qp * J, i+1)

  Q.block_until_ready()
  elapsed_time = time.time() - start
  print("elapsed_time:", elapsed_time)

# global variables
gamma = 1.4e0
dt    = 0.01e0
rho0  = 0.4e0 
M0    = 1.e0

main(nx,ny,nz,gamma,R,dt,rho0,M0)

