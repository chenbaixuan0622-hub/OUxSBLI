program main
  use, intrinsic :: iso_fortran_env
  use mod_globals, only : id_RungeKutta, nx, ny, nz, x, y, z, xix, etay, Jacobian, T, Q, T0
  use set
  use set_coordinate
  use calc_time_dev
  implicit none
  real(8) t_start, t_end
  real(8) Vin(ny,2)

  call set_grid(x,y,z)
  call set_xix(x,xix)
  call set_etay(y,etay)
  call set_Jacobian(x,y,xix,etay,Jacobian)
  call set_init(Q,Vin)
  T0 = T

  call cpu_time(t_start)
  call RungeKutta(id_RungeKutta,x,y,z,xix,etay,Jacobian,T0,Q,Vin)
  call cpu_time(t_end)
  print *, "elapsed time:", t_end - t_start
end program main

