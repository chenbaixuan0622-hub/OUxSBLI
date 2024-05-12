program main
  use, intrinsic :: iso_fortran_env
  use mod_globals, only : id_RungeKutta, nx, ny, nz, x, y, z, xix, etay, Jacobian, T, Q, T0
  use set
  use set_coordinate
  use calc_time_dev
  implicit none
  real(8) t_start, t_end
  real(8) dx(nx), dy(ny)
  real(8) Vin(ny,2)

  call set_grid(x,y,z,dx,dy)
  call set_xix(dx,xix)
  call set_etay(dy,etay)
  call set_Jacobian(dx,dy,Jacobian)
  call set_init(x,y,z,Q,Vin)
  T0 = T

  call cpu_time(t_start)
  call RungeKutta(id_RungeKutta,x,dx,xix,y,dy,etay,z,Jacobian,T0,Q,Vin)
  call cpu_time(t_end)
  print *, "elapsed time:", t_end - t_start
end program main

