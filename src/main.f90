program main
  use, intrinsic :: iso_fortran_env
  use mod_globals, only : id_RungeKutta, nx, ny, nz, x, y, z, xix, etay, Jacobian, T, Q, T0
  use set
  use set_coordinate
  use calc_time_dev
  implicit none
  real(8) t_start, t_end
  real(8), allocatable :: dxs(:), dys(:)
  real(8) Vin(ny,2)

  allocate(dxs(nx),dys(ny))

  call set_grid(x,y,z,dxs,dys)
  call set_xix(dxs,xix)
  call set_etay(dys,etay)
  call set_Jacobian(dxs,dys,Jacobian)
  call set_init(Q,Vin)
  T0 = T

  deallocate(dxs,dys)

  call cpu_time(t_start)
  call RungeKutta(id_RungeKutta,x,y,z,xix,etay,Jacobian,T0,Q,Vin)
  call cpu_time(t_end)
  print *, "elapsed time:", t_end - t_start
end program main

