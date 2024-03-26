program main
  use, intrinsic :: iso_fortran_env
  use mod_globals, only : Lx, Ly, Lz, T, Q, T0
  use set
  use calc_time_dev
  implicit none
  real(8) t_start, t_end
  
  call set_init(Q)
  T0 = T

  call cpu_time(t_start)
  call RungeKutta(T0,Q)
  call cpu_time(t_end)
  print *, "elapsed time:", t_end - t_start
end program main

