program main
  use, intrinsic :: iso_fortran_env
  use mod_globals, only : id_hybrid, Lx, Ly, Lz, T, Q, T0
  use set
  use calc_time_dev
  use calc_time_dev_hybrid, RungeKutta_hybrid => RungeKutta
  implicit none
  real(8) t_start, t_end
  
  call set_init(Q)
  T0 = T

  call cpu_time(t_start)
  if (id_hybrid /= 1) then
    call RungeKutta(T0,Q)
  else
    call RungeKutta_hybrid(T0,Q)
  endif
  call cpu_time(t_end)
end program main

