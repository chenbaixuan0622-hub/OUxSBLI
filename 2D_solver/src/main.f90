program main
  use, intrinsic :: iso_fortran_env
  use mod_globals, only : id_visc, Lx, Ly, nx, ny, nt, np, dx, dy, dt, gamma, T
  use set
  use calc_time_dev
  implicit none
  real(8) mu, kappa, Cv, Cp
  real(8), allocatable :: Q(:,:,:)
  real(8), allocatable :: T0(:,:)
  real(8) t_start, t_end
  
  ! calc physical properties
  mu = (1.4592d-6 * T ** (1.5d0)) / (109.1d0 + T)
  kappa = (2.334d-3 * T ** (1.5d0)) / (164.54d0 + T)
  Cp = 1030.5d0 - 0.19975d0 * T + 3.9734d-4 * T ** 2
  Cv = Cp / gamma
  open(1, file="output.d")
  write(1,"('physical properties')")
  write(1,"('gamma =', f7.4)") gamma
  if (id_visc == 1) then
    write(1,"('Navier-Stokes solver was chosen')")
    write(1,"('mu    =', e12.4, '[Pa s]')") mu
    write(1,"('kappa =', e12.4, '[W/(m K)]')") kappa
    write(1,"('Cv    =', f10.4, '[J/(kg K)]')") Cv
    write(1,"('Cp    =', f10.4, '[J/(kg K)]')") Cp
  else
    write(1,"('Euler solver was chosen')")
  endif
  write(1,"('\n')")

  write(1,"('mesh info')")
  write(1,"('Lx =', f9.4, ' was devided by', i4, ' dx =', e12.4)") Lx, (nx-1), dx
  write(1,"('Ly =', f9.4, ' was devided by', i4, ' dy =', e12.4)") Ly, (ny-1), dy
  write(1,"('\n')")
  close(1)

  ! parallel computing
  allocate(Q(nx,ny,4))
  allocate(T0(nx,ny))
  
  call set_init(Q)
  T0 = T

  call cpu_time(t_start)
  call RungeKutta(T0,Q)
  call cpu_time(t_end)
  
  open(1,file="output.d",position='append')
  write(1,"('time info')")
  write(1,"('dt                       =', e12.4)") dt
  write(1,"('end time                 =', e12.4)") nt * np * dt
  write(1,"('elapsed time             =', i10, '[s]')") int(t_end-t_start)
  close(1)
  deallocate(Q,T0)
end program main

