program main
  use, intrinsic :: iso_fortran_env
  use set_init
  use calc_time_dev
  implicit none
  integer nx, ny, nz, nt, np
  real(8) gamma, T, mu, kappa, Cv, Cp, Cs, U_ref, RHO, L, M, pi, Lx, Ly, Lz, dx, dy, dz, dt, non_dt
  real(8), allocatable :: Q(:,:,:,:)
  real(8) t0, t1
  
  ! read file
  open(28,file='input.d', action='read')
  read(28,*) nx
  read(28,*) ny
  read(28,*) nz
  read(28,*) nt
  read(28,*) np
  read(28,*) gamma
  read(28,*) T
  read(28,*) RHO
  read(28,*) L
  read(28,*) M
  read(28,*) Cs
  close(28)

  ! set grid
  pi = acos(-1.0d0)
  Lx = 2.d0 * pi
  Ly = 2.d0 * pi
  Lz = 2.d0 * pi
  dx = Lx / dble(nx - 1)
  dy = Ly / dble(ny - 1)
  dz = Lz / dble(nz - 1)
  
  ! time
  dt = 0.01d0
  U_ref = M

  ! calc physical properties
  mu = (1.4592d-6 * T ** (1.5d0)) / (109.1d0 + T)
  kappa = (2.334d-3 * T ** (1.5d0)) / (164.54d0 + T)
  Cp = 1030.5d0 - 0.19975d0 * T + 3.9734d-4 * T ** 2
  Cv = Cp / gamma
  open(1, file="output.d")
  write(1,"('physical properties')")
  write(1,"('gamma =', f7.4)") gamma
  write(1,"('mu    =', e12.4, '[Pa s]')") mu
  write(1,"('kappa =', e12.4, '[W/(m K)]')") kappa
  write(1,"('Cv    =', f10.4, '[J/(kg K)]')") Cv
  write(1,"('Cp    =', f10.4, '[J/(kg K)]')") Cp
  write(1,"('Re    =', f10.4)") RHO * U_ref * L / mu
  write(1,"('\n')")

  write(1,"('mesh info')")
  write(1,"('Lx =', f9.4, ' was devided by', i4, ' dx =', e12.4)") Lx, (nx-1), dx
  write(1,"('Ly =', f9.4, ' was devided by', i4, ' dy =', e12.4)") Ly, (ny-1), dy
  write(1,"('Lz =', f9.4, ' was devided by', i4, ' dz =', e12.4)") Lz, (nz-1), dz
  write(1,"('\n')")
  close(1)

  ! parallel computing
  allocate(Q(nx,ny,nz,5))
  
  ! set initial condition
  call TaylorGreen(nx,ny,nz,gamma,RHO,L,M,Q)

  call cpu_time(t0)
  call RungeKutta(nx,ny,nz,nt,np,dx,dy,dz,dt,gamma,mu,kappa,Cv,Cp,Cs,Q)
  call cpu_time(t1)
  non_dt = dt / (L / U_ref)
  open(1,file="output.d",position='append')
  write(1,"('time info')")
  write(1,"('dt                       =', e12.4)") dt
  write(1,"('non-dimentional dt       =', e12.4)") non_dt
  write(1,"('Courant number           =', e12.4)") U_ref * dt / dx
  write(1,"('end time                 =', e12.4)") nt * np * dt
  write(1,"('non-dimentional end time =', e12.4)") nt * np * non_dt  
  write(1,"('elapsed time             =', i10, '[s]')") int(t1-t0)
  close(1)
  deallocate(Q)
end program main
