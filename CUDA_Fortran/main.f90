program main
  use, intrinsic :: iso_fortran_env
  use set_init
  use calc_time_dev
  implicit none
  integer nx, ny, nz, nt, np
  real(8) gamma, T, mu, kappa, Cp, RHO, L, M, pi, dx, dy, dz, dt
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
  close(28)

  ! set grid
  pi = acos(-1.0d0)
  dx = 2.0d0 * pi / dble(nx - 1)
  dy = 2.0d0 * pi / dble(ny - 1)
  dz = 2.0d0 * pi / dble(nz - 1)
  
  ! time
  dt = 0.01d0

  ! calc physical properties
  mu = (1.4592d-6 * T ** (1.5d0)) / (109.1d0 + T)
  kappa = (2.334d-3 * T ** (1.5d0)) / (164.54d0 + T)
  Cp = 1030.5d0 - 0.19975d0 * T + 3.9734d-4 * T ** 2
  write(*,*) 'mu =',mu,' [Pa s]'
  write(*,*) 'kappa =',kappa,' [W/(m K)]'
  write(*,*) 'Cp =',Cp,' [J/(kg K)]'

  ! parallel computing
  allocate(Q(nx,ny,nz,5))
  
  ! set initial condition
  call TaylorGreen(nx,ny,nz,gamma,RHO,L,M,Q)

  call cpu_time(t0)
  call RungeKutta(nx,ny,nz,nt,np,dx,dy,dz,dt,gamma,mu,kappa,Cp,Q)
  call cpu_time(t1)
  print *, 'elapsed time:', t1-t0

  deallocate(Q)
end program main