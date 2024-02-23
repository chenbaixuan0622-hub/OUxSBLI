program main
  use, intrinsic :: iso_fortran_env
  use set_parallel
  use set_init
  use calc_time_dev
  implicit none
  integer nx, ny, nz, nyp, nzp, nt, np, ni, nj, nk
  integer blockDim(3)
  real(8) gamma, T, mu, kappa, Cp, RHO, L, M, pi, dx, dy, dz, dt
  real(8), allocatable :: Q(:,:,:,:,:)
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
  
  if (mod(ny,2) /= 0 .or. mod(nz,2) /= 0) then
    write(*,*) 'ny and nz should be even number'
    stop
  endif

  ! set grid
  pi = acos(-1.0d0)
  dx = 2.0d0 * pi / dble(nx - 1)
  dy = 2.0d0 * pi / dble(ny - 1)
  dz = 2.0d0 * pi / dble(nz - 1)
  
  dx = 1.d0 / dx
  dy = 1.d0 / dy
  dz = 1.d0 / dz

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
  call set_index(ny,nz,nyp,nzp)
  allocate(Q(nx,nyp,nzp,5,4))

  ! define block size
  ! 8 * 8 * 8 = 512 threads
  blockDim(1) = 8
  blockDim(2) = 8
  blockDim(3) = 8
  ni = (nx - 4) / blockDim(1)
  nj = (nyp - 4) / blockDim(2)
  nk = (nzp - 4) / blockDim(3)
  if (mod(nx-4,ni)/=0 .or. mod(nyp-4,nj)/=0 .or. mod(nzp-4,nk)/=0 .or. ni<5 .or. nj<5 .or. nk<5) then
    write(*,*) 'Cannot split properly'
    stop
  endif
  
  ! set initial condition
  call TaylorGreen(nx,ny,nz,nyp,nzp,gamma,RHO,L,M,Q)

  call cpu_time(t0)
  call RungeKutta(nx,ny,nz,nyp,nzp,ni,nj,nk,blockDim,nt,np,dx,dy,dz,dt,gamma,mu,kappa,Cp,Q)
  call cpu_time(t1)
  print *, 'elapsed time:', t1-t0

  deallocate(Q)
end program main