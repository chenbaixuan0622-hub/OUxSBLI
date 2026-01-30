program main1d
  use, intrinsic :: iso_fortran_env
  use mpi
  use mod_globals, only : nx
  use set
  use calc_time_dev
  implicit none
  real(8) t_start, t_end
  real(8), allocatable :: x(:) 
  real(8), allocatable, pinned :: Q(:,:)
  ! MPI
  integer ierr, nranks, myrank

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  print *, "my rank is", myrank

  allocate(Q(nx,3), x(nx))

  call set_grid(nx, x)
  call set_init(nx, x, Q)

  call cpu_time(t_start)
  call RungeKutta(myrank, nx, x, Q)
  call cpu_time(t_end)

  print *, "elapsed time:", t_end - t_start

  deallocate(Q, x)
  call MPI_FINALIZE(ierr)
end program main1d

