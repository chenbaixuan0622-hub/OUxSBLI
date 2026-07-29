program main
  use, intrinsic :: iso_fortran_env
  use cudafor
  use mpi
  use mod_globals, only : nx, threadsE, threadsEv, threads, blocksE, blocksEv, blocks
  use set
  use calc_time_dev
  use print_1d
  implicit none
  integer mygpu, s, m
  real(8) t_start, t_end
  real(8), allocatable :: x(:), Q(:,:)
  ! MPI (single-rank; kept for structural consistency with the 2D/3D solvers)
  integer nranks, myrank, ierr

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  mygpu = 0

  blocksE  = dim3((nx-1+threadsE%x -1)/threadsE%x,  1, 1)
  blocksEv = dim3((nx-1+threadsEv%x-1)/threadsEv%x, 1, 1)
  blocks   = dim3((nx-2+threads%x  -1)/threads%x,   1, 1)

  allocate(Q(nx,3), x(nx))
  call set_grid(myrank, nx, x)
  call set_init(myrank, nx, x, Q)

  call cpu_time(t_start)
  call RungeKutta(myrank, mygpu, nx, Q)
  call cpu_time(t_end)

  if (t_end - t_start <= 60.d0) then
    s = int(t_end - t_start)
    print *, "calculation time:", s, " [sec]"
  else
    m = int(t_end - t_start) / 60
    s = int(t_end - t_start) - 60 * m
    print *, "calculation time:", m, " [min] ", s, " [sec]"
  endif

  call write_Q_dat(nx, x, Q)

  deallocate(Q, x)
  call MPI_FINALIZE(ierr)
end program main
