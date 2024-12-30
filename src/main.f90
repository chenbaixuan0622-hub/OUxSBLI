program main
  use, intrinsic :: iso_fortran_env
  use mpi
  use mod_globals, only : id_RungeKutta, id_recal, nx, ny, nz, Lx, Ly, Lz
  use set
  use set_coordinate
  use calc_time_dev
  implicit none
  integer i, j, l, m, s
  real(8) t_start, t_end
  real(8), allocatable :: x(:), dx(:), y(:), dy(:), z(:), dz(:), Jacobian(:), Q(:,:,:,:)
  character(len=40) filename
  ! MPI
  integer nranks, myrank, ierr, ireq, istat(MPI_STATUS_SIZE)

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

  print *, "my rank is", myrank

  allocate(Q(nx,ny,nz,5),x(nx),dx(nx),y(ny),dy(ny),z(nz),dz(nz),Jacobian(ny))

  ! set grid information
  if (mod(myrank,2) == 0) then
    call set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    call set_Jacobian_y(nx, ny, nz, dx, dy, dz, Jacobian)
    if (kind(id_recal) == 4) then
      write(filename, "(a, i5.5, a)") "recal/Q", int(myrank/2+1), ".dat"
      write(*,*) "simulation restarted"
      open(10,file=filename,action="read",form="unformatted",access="stream")
      read(10) Q
      close(10)
    elseif (kind(id_recal) == 2) then
      write(*,*) "set initial condition"
      call set_init(myrank,nx,ny,nz,x,y,z,Q)
    else
      write(*,*) "wrong paramater was found"
    endif
  else
    call set_grid(myrank-1, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    call set_Jacobian_y(nx, ny, nz, dx, dy, dz, Jacobian)
  endif

  call cpu_time(t_start)
  call RungeKutta(id_RungeKutta, myrank, nx, ny, nz, x, dx, y, dy, z, dz, Jacobian, Q)
  call cpu_time(t_end)

  if (mod(myrank,2) == 0) then
    ! save data
    do l = 1, nz
      do j = 1, ny
        do i = 1, nx
          Q(i,j,l,:) = Jacobian(j) * Q(i,j,l,:)
    enddo;enddo;enddo
    write(filename, "(a, i5.5, a)") "recal/Q", int(myrank/2+1), ".dat"
    open(10,file=filename,status="replace",action="write",form="unformatted",access="stream")
    write(10) Q
    close(10)
    if (t_end - t_start <= 60.d0) then
      s = int(t_end - t_start)
      print *, "elapsed time:", s, " [sec]"
    else
      m = int(t_end - t_start) / 60
      s = int(t_end - t_start) - 60 * m
      print *, "elapsed time:", m, " [min] ", s, " [sec]"
    endif
  endif

  deallocate(Q,x,dx,y,dy,z,dz,Jacobian)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_FINALIZE(ierr)
end program main

