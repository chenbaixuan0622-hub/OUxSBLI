program main2
  use, intrinsic :: iso_fortran_env
  use mpi
  use nvtx
  use mod_globals, only : id_recal, nx1, ny1, nz1, nx2, ny2, nz2
  use set
  use set_coordinate
  use calc_time_dev2
  implicit none
  integer i, j, nx, ny, nz
  real(8) t_start, t_end
  real(8), allocatable :: x(:), xix(:), dx(:), y(:), etay(:), dy(:), z(:), Jacobian(:,:)
  real(8), allocatable, pinned :: Q(:,:,:,:)
  ! MPI
  integer ierr, nranks, myrank, status(MPI_STATUS_SIZE)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MPI rank  ! 0 ! Q1  ! calc  !
  !           ! 1 !     ! IO    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !           ! 2 ! Q2  ! calc  !
  !           ! 3 !     ! IO    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  if (nranks /= 4) then
    print *, "num of rank must be 4"
    stop
  endif

  print *, "my rank is", myrank

  if (myrank <= 1) then
    nx = nx1
    ny = ny1
    nz = nz1
    allocate(Q(nx,ny,nz,5),x(nx),xix(nx-1),dx(nx-1),y(ny),etay(ny-1),dy(ny-1),z(nz),Jacobian(nx,ny))
  elseif (myrank >= 2) then
    nx = nx2
    ny = ny2
    nz = nz2
    allocate(Q(nx,ny,nz,5),x(nx),xix(nx-1),dx(nx-1),y(ny),etay(ny-1),dy(ny-1),z(nz),Jacobian(nx,ny))
  endif

  ! set grid information
  if (mod(myrank,2) == 0) then
    call set_grid(myrank,nx,ny,nz,x,y,z,dx,dy)
    xix(:)  = 1.d0 / dx(:)
    etay(:) = 1.d0 / dy(:)
    call set_Jacobian(nx,ny,dx,dy,Jacobian)
    call MPI_SEND(x,        nx,         MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
    call MPI_SEND(dx,       nx-1,       MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
    call MPI_SEND(xix,      nx-1,       MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
    call MPI_SEND(y,        ny,         MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
    call MPI_SEND(dy,       ny-1,       MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
    call MPI_SEND(etay,     ny-1,       MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
    call MPI_SEND(z,        nz,         MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
    call MPI_SEND(Jacobian, nx*ny,      MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
  else
    call MPI_RECV(x,        nx,         MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, status, ierr)
    call MPI_RECV(dx,       nx-1,       MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, status, ierr)
    call MPI_RECV(xix,      nx-1,       MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, status, ierr)
    call MPI_RECV(y,        ny,         MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, status, ierr)
    call MPI_RECV(dy,       ny-1,       MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, status, ierr)
    call MPI_RECV(etay,     ny-1,       MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, status, ierr)
    call MPI_RECV(z,        nz,         MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, status, ierr)
    call MPI_RECV(Jacobian, nx*ny,      MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, status, ierr)
  endif

  if (mod(myrank,2) == 0) then
    if (kind(id_recal) == 4) then
      write(*,*) "simulation restarted"
      if (myrank == 0) then
        open(10,file="recal/Q1.dat",action="read",form="unformatted",access="stream")
      else
        open(10,file="recal/Q2.dat",action="read",form="unformatted",access="stream")
      endif
      read(10) Q
      close(10)
    elseif (kind(id_recal) == 2) then
      write(*,*) "set initial condition"
      call set_init(myrank,nx,ny,nz,x,y,z,Q)
    else
      write(*,*) "wrong paramater was found"
    endif
  endif

  ! share initial condition
  if (mod(myrank,2) == 0) then
    call MPI_SEND(Q, nx*ny*nz*5, MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
  else
    call MPI_RECV(Q, nx*ny*nz*5, MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, status, ierr)
  endif

  call cpu_time(t_start)
  call RungeKutta(myrank,nx,ny,nz,x,dx,xix,y,dy,etay,z,Jacobian,Q)
  call cpu_time(t_end)

  if (mod(myrank,2) == 0) then
    ! save data
    do j = 1, ny
      do i = 1, nx
        Q(i,j,:,:) = Jacobian(i,j) * Q(i,j,:,:)
    enddo;enddo
    if (myrank == 0) then
      open(10,file="recal/Q1.dat",status="replace",action="write",form="unformatted",access="stream")
    else
      open(10,file="recal/Q2.dat",status="replace",action="write",form="unformatted",access="stream")
    endif
    write(10) Q
    close(10)
    print *, "elapsed time:", t_end - t_start
  endif

  deallocate(Q,x,xix,dx,y,etay,dy,z,Jacobian)
  call MPI_FINALIZE(ierr)
end program main2

