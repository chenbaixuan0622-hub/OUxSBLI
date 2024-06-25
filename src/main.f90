module mod_allocate
  implicit none
  interface alloc
    module procedure alloc2D, alloc3D
  end interface
contains
  subroutine alloc2D(Q,z)
    use mod_globals, only : nx, ny
    real(8), intent(out), allocatable :: Q(:,:,:)
    real(8), intent(out), allocatable :: z(:)
    allocate(Q(nx,ny,4),z(1))
  end subroutine alloc2D

  subroutine alloc3D(Q,z)
    use mod_globals, only : nx, ny, nz
    real(8), intent(out), allocatable :: Q(:,:,:,:)
    real(8), intent(out), allocatable :: z(:)
    allocate(Q(nx,ny,nz,5),z(nz))
  end subroutine alloc3D
end module mod_allocate

program main
  use, intrinsic :: iso_fortran_env
  use mpi
  use nvtx
  use mod_allocate
  use mod_globals, only : id_recal, id_RungeKutta, nx, ny, nz, Q
  use set
  use set_coordinate
  use calc_time_dev
  implicit none
  integer i, j
  real(8) t_start, t_end
  real(8), allocatable :: x(:), xix(:), dx(:), y(:), etay(:), dy(:), z(:), Jacobian(:,:)
  ! MPI
  integer ierr, nranks, myrank

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  print *, "my rank is", myrank

  allocate(x(nx),xix(nx-1),dx(nx-1),y(ny),etay(ny-1),dy(ny-1),Jacobian(nx,ny))

  call alloc(Q,z)

  call set_grid(nx,ny,nz,x,y,z,dx,dy)
  call set_xix(nx,dx,xix)
  call set_etay(ny,dy,etay)

  if (myrank == 0) then
    if (kind(id_recal) == 4) then
      write(*,*) "simulation restarted"
      open(10,file="recal/Q.dat",action="read",form="unformatted",access="stream")
      read(10) Q
      close(10)
    elseif (kind(id_recal) == 2) then
      write(*,*) "set initial condition"
      call set_init(nx,ny,nz,x,y,z,Q)
    else
      write(*,*) "wrong paramater was found"
    endif
  endif
  call set_Jacobian(nx,ny,dx,dy,Jacobian)

  call cpu_time(t_start)
  call RungeKutta(id_RungeKutta,myrank,nx,ny,nz,x,dx,xix,y,dy,etay,z,Jacobian,Q)
  call cpu_time(t_end)

  if (myrank == 0) then
    ! save data
    do j = 1, ny
      do i = 1, nx
        Q(i,j,:,:) = Jacobian(i,j) * Q(i,j,:,:)
    enddo;enddo
    open(10,file="recal/Q.dat",status="replace",action="write",form="unformatted",access="stream")
    write(10) Q
    close(10)
    print *, "elapsed time:", t_end - t_start
  endif

  deallocate(Q,x,xix,dx,y,etay,dy,z,Jacobian)
  call MPI_FINALIZE(ierr)
end program main

