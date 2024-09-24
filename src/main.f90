program main
  use, intrinsic :: iso_fortran_env
  use mpi
  use nvtx
  use mod_globals, only : id_recal, id_RungeKutta, nx, ny, nz
  use set
  use set_coordinate
  use calc_time_dev
  implicit none
  integer i, j, l
  real(8) t_start, t_end
  real(8), allocatable :: x(:), xix(:), dx(:), y(:), etay(:), dy(:), z(:), zetaz(:), dz(:), Jacobian(:,:,:)
  real(8), allocatable, pinned :: Q(:,:,:,:)
  ! MPI
  integer ierr, nranks, myrank

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  print *, "my rank is", myrank

  allocate(Q(nx,ny,nz,5),x(nx),xix(nx-1),dx(nx-1),y(ny),etay(ny-1),dy(ny-1),z(nz),zetaz(nz-1),dz(nz-1),Jacobian(nx,ny,nz))

  call set_grid(nx,ny,nz,x,y,z,dx,dy,dz)
  call set_xix(nx,dx,xix)
  call set_etay(ny,dy,etay)
  call set_zetaz(nz,dz,zetaz)

  if (myrank == 0) then
    if (kind(id_recal) == 4) then
      open(10,file="recal/Q.dat",action="read",form="unformatted",access="stream")
      read(10) Q
      close(10)
      write(*,*) "simulation restarted"
    elseif (kind(id_recal) == 2) then
      write(*,*) "set initial condition"
      call set_init(nx,ny,nz,x,y,z,Q)
    else
      write(*,*) "wrong paramater was found"
    endif
  endif
  call set_Jacobian(nx,ny,nz,dx,dy,dz,Jacobian)

  call cpu_time(t_start)
  call RungeKutta(id_RungeKutta,myrank,nx,ny,nz,x,dx,xix,y,dy,etay,z,dz,zetaz,Jacobian,Q)
  call cpu_time(t_end)

  if (myrank == 0) then
    ! save data
    do l = 1, nz
      do j = 1, ny
        do i = 1, nx
          Q(i,j,l,:) = Jacobian(i,j,l) * Q(i,j,l,:)
    enddo;enddo;enddo
    open(10,file="recal/Q.dat",status="replace",action="write",form="unformatted",access="stream")
    write(10) Q
    close(10)
    print *, "elapsed time:", t_end - t_start
  endif

  deallocate(Q,x,xix,dx,y,etay,dy,z,zetaz,dz,Jacobian)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_FINALIZE(ierr)
end program main

