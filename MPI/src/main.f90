program main
  use, intrinsic :: iso_fortran_env
  use mod_allocate
  use mod_globals, only : id_recal, id_RungeKutta, nx, ny, nz
  use set
  use set_coordinate
  use calc_time_dev
  implicit none
  integer i, j
  real(8) t_start, t_end
  real(8), allocatable :: x(:), xix(:), dx(:), y(:), etay(:), dy(:), z(:), Jacobian(:,:), Q(:,:,:,:)
  allocate(x(nx),xix(nx-1),dx(nx-1),y(ny),etay(ny-1),dy(ny-1),z(nx),Jacobian(nx,ny),Q(nx,ny,nz,5))

  call set_grid(nx,ny,nz,x,y,z,dx,dy)
  call set_xix(nx,dx,xix)
  call set_etay(ny,dy,etay)
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
  call set_Jacobian(nx,ny,dx,dy,Jacobian)

  call cpu_time(t_start)
  call RungeKutta(id_RungeKutta,nx,ny,nz,x,dx,xix,y,dy,etay,z,Jacobian,Q)
  call cpu_time(t_end)

  ! save data
  do j = 1, ny
    do i = 1, nx
      Q(i,j,:,:) = Jacobian(i,j) * Q(i,j,:,:)
  enddo;enddo
  open(10,file="recal/Q.dat",status="replace",action="write",form="unformatted",access="stream")
  write(10) Q
  close(10)
  print *, "elapsed time:", t_end - t_start

  deallocate(Q,x,xix,dx,y,etay,dy,z,Jacobian)
end program main

