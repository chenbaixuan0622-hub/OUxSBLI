program main
  use, intrinsic :: iso_fortran_env
  use mod_globals, only : nx, ny, nz, id_recal
  use set
  use set_coordinate
  use calc_time_dev
  implicit none
  integer i, j, l
  real(8) t_start, t_end
  real(8), allocatable :: Q(:,:,:,:), x(:), xix(:), dx(:), y(:), etay(:), dy(:), z(:), zetaz(:), dz(:), Jacobian(:,:,:)

  allocate(Q(nx,ny,nz,5),x(nx),xix(nx-1),dx(nx-1),y(ny),etay(ny-1),dy(ny-1),z(nz),zetaz(nz-1),dz(nz-1),Jacobian(nx,ny,nz))

  call set_grid(nx,ny,nz,x,y,z,dx,dy,dz)
  call set_xix(nx,dx,xix)
  call set_etay(ny,dy,etay)
  call set_zetaz(nz,dz,zetaz)

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
  call set_Jacobian(nx,ny,nz,dx,dy,dz,Jacobian)

  call cpu_time(t_start)
  call RungeKutta(nx,ny,nz,x,dx,xix,y,dy,etay,z,dz,zetaz,Jacobian,Q)
  call cpu_time(t_end)

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

  deallocate(Q,x,xix,dx,y,etay,dy,z,zetaz,dz,Jacobian)
end program main

