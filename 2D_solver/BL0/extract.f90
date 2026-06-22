program extract
  implicit none
  integer j, nx, ny, nxe, filesize, ios
  real(8) :: tmp(5)
  real(8), allocatable :: x(:), y(:), Q(:,:,:)
  character(len=40) filename

  write(filename, "(a)") "./recal/x.dat"
  open(10, file=filename, action="read", form="unformatted", access="stream", status="old", iostat=ios) 
  if (ios /= 0) then
    print *, "Error opening file x.dat."
  endif
  inquire(unit=10, size=filesize)
  nx = filesize / 8
  allocate(x(nx))
  read(10) x
  close(10)

  write(filename, "(a)") "./recal/y.dat"
  open(10, file=filename, action="read", form="unformatted", access="stream", status="old", iostat=ios) 
  if (ios /= 0) then
    print *, "Error opening file y.dat."
  endif
  inquire(unit=10, size=filesize)
  ny = filesize / 8
  allocate(y(ny))
  read(10) y
  close(10)

  print *, "Grid size: nx=", nx, ", ny=", ny
  allocate(Q(nx,4,ny))
  
  write(filename, "(a)") "./recal/Q00001.dat"
  open(10, file=filename, action="read", form="unformatted", access="stream", status="old", iostat=ios) 
  if (ios /= 0) then
    print *, "Error opening file Q.dat."
  endif
  read(10) Q
  close(10)

  ! extract surface
  print *, "Choose extract line nxe:"
  read(*,*) nxe

  write(filename, "(a, i5.5, a)") "./Qx", int(nxe), ".dat"
  open(10, file=filename, action="write", form="formatted", status="replace") 
  write(10, "(a)") "y rho rhou rhov e"
  do j = 1, ny
    tmp = [y(j), Q(nxe,1,j), Q(nxe,2,j), Q(nxe,3,j), Q(nxe,4,j)]
    write(10, *) tmp
  enddo
  close(10)

  deallocate(x, y, Q)
end program extract
