module print
  use mod_globals, only : gamma
  implicit none
contains
  subroutine print_vtk(step,nx,ny,dx,dy,Q)
    integer, intent(in) :: step, nx, ny
    real(8), intent(in) :: dx, dy, Q(nx,ny,4)
    integer i, j
    real(8), dimension(nx,ny) :: rho, u, v, p
    character(len=40) filename
    rho = Q(:,:,1)
    u = Q(:,:,2) / rho
    v = Q(:,:,3) / rho
    p = (gamma - 1.d0) * (Q(:,:,4) - 0.5d0 * rho * (u ** 2 + v ** 2 ))
    
    write(filename, "(a, i5.5,a)") "data/Q",int(step),".vtk"
    open(10,file=filename)
    write(10,"('# vtk DataFile Version 3.0')")
    write(10,"('Q')")
    write(10,"('ASCII')")
    write(10,"('DATASET STRUCTURED_GRID')")
    write(10,"('DIMENSIONS',2(1x,i4))") nx, ny
    write(10,"('POINTS',i9,' float')") nx * ny
    do j = 1, ny
      do i = 1, nx
        write(10,"(2(f9.4,1x))") (i-1)*dx, (j-1)*dy
      enddo
    enddo 

    write(10,"('POINT_DATA',i9)") nx * ny
    write(10,"('VECTORS Velocity float')")
    do j = 1, ny
      do i = 1, nx
        write(10,"(2(f9.4,1x))") u(i,j), v(i,j)
      enddo
    enddo

    write(10,"('SCALARS rho float')")
    write(10,"('LOOKUP_TABLE default')")
    do j = 1, ny
      do i = 1, nx
        write(10,"(f9.4,1x)") rho(i,j)
      enddo
    enddo
    
    write(10,"('SCALARS P float')")
    write(10,"('LOOKUP_TABLE default')")
    do j = 1, ny
      do i = 1, nx
        write(10,"(f9.4,1x)") p(i,j)
      enddo
    enddo
    close(10)
  end subroutine print_vtk
end module print

