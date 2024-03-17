module print
  use mod_globals, only : nx, ny, nz, dx, dy, dz, gamma
  implicit none
  interface print_vtk
    module procedure print_vtk_2D, print_vtk_3D
  end interface

contains
  subroutine print_header(ni,nj,nk,di,dj,dk)
    integer, intent(in) :: ni, nj, nk
    real(8), intent(in) :: di, dj, dk
    integer i, j, k
    write(10,"('# vtk DataFile Version 3.0')")
    write(10,"('Q')")
    write(10,"('ASCII')")
    write(10,"('DATASET STRUCTURED_GRID')")
    write(10,"('DIMENSIONS',3(1x,i4))") ni, nj, nk
    write(10,"('POINTS',i9,' float')") ni * nj * nk
    write(10,"(3(f9.4,1x))") ((((i-1)*di, (j-1)*dj, (k-1)*dk,i=1,ni),j=1,nj),k=1,nk)

    write(10,"('POINT_DATA',i9)") ni * nj * nk
    write(10,"('VECTORS Velocity float')")
  end subroutine print_header

  subroutine print_vtk_2D(step,Q)
    integer, intent(in) :: step
    real(8), intent(in) :: Q(nx,ny,4)
    integer i, j
    real(8), dimension(nx,ny) :: rho, u, v, p
    character(len=40) filename
    rho = Q(:,:,1)
    u = Q(:,:,2) / rho
    v = Q(:,:,3) / rho
    p = (gamma - 1.d0) * (Q(:,:,4) - 0.5d0 * rho * (u**2 + v**2))
    
    write(filename, "(a, i5.5,a)") "data/Q",int(step),".vtk"
    open(10,file=filename)
    call print_header(nx,ny,1,dx,dy,0.d0)

    write(10,"(3(f9.4,1x))") ((u(i,j), v(i,j), 0.d0,i=1,nx),j=1,ny)

    write(10,"('SCALARS rho float')")
    write(10,"('LOOKUP_TABLE default')")
    write(10,"(f9.4,1x)") ((rho(i,j),i=1,nx),j=1,ny)
    
    write(10,"('SCALARS P float')")
    write(10,"('LOOKUP_TABLE default')")
    write(10,"(f9.4,1x)") ((p(i,j),i=1,nx),j=1,ny)
    close(10)
  end subroutine print_vtk_2D
  
  subroutine print_vtk_3D(step,Q)
    integer, intent(in) :: step
    real(8), intent(in) :: Q(nx,ny,nz,5)
    integer i, j, k
    real(8), dimension(nx,ny,nz) :: rho, u, v, w, p
    character(len=40) filename
    rho = Q(:,:,:,1)
    u = Q(:,:,:,2) / rho
    v = Q(:,:,:,3) / rho
    w = Q(:,:,:,4) / rho
    p = (gamma - 1.d0) * (Q(:,:,:,5) - 0.5d0 * rho * (u**2 + v**2 + w**2))
    
    write(filename, "(a, i5.5,a)") "data/Q",int(step),".vtk"
    open(10,file=filename)
    call print_header(nx,ny,nz,dx,dy,dz)

    write(10,"(3(f9.4,1x))") (((u(i,j,k), v(i,j,k), w(i,j,k),i=1,nx),j=1,ny),k=1,nz)

    write(10,"('SCALARS rho float')")
    write(10,"('LOOKUP_TABLE default')")
    write(10,"(f9.4,1x)") (((rho(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    
    write(10,"('SCALARS P float')")
    write(10,"('LOOKUP_TABLE default')")
    write(10,"(f9.4,1x)") (((p(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    close(10)
  end subroutine print_vtk_3D
end module print

