module print
  use mod_globals, only : nt, nx, ny, nz, dx, dy, dz, dt, gamma
  implicit none
  interface print_vtk
    module procedure print_vtk_2D, print_vtk_3D
  end interface

contains
  subroutine print_entropy(rho0,p0,rho,p,T)
    real(8), intent(in), dimension(nx,ny,nz) :: rho0, p0, rho, p, T
    real(8) Cp, Cv
    real(8) :: R = 287.d0
    real(8) :: ds = 0.d0
    integer i, j, k
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          Cp = 1030.5d0 - 0.19975d0 * T(i,j,k) + 3.9734 * T(i,j,k)**2
          Cv = Cp - R
          ds = ds + Cv * log(p(i,j,k) / p0(i,j,k)) - Cp * log(rho(i,j,k)/ rho0(i,j,k))
        enddo
      enddo
    enddo
  end subroutine print_entropy

  subroutine print_KE(step,rho,u,v,w)
    integer, intent(in) :: step
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w
    real(8) ke, t
    ke = sum(rho * (u**2 + v**2 + w**2))
    t = nt * step * dt
    open(10,file="data/kinetic_energy.d", position="append")
    write(10,"(2(f9.4,1x))") t, ke
    close(10)
  end subroutine print_KE

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

  subroutine print_vtk_2D(step,Q,T)
    integer, intent(in) :: step
    real(8), intent(in) :: Q(nx,ny,4)
    real(8), intent(in), optional :: T(nx,ny)
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
    write(10,"(f11.4,1x)") ((p(i,j),i=1,nx),j=1,ny)

    write(10,"('SCALARS T float')")
    write(10,"('LOOKUP_TABLE default')")
    write(10,"(f9.4,1x)") ((T(i,j),i=1,nx),j=1,ny)
    close(10)
  end subroutine print_vtk_2D
  
  subroutine print_vtk_3D(step,Q,T)
    integer, intent(in) :: step
    real(8), intent(in) :: Q(nx,ny,nz,5)
    real(8), intent(in), optional :: T(nx,ny,nz)
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
    write(10,"(f11.4,1x)") (((p(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    
    write(10,"('SCALARS T float')")
    write(10,"('LOOKUP_TABLE default')")
    write(10,"(f9.4,1x)") (((T(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    close(10)

    !call print_entropy(step,rho0,p0,rho,p,T)
    call print_KE(step,rho,u,v,w)
  end subroutine print_vtk_3D
end module print

