module print
  use mod_globals, only : nt, nx, ny, nz, dx, dy, dz, dt, gamma
  implicit none
  interface print_vtk
    subroutine print_vtk_2D(step,x,y,Jacobian,Q,T)
      integer, intent(in) :: step
      real(8), intent(in), dimension(nx,ny) :: x, y, Jacobian
      real(8), intent(in) :: Q(nx,ny,4)
      real(8), intent(in), optional :: T(nx,ny)
    end subroutine print_vtk_2D

    subroutine print_vtk_3D(step,x,y,z,Jacobian,Q,T,ke0,entropy0,mut)
      integer, intent(in) :: step
      real(8), intent(in), dimension(nx,ny,nz) :: x, y, z, Jacobian
      real(8), intent(in) :: Q(nx,ny,nz,5)
      real(8), intent(in) :: T(nx,ny,nz)
      real(8), intent(inout) :: ke0, entropy0
      real(8), intent(in), optional :: mut(nx,ny,nz)
    end subroutine print_vtk_3D
  end interface
contains
  subroutine print_entropy(step,rho,p,rhos0)
    integer, intent(in) :: step
    real(8), intent(in), dimension(nx,ny,nz) :: rho, p
    real(8), intent(inout) :: rhos0
    real(8) rhos, t
    rhos = sum(rho * log(p * rho ** (-gamma)))
    if (step == 0) then
      rhos0 = rhos
    endif
    t = nt * step * dt
    open(10,file="data/entropy.d", position="append")
    write(10,"(2(f9.4,1x))") t, (rhos0 - rhos) / rhos0
    close(10)
  end subroutine print_entropy

  subroutine print_KE(step,rho,u,v,w,ke0)
    integer, intent(in) :: step
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w
    real(8), intent(inout) :: ke0
    real(8) ke, t
    ke = sum(0.5d0 * rho * (u**2 + v**2 + w**2))
    if (step == 0) then
      ke0 = ke
    endif
    t = nt * step * dt
    open(10,file="data/kinetic_energy.d", position="append")
    write(10,"(2(f9.4,1x))") t, ke / ke0
    close(10)
  end subroutine print_KE

  subroutine print_boundary_layer(u)
    real(8), intent(in) :: u(ny)
    integer j
    open(10,file="data/boundary_layer.d",action="write")
    do j = 2, ny
      write(10,"(2(f12.7,1x))") dble(j-2)*dy, u(j)/u(ny)
    enddo
    close(10)
  end subroutine print_boundary_layer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_header_2D(ni,nj,x,y)
    integer, intent(in) :: ni, nj
    real(8), intent(in), dimension(nx,ny) :: x, y
    integer i, j
    write(10,"('# vtk DataFile Version 3.0')")
    write(10,"('Q')")
    write(10,"('ASCII')")
    write(10,"('DATASET STRUCTURED_GRID')")
    write(10,"('DIMENSIONS',3(1x,i4))") ni, nj, 1
    write(10,"('POINTS',i9,' float')") ni * nj
    write(10,"(3(f12.7,1x))") ((x(i,j), y(i,j), 0.d0,i=1,ni),j=1,nj)

    write(10,"('POINT_DATA',i9)") ni * nj
    write(10,"('VECTORS Velocity float')")
  end subroutine print_header_2D

  subroutine print_header_3D(ni,nj,nk,x,y,z)
    integer, intent(in) :: ni, nj, nk
    real(8), intent(in), dimension(nx,ny,nz) :: x, y, z
    integer i, j, k
    write(10,"('# vtk DataFile Version 3.0')")
    write(10,"('Q')")
    write(10,"('ASCII')")
    write(10,"('DATASET STRUCTURED_GRID')")
    write(10,"('DIMENSIONS',3(1x,i4))") ni, nj, nk
    write(10,"('POINTS',i9,' float')") ni * nj * nk
    write(10,"(3(f12.7,1x))") (((x(i,j,k), y(i,j,k), z(i,j,k),i=1,ni),j=1,nj),k=1,nk)

    write(10,"('POINT_DATA',i9)") ni * nj * nk
    write(10,"('VECTORS Velocity float')")
  end subroutine print_header_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_vtk_2D(step,x,y,Jacobian,Q,T)
    integer, intent(in) :: step
    real(8), intent(in), dimension(nx,ny) :: x, y, Jacobian
    real(8), intent(inout) :: Q(nx,ny,4)
    real(8), intent(in), optional :: T(nx,ny)
    integer i, j
    real(8), dimension(nx,ny) :: rho, u, v, p
    character(len=40) filename
    do i = 1, 4
      Q(:,:,i) = Q(:,:,i) * Jacobian(:,:)
    enddo
    rho = Q(:,:,1)
    u = Q(:,:,2) / rho
    v = Q(:,:,3) / rho
    p = (gamma - 1.d0) * (Q(:,:,4) - 0.5d0 * rho * (u**2 + v**2))
    
    write(filename, "(a, i5.5,a)") "data/Q",int(step),".vtk"
    open(10,file=filename)
    call print_header_2D(nx,ny,x,y)

    write(10,"(3(f10.4,1x))") ((u(i,j), v(i,j), 0.d0,i=1,nx),j=1,ny)

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
  
    call print_boundary_layer(u(int(0.5*nx),:))
  end subroutine print_vtk_2D
  
  subroutine print_vtk_3D(step,x,y,z,Jacobian,Q,T,ke0,rhos0,mut)
    integer, intent(in) :: step
    real(8), intent(in), dimension(nx,ny,nz) :: x, y, z, Jacobian
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    real(8), intent(in) :: T(nx,ny,nz)
    real(8), intent(inout) :: ke0, rhos0
    real(8), intent(in), optional :: mut(nx,ny,nz)
    integer i, j, k
    real(8), dimension(nx,ny,nz) :: rho, u, v, w, p, nut
    character(len=40) filename
    do i = 1, 5
      Q(:,:,:,i) = Q(:,:,:,i) * Jacobian(:,:,:)
    enddo
    rho = Q(:,:,:,1)
    u = Q(:,:,:,2) / rho
    v = Q(:,:,:,3) / rho
    w = Q(:,:,:,4) / rho
    p = (gamma - 1.d0) * (Q(:,:,:,5) - 0.5d0 * rho * (u**2 + v**2 + w**2))
    
    write(filename, "(a, i5.5,a)") "data/Q",int(step),".vtk"
    open(10,file=filename)
    call print_header_3D(nx,ny,nz,x,y,z)

    write(10,"(3(f10.4,1x))") (((u(i,j,k), v(i,j,k), w(i,j,k),i=1,nx),j=1,ny),k=1,nz)

    write(10,"('SCALARS rho float')")
    write(10,"('LOOKUP_TABLE default')")
    write(10,"(f9.4,1x)") (((rho(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    
    write(10,"('SCALARS P float')")
    write(10,"('LOOKUP_TABLE default')")
    write(10,"(f11.4,1x)") (((p(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    
    write(10,"('SCALARS T float')")
    write(10,"('LOOKUP_TABLE default')")
    write(10,"(f9.4,1x)") (((T(i,j,k),i=1,nx),j=1,ny),k=1,nz)

    if (present(mut)) then
      nut(:,:,:) = mut(:,:,:) / rho(:,:,:)
      write(10,"('SCALARS nut float')")
      write(10,"('LOOKUP_TABLE default')")
      write(10,"(f12.10,1x)") (((nut(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    endif
    close(10)

    call print_entropy(step,rho,p,rhos0)
    call print_KE(step,rho,u,v,w,ke0)
  end subroutine print_vtk_3D
end module print

