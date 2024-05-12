module print
  use mod_globals, only : nt, nx, ny, nz, dt, gamma
  implicit none
  interface print_vtk
    subroutine print_vtk_2D(step,x,y,Jacobian,Q,T)
      integer, intent(in)           :: step
      real(8), intent(in)           :: x(nx), y(ny), Jacobian(nx,ny)
      real(8), intent(in)           :: Q(nx,ny,4)
      real(8), intent(in), optional :: T(nx,ny)
    end subroutine print_vtk_2D

    subroutine print_vtk_3D(step,x,y,z,Jacobian,Q,T,ke0,entropy0,mut)
      integer, intent(in)           :: step
      real(8), intent(in)           :: x(nx), y(ny), z(nz), Jacobian(nx,ny)
      real(8), intent(in)           :: Q(nx,ny,nz,5)
      real(8), intent(in)           :: T(nx,ny,nz)
      real(8), intent(inout)        :: ke0, entropy0
      real(8), intent(in), optional :: mut(nx,ny,nz)
    end subroutine print_vtk_3D
  end interface
contains
  function mean(a) result(ans)
    real(8), intent(in) :: a(:,:,:)
    real(8) ans
    ans = sum(a) / size(a)
  end function mean

  subroutine calc_vorticity(x,y,z,u,v,w,omegax,omegay,omegaz)
    real(8), intent(in)                             :: x(nx), y(ny), z(nz)
    real(8), intent(in), dimension(nx,ny,nz)        :: u, v, w
    real(8), intent(out), dimension(nx-2,ny-2,nz-2) :: omegax, omegay, omegaz
    integer i, j, k
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          omegax(i-1,j-1,k-1) = (-w(i,j-1,k) + w(i,j+1,k)) / (-y(j-1) + y(j+1)) &
          & - (-v(i,j,k-1) + v(i,j,k+1)) / (-z(k-1) + z(k+1))
          omegay(i-1,j-1,k-1) = (-u(i,j,k-1) + u(i,j,k+1)) / (-z(k-1) + z(k+1)) &
          & - (-w(i-1,j,k) + w(i+1,j,k)) / (-x(i-1) + x(i+1))
          omegaz(i-1,j-1,k-1) = (-v(i-1,j,k) + v(i+1,j,k)) / (-x(i-1) + x(i+1)) &
          & - (-u(i,j-1,k) + u(i,j+1,k)) / (-y(j-1) + y(j+1))
    enddo;enddo;enddo
  end subroutine calc_vorticity

  subroutine print_entropy(step,rho,p,entropy0)
    integer, intent(in)                       :: step
    real(8), intent(in), dimension(nx,ny,nz)  :: rho, p
    real(8), intent(inout)                    :: entropy0
    real(8) entropy, t
    entropy = sum(rho * log(p * rho ** (-gamma)))
    if (step == 0) then
      entropy0 = entropy
    endif
    t = nt * step * dt
    open(10,file="data/entropy.d", position="append")
    write(10,"(2(f9.4,1x))") t, (entropy0 - entropy) / entropy0
    close(10)
  end subroutine print_entropy

  subroutine print_KE(step,rho,u,v,w,ke0)
    integer, intent(in)                       :: step
    real(8), intent(in), dimension(nx,ny,nz)  :: rho, u, v, w
    real(8), intent(inout)                    :: ke0
    real(8) ke, t
    ke = mean(0.5d0 * rho * (u**2 + v**2 + w**2))
    if (step == 0) then
      ke0 = ke
    endif
    t = nt * step * dt
    open(10,file="data/kinetic_energy.d", position="append")
    !write(10,"(2(f9.4,1x))") t, ke
    write(10,"(2e12.4)") t, ke
    close(10)
  end subroutine print_KE

  subroutine print_enstrophy(step,x,y,z,rho,u,v,w)
    integer, intent(in)                       :: step
    real(8), intent(in)                       :: x(nx), y(ny), z(nz)
    real(8), intent(in), dimension(nx,ny,nz)  :: rho, u, v, w
    real(8) enstrophy, t
    real(8), dimension(nx-2,ny-2,nz-2) :: omegax, omegay, omegaz
    t = nt * step * dt
    call calc_vorticity(x,y,z,u,v,w,omegax,omegay,omegaz)
    enstrophy = mean(0.5d0 * rho(2:nx-1,2:ny-1,2:nz-1) * (omegax**2 + omegay**2 + omegaz**2))
    open(10,file="data/enstrophy.d", position="append")
    !write(10,"(2(f9.4,1x))") t, enstrophy
    write(10,"(2e12.4)") t, enstrophy
    close(10)
  end subroutine print_enstrophy

  subroutine print_boundary_layer(y,u)
    real(8), intent(in), dimension(ny) :: y, u
    integer j
    open(10,file="data/boundary_layer.d",action="write")
    do j = 1, ny
      write(10,"(2(f12.7,1x))") y(j), u(j)/u(ny)
    enddo
    close(10)
  end subroutine print_boundary_layer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_header(ni,nj,nk,x,y,z)
    integer, intent(in) :: ni, nj, nk
    real(8), intent(in) :: x(ni), y(nj), z(nk)
    integer i, j, k
    character :: lf*1, str1*8, str2*8, str3*8, str4*8
    lf = char(10)
    write(str1(1:8),'(i8)') ni
    write(str2(1:8),'(i8)') nj
    write(str3(1:8),'(i8)') nk
    write(str4(1:8),'(i8)') ni * nj * nk

    write(10) '# vtk DataFile Version 3.0'//lf
    write(10) 'Q'//lf
    write(10) 'BINARY'//lf
    write(10) 'DATASET RECTILINEAR_GRID'//lf
    write(10) 'DIMENSIONS'//str1//str2//str3//lf
    write(10) 'X_COORDINATES'//str1//'  float'//lf
    write(10) real(x), lf
    write(10) lf//'Y_COORDINATES'//str2//'  float'//lf
    write(10) real(y), lf
    write(10) lf//'Z_COORDINATES'//str3//'  float'//lf
    write(10) real(z), lf

    write(10) lf//'POINT_DATA'//str4//lf
  end subroutine print_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_vtk_2D(step,x,y,Jacobian,Q,T)
    integer, intent(in)                   :: step
    real(8), intent(in)                   :: x(nx), y(ny)
    real(8), intent(in), dimension(nx,ny) :: Jacobian
    real(8), intent(in)                   :: Q(nx,ny,4)
    real(8), intent(in), optional         :: T(nx,ny)
    integer i, j
    real(8), dimension(nx,ny) :: rho, u, v, p
    real(8) :: z(1) = 0.d0
    character(len=40) filename
    character :: lf*1
    lf = char(10)
    do j = 1, ny
      do i = 1, nx
        rho(i,j) = Jacobian(i,j) * Q(i,j,1)
        u(i,j) = Jacobian(i,j) * Q(i,j,2) / rho(i,j)
        v(i,j) = Jacobian(i,j) * Q(i,j,3) / rho(i,j)
        p(i,j) = (gamma - 1.d0) * (Jacobian(i,j) * Q(i,j,4) - 0.5d0 * rho(i,j) * (u(i,j)**2 + v(i,j)**2))
    enddo;enddo

    write(filename, "(a, i5.5,a)") "data/Q",int(step),".vtk"
    open(10,file=filename,status="replace",action="write",form="unformatted",access="stream",convert="BIG_ENDIAN")
    call print_header(nx,ny,1,x,y,z)

    write(10) 'VECTORS Velocity float'//lf
    do j = 1, ny
      do i = 1, nx
        write(10) real(u(i,j)), real(v(i,j)), 0.e0
    enddo;enddo

    write(10) lf//'SCALARS rho float'//lf
    write(10) 'LOOKUP_TABLE default'//lf
    write(10) real(rho), lf
    
    write(10) lf//'SCALARS P float'//lf
    write(10) 'LOOKUP_TABLE default'//lf
    write(10) real(p), lf

    write(10) lf//'SCALARS T float'//lf
    write(10) 'LOOKUP_TABLE default'//lf
    write(10) real(T), lf
    close(10)
  
    call print_boundary_layer(y,u(int(0.5*nx),:))
  end subroutine print_vtk_2D
  
  subroutine print_vtk_3D(step,x,y,z,Jacobian,Q,T,ke0,entropy0,mut)
    integer, intent(in)                   :: step
    real(8), intent(in)                   :: x(nx), y(ny), z(nz)
    real(8), intent(in), dimension(nx,ny) :: Jacobian
    real(8), intent(in)                   :: Q(nx,ny,nz,5)
    real(8), intent(in)                   :: T(nx,ny,nz)
    real(8), intent(inout)                :: ke0, entropy0
    real(8), intent(in), optional         :: mut(nx,ny,nz)
    integer i, j, k, l
    real(8), dimension(nx,ny,nz) :: rho, u, v, w, p, nut
    character(len=40) filename
    character :: lf*1
    lf = char(10)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rho(i,j,k) = Jacobian(i,j) * Q(i,j,k,1)
          u(i,j,k) = Jacobian(i,j) * Q(i,j,k,2) / rho(i,j,k)
          v(i,j,k) = Jacobian(i,j) * Q(i,j,k,3) / rho(i,j,k)
          w(i,j,k) = Jacobian(i,j) * Q(i,j,k,4) / rho(i,j,k)
          p(i,j,k) = (gamma - 1.d0) * (Jacobian(i,j) * Q(i,j,k,5) - 0.5d0 * rho(i,j,k) * (u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2))
    enddo;enddo;enddo
    
    write(filename, "(a, i5.5,a)") "data/Q",int(step),".vtk"
    open(10,file=filename,status="replace",action="write",form="unformatted",access="stream",convert="BIG_ENDIAN")
    call print_header(nx,ny,nz,x,y,z)

    write(10) 'VECTORS Velocity float'//lf
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          write(10) real(u(i,j,k)), real(v(i,j,k)), real(w(i,j,k))
    enddo;enddo;enddo

    write(10) lf//'SCALARS rho float'//lf
    write(10) 'LOOKUP_TABLE default'//lf
    write(10) real(rho), lf
    
    write(10) lf//'SCALARS P float'//lf
    write(10) 'LOOKUP_TABLE default'//lf
    write(10) real(p), lf

    write(10) lf//'SCALARS T float'//lf
    write(10) 'LOOKUP_TABLE default'//lf
    write(10) real(T), lf

    if (present(mut)) then
      nut(:,:,:) = mut(:,:,:) / rho(:,:,:)
      write(10) lf//'SCALARS nut float'//lf
      write(10) 'LOOKUP_TABLE default'//lf
      write(10) real(nut), lf
    endif
    close(10)

    call print_entropy(step,rho,p,entropy0)
    call print_KE(step,rho,u,v,w,ke0)
    call print_enstrophy(step,x,y,z,rho,u,v,w)
  end subroutine print_vtk_3D
end module print

