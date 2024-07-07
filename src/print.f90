module print
  use mod_globals, only : nt, dt, gamma, R
  implicit none
  interface
    subroutine print_entropy(step,nx,ny,nz,rho,p,entropy0,myrank)
      integer, intent(in)                       :: step, nx, ny, nz
      real(4), intent(in), dimension(nx,ny,nz)  :: rho, p
      real(4), intent(inout)                    :: entropy0
      integer, intent(in), optional             :: myrank
    end subroutine print_entropy
  end interface
  
  interface
    subroutine print_KE(step,nx,ny,nz,rho,u,v,w,ke0,myrank)
      integer, intent(in)                       :: step, nx, ny, nz
      real(4), intent(in), dimension(nx,ny,nz)  :: rho, u, v, w
      real(4), intent(inout)                    :: ke0
      integer, intent(in), optional             :: myrank
    end subroutine print_KE
  end interface

  interface
    subroutine print_enstrophy(step,nx,ny,nz,x,y,z,rho,u,v,w,myrank)
      integer, intent(in)                       :: step, nx, ny, nz
      real(4), intent(in)                       :: x(nx), y(ny), z(nz)
      real(4), intent(in), dimension(nx,ny,nz)  :: rho, u, v, w
      integer, intent(in), optional             :: myrank
    end subroutine print_enstrophy
  end interface

  interface
    subroutine print_rms(step,nx,ny,nz,u,v,w,p,myrank)
      integer, intent(in)                       :: step, nx, ny, nz
      real(4), intent(in), dimension(nx,ny,nz)  :: u, v, w, p
      integer, intent(in), optional             :: myrank
    end subroutine print_rms
  end interface

  interface
    subroutine print_xml(ni,nj,nk,dimension,x,y,z,rho1d,p1d,T1d,M1d,v1d,Q1d)
      integer, intent(in)                                 :: ni, nj, nk, dimension
      real(4), intent(in)                                 :: x(ni), y(nj), z(nk)
      real(4), intent(in), dimension(ni*nj*nk)            :: rho1d, p1d, T1d, M1d
      real(4), intent(in), dimension(dimension*ni*nj*nk)  :: v1d
      real(4), intent(in), dimension(ni*nj*nk), optional  :: Q1d
    end subroutine print_xml
  end interface

  interface print_vtk
    subroutine print_vtk_2D(step,nx,ny,x,y,Q,T)
      integer, intent(in)           :: step, nx, ny
      real(4), intent(in)           :: x(nx), y(ny)
      real(4), intent(in)           :: Q(nx,ny,4)
      real(4), intent(in), optional :: T(nx,ny)
    end subroutine print_vtk_2D

    subroutine print_vtk_3D(step,nx,ny,nz,x,y,z,Jacobian,QJ,ke0,entropy0,myrank)
      integer, intent(in)           :: step, nx, ny, nz
      real(4), intent(in)           :: x(nx), y(ny), z(nz), Jacobian(nx,ny,nz)
      real(4), intent(in)           :: QJ(nx,ny,nz,5)
      real(4), intent(inout)        :: ke0, entropy0
      integer, intent(in), optional :: myrank
    end subroutine print_vtk_3D
  end interface

  interface mean
    module procedure mean1D, mean2D, mean3D
  end interface
contains
  function mean1D(a) result(ans)
    real(4), intent(in) :: a(:)
    real(4) ans
    ans = sum(a) / size(a)
  end function mean1D

  function mean2D(a) result(ans)
    real(4), intent(in) :: a(:,:)
    real(4) ans
    ans = sum(a) / size(a)
  end function mean2D

  function mean3D(a) result(ans)
    real(4), intent(in) :: a(:,:,:)
    real(4) ans
    ans = sum(a) / size(a)
  end function mean3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine simple_bc(nx,ny,nz,A)
    integer, intent(in)    :: nx, ny, nz
    real(4), intent(inout) :: A(nx,ny,nz)
    ! x direction
    A(1,:,:)  = A(2,:,:)
    A(nx,:,:) = A(nx-1,:,:)
    ! y direction
    A(:,1,:)  = A(:,2,:)
    A(:,ny,:) = A(:,ny-1,:)
    ! z direction
    A(:,:,1)  = A(:,:,2)
    A(:,:,nz) = A(:,:,nz-1)
  end subroutine simple_bc

  function divergence(dudx,dvdy,dwdz) result(ans)
    real(4), intent(in) :: dudx, dvdy, dwdz
    real(4) ans
    ans = dudx + dvdy + dwdz
  end function divergence

  function vorticity(dudy,dudz,dvdx,dvdz,dwdx,dwdy) result(ans)
    real(4), intent(in) :: dudy, dudz, dvdx, dvdz, dwdx, dwdy
    real(4) ans(3)
    ans(:) = (/dwdy - dvdz, dudz - dwdx, dvdx - dudy/)
  end function vorticity

  function Qcriterion(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz) result(ans)
    real(4), intent(in) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    real(4) ans
    ans = (dudx * dvdy + dvdy * dwdz + dwdz * dudx) &
        - (dudy * dvdx + dvdz * dwdy + dwdx * dudz)
  end function Qcriterion
  
  subroutine calc_strain_tensor(nx,ny,nz,x,y,z,u,v,w,div,omega,Q)
    integer, intent(in)                         :: nx, ny, nz
    real(4), intent(in)                         :: x(nx), y(ny), z(nz)
    real(4), intent(in), dimension(nx,ny,nz)    :: u, v, w
    real(4), intent(out), dimension(nx,ny,nz)   :: div, Q
    real(4), intent(out), dimension(nx,ny,nz,3) :: omega
    integer i, j, k
    real(4) dudx, dvdx, dwdx, dudy, dvdy, dwdy, dudz, dvdz, dwdz
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          if (3 <= i .and. i <= nx-2) then
            dudx = (u(i-2,j,k) - 8.d0 * u(i-1,j,k) + 8.d0 * u(i+1,j,k) - u(i+2,j,k)) / (6.d0 * (-x(i-1) + x(i+1)))
            dvdx = (v(i-2,j,k) - 8.d0 * v(i-1,j,k) + 8.d0 * v(i+1,j,k) - v(i+2,j,k)) / (6.d0 * (-x(i-1) + x(i+1)))
            dwdx = (w(i-2,j,k) - 8.d0 * w(i-1,j,k) + 8.d0 * w(i+1,j,k) - w(i+2,j,k)) / (6.d0 * (-x(i-1) + x(i+1)))
          else
            dudx = (-u(i-1,j,k) + u(i+1,j,k)) / (-x(i-1) + x(i+1))
            dvdx = (-v(i-1,j,k) + v(i+1,j,k)) / (-x(i-1) + x(i+1))
            dwdx = (-w(i-1,j,k) + w(i+1,j,k)) / (-x(i-1) + x(i+1))
          endif
          if (3 <= j .and. j <= ny-2) then
            dudy = (u(i,j-2,k) - 8.d0 * u(i,j-1,k) + 8.d0 * u(i,j+1,k) - u(i,j+2,k)) / (6.d0 * (-y(j-1) + y(j+1)))
            dvdy = (v(i,j-2,k) - 8.d0 * v(i,j-1,k) + 8.d0 * v(i,j+1,k) - v(i,j+2,k)) / (6.d0 * (-y(j-1) + y(j+1)))
            dwdy = (w(i,j-2,k) - 8.d0 * w(i,j-1,k) + 8.d0 * w(i,j+1,k) - w(i,j+2,k)) / (6.d0 * (-y(j-1) + y(j+1)))
          else
            dudy = (-u(i,j-1,k) + u(i,j+1,k)) / (-y(j-1) + y(j+1))
            dvdy = (-v(i,j-1,k) + v(i,j+1,k)) / (-y(j-1) + y(j+1))
            dwdy = (-w(i,j-1,k) + w(i,j+1,k)) / (-y(j-1) + y(j+1))
          endif
          if (3 <= k .and. k <= nz-2) then
            dudz = (u(i,j,k-2) - 8.d0 * u(i,j,k-1) + 8.d0 * u(i,j,k+1) - u(i,j,k+2)) / (6.d0 * (-z(k-1) + z(k+1)))
            dvdz = (v(i,j,k-2) - 8.d0 * v(i,j,k-1) + 8.d0 * v(i,j,k+1) - v(i,j,k+2)) / (6.d0 * (-z(k-1) + z(k+1)))
            dwdz = (w(i,j,k-2) - 8.d0 * w(i,j,k-1) + 8.d0 * w(i,j,k+1) - w(i,j,k+2)) / (6.d0 * (-z(k-1) + z(k+1)))
          else
            dudz = (-u(i,j,k-1) + u(i,j,k+1)) / (-z(k-1) + z(k+1))
            dvdz = (-v(i,j,k-1) + v(i,j,k+1)) / (-z(k-1) + z(k+1))
            dwdz = (-w(i,j,k-1) + w(i,j,k+1)) / (-z(k-1) + z(k+1))
          endif
          div(i,j,k) = divergence(dudx,dvdy,dwdz)
          omega(i,j,k,:) = vorticity(dudy,dudz,dvdx,dvdz,dwdx,dwdy)
          Q(i,j,k) = Qcriterion(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
    enddo;enddo;enddo
    ! set boundary condition
    call simple_bc(nx,ny,nz,div)
    call simple_bc(nx,ny,nz,omega(:,:,:,1))
    call simple_bc(nx,ny,nz,omega(:,:,:,2))
    call simple_bc(nx,ny,nz,omega(:,:,:,3))
    call simple_bc(nx,ny,nz,Q)
  end subroutine calc_strain_tensor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_entropy(step,nx,ny,nz,rho,p,entropy0,myrank)
    integer, intent(in)                       :: step, nx, ny, nz
    real(4), intent(in), dimension(nx,ny,nz)  :: rho, p
    real(4), intent(inout)                    :: entropy0
    integer, intent(in), optional             :: myrank
    real(4) entropy, t
    character(len=40) filename
    entropy = sum(rho(3:nx-2,3:ny-2,3:nz-2) * &
              log(p(3:nx-2,3:ny-2,3:nz-2) * rho(3:nx-2,3:ny-2,3:nz-2) ** (-gamma)))
    if (step == 0) then
      entropy0 = entropy
    endif
    t = nt * step * dt
    if (present(myrank)) then
      write(filename, "(a, i1.1, a)") "data/",int(myrank),"/entropy.d"
      open(10,file=filename, position="append")
    else
      open(10,file="data/entropy.d", position="append")
    endif
    write(10,"(2(f9.4,1x))") t, (entropy0 - entropy) / entropy0
    close(10)
  end subroutine print_entropy

  subroutine print_KE(step,nx,ny,nz,rho,u,v,w,ke0,myrank)
    integer, intent(in)                       :: step, nx, ny, nz
    real(4), intent(in), dimension(nx,ny,nz)  :: rho, u, v, w
    real(4), intent(inout)                    :: ke0
    integer, intent(in), optional             :: myrank
    real(4) ke, t
    character(len=40) filename
    ke = mean(0.5e0 * rho(3:nx-2,3:ny-2,3:nz-2) * &
         (u(3:nx-2,3:ny-2,3:nz-2)**2 + v(3:nx-2,3:ny-2,3:nz-2)**2 + w(3:nx-2,3:ny-2,3:nz-2)**2))
    if (step == 0) then
      ke0 = ke
    endif
    t = nt * step * dt
    if (present(myrank)) then
      write(filename, "(a, i1.1, a)") "data/",int(myrank),"/kinetic_energy.d"
      open(10,file=filename, position="append")
    else
      open(10,file="data/kinetic_energy.d", position="append")
    endif
    !write(10,"(2(f9.4,1x))") t, ke
    write(10,"(2e12.4)") t, ke / ke0
    close(10)
  end subroutine print_KE

  subroutine print_enstrophy(step,nx,ny,nz,x,y,z,rho,omega,myrank)
    integer, intent(in)                        :: step, nx, ny, nz
    real(4), intent(in)                        :: x(nx), y(ny), z(nz)
    real(4), intent(in), dimension(nx,ny,nz)   :: rho
    real(4), intent(in), dimension(nx,ny,nz,3) :: omega
    integer, intent(in), optional              :: myrank
    real(4) enstrophy, t
    character(len=40) filename
    t = nt * step * dt
    enstrophy = mean(0.5e0 * rho(:,:,:) * (omega(:,:,:,1)**2 + omega(:,:,:,2)**2 + omega(:,:,:,3)**2))
    if (present(myrank)) then
      write(filename, "(a, i1.1, a)") "data/",int(myrank),"/enstrophy.d"
      open(10,file=filename, position="append")
    else
      open(10,file="data/enstrophy.d", position="append")
    endif
    !write(10,"(2(f9.4,1x))") t, enstrophy
    write(10,"(2e12.4)") t, enstrophy
    close(10)
  end subroutine print_enstrophy

  subroutine print_rms(step,nx,ny,nz,u,v,w,p,myrank)
    integer, intent(in)                       :: step, nx, ny, nz
    real(4), intent(in), dimension(nx,ny,nz)  :: u, v, w, p
    integer, intent(in), optional             :: myrank
    real(4) t, umean, vmean, wmean, pmean, urms, vrms, wrms, prms, NxNyNz
    integer i, j, k
    character(len=40) filename
    t = nt * step * dt
    NxNyNz = real(nx * ny * nz)
    umean = sum(u) / NxNyNz
    vmean = sum(v) / NxNyNz
    wmean = sum(w) / NxNyNz
    pmean = sum(p) / NxNyNz

    urms = sqrt(sum(u**2) / NxNyNz - umean**2)
    vrms = sqrt(sum(v**2) / NxNyNz - vmean**2)
    wrms = sqrt(sum(w**2) / NxNyNz - wmean**2)
    prms = sqrt(sum(p**2) / NxNyNz - pmean**2)

    if (present(myrank)) then
      write(filename, "(a, i1.1, a)") "data/",int(myrank),"/rms.d"
      open(10,file=filename, position="append")
    else
      open(10,file="data/rms.d", position="append")
    endif
    write(10,"(5e12.4)") t, urms, vrms, wrms, prms
    close(10)
  end subroutine print_rms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function mu(T) result(ans)
    real(4), intent(in), value :: T
    real(4) :: ans
    real(4) :: mu0 = 1.716e-5
    real(4) :: T0 = 273.2e0
    real(4) :: S = 111.e0
    ans = mu0 * ((T0 + S) / (T + S)) * (T / T0)**1.5
  end function mu
  
  subroutine print_turbulent_boundary_layer(step,nx,ny,nz,dy,y,T,u,rho)
    integer, intent(in), value                :: step, nx, ny, nz
    real(4), intent(in), value                :: dy
    real(4), intent(in), dimension(ny)        :: y
    real(4), intent(in), dimension(nx,nz)     :: T
    real(4), intent(in), dimension(nx,ny,nz)  :: u, rho
    integer j
    real(4) rhow, nuw, dudy, tw, ut, uvd
    real(4), dimension(ny) :: yplus, uplus
    rhow = mean(rho(:,1,:))
    nuw = mu(mean(T(:,:))) / rhow
    dudy = dy * mean(-u(:,1,:) + u(:,2,:))
    tw = rhow * nuw * dudy
    ut = sqrt(tw / rhow)
    open(10,file="data/yplus.d",action="write")
    yplus(1) = ut * y(1) / nuw
    uplus(1) = mean(u(:,1,:)) / ut
    do j = 2, ny
      yplus(j) = ut * y(j) / nuw
      ! van Driest transformation
      uvd = mean(u(:,j-1,:)) + sqrt(mean(rho(:,j,:)) / rhow) * (-mean(u(:,j-1,:)) + mean(u(:,j,:)))
      uplus(j) = uvd / ut
      write(10,"(2e12.4)") yplus(j), uplus(j)
    enddo
    close(10)
    ! tau
    open(10,file="data/tau.d", position="append")
    write(10,"(2(f9.4,1x))") nt * dt * step, tw
    close(10)
  end subroutine print_turbulent_boundary_layer

  subroutine print_boundary_layer(nx,ny,nz,y,u)
    integer, intent(in), value                :: nx, ny, nz
    real(4), intent(in), dimension(ny)        :: y
    real(4), intent(in), dimension(nx,ny,nz)  :: u
    integer j
    open(10,file="data/boundary_layer.d",action="write")
    do j = 1, ny
      write(10,"(2(f12.7,1x))") y(j), mean(u(:,j,:))/mean(u(:,ny,:))
    enddo
    close(10)
  end subroutine print_boundary_layer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_xml(ni,nj,nk,dimension,x,y,z,rho1d,p1d,T1d,M1d,v1d,Q1d)
    integer, intent(in)                                 :: ni, nj, nk, dimension
    real(4), intent(in)                                 :: x(ni), y(nj), z(nk)
    real(4), intent(in), dimension(ni*nj*nk)            :: rho1d, p1d, T1d, M1d
    real(4), intent(in), dimension(dimension*ni*nj*nk)  :: v1d
    real(4), intent(in), dimension(ni*nj*nk), optional  :: Q1d
    integer(4) byte_x, byte_y, byte_z, byte_rho, byte_p, byte_T, byte_M, byte_v, byte_Q
    character :: lf*1, str1*4, str2*4, str3*4, str4*1
    character :: offset1*10, offset2*10, offset3*10, offset4*10, offset5*10, offset6*10, offset7*10, offset8*10
    lf = char(10)
    write(str1(1:4),'(i4)') ni-1
    write(str2(1:4),'(i4)') nj-1
    write(str3(1:4),'(i4)') nk-1
    write(str4(1:1),'(i1)') dimension
    byte_x   = 4 + 4 * ni
    byte_y   = 4 + 4 * nj
    byte_z   = 4 + 4 * nk
    byte_rho = 4 + 4 * (ni * nj * nk)
    byte_p   = byte_rho
    byte_T   = byte_rho
    byte_M   = byte_rho
    byte_v   = 4 + 4 * (dimension * ni * nj * nk)
    byte_Q   = byte_rho
    write(offset1(1:10),'(i10)') byte_x
    write(offset2(1:10),'(i10)') byte_x + byte_y
    write(offset3(1:10),'(i10)') byte_x + byte_y + byte_z
    write(offset4(1:10),'(i10)') byte_x + byte_y + byte_z + byte_rho
    write(offset5(1:10),'(i10)') byte_x + byte_y + byte_z + byte_rho + byte_p
    write(offset6(1:10),'(i10)') byte_x + byte_y + byte_z + byte_rho + byte_p + byte_T
    write(offset7(1:10),'(i10)') byte_x + byte_y + byte_z + byte_rho + byte_P + byte_T + byte_M
    write(offset8(1:10),'(i10)') byte_x + byte_y + byte_z + byte_rho + byte_P + byte_T + byte_M + byte_v

    write(10) '<?xml version="1.0"?>'//lf
    write(10) '<VTKFile type="RectilinearGrid" version="1.0" byte_order="LittleEndian">'//lf
    write(10) '  <RectilinearGrid WholeExtent="0 '//str1//' 0 '//str2//' 0 '//str3//'">'//lf
    write(10) '    <Piece Extent="0 '//str1//' 0 '//str2//' 0 '//str3//'">'//lf
    write(10) '      <Coordinates>'//lf
    write(10) '        <DataArray type="Float32" format="appended" offset="0"/>'//lf
    write(10) '        <DataArray type="Float32" format="appended" offset="'//offset1//'"/>'//lf
    write(10) '        <DataArray type="Float32" format="appended" offset="'//offset2//'"/>'//lf
    write(10) '      </Coordinates>'//lf
    write(10) '      <PointData>'//lf
    write(10) '        <DataArray type="Float32" format="appended" offset="'//offset3//'"&
    & Name="rho" NumberOfComponents="1"/>'//lf
    write(10) '        <DataArray type="Float32" format="appended" offset="'//offset4//'"&
    & Name="p" NumberOfComponents="1"/>'//lf
    write(10) '        <DataArray type="Float32" format="appended" offset="'//offset5//'"&
    & Name="T" NumberOfComponents="1"/>'//lf
    write(10) '        <DataArray type="Float32" format="appended" offset="'//offset6//'"&
    & Name="M" NumberOfComponents="1"/>'//lf
    write(10) '        <DataArray type="Float32" format="appended" offset="'//offset7//'"&
    & Name="velocity" NumberOfComponents="'//str4//'"/>'//lf
    if (present(Q1d)) then
    write(10) '        <DataArray type="Float32" format="appended" offset="'//offset8//'"&
    & Name="Qcriterion" NumberOfComponents="1"/>'//lf
    endif
    write(10) '      </PointData>'//lf
    write(10) '    </Piece>'//lf
    write(10) '  </RectilinearGrid>'//lf
    write(10) '  <AppendedData encoding="raw">'//lf
    if (present(Q1d)) then
      write(10) '  _', byte_x, x, byte_y, y, byte_z, z,&
      & byte_rho, rho1d, byte_p, p1d, byte_T, T1d, byte_M, M1d, byte_v, v1d, byte_Q, Q1d, lf
    else
      write(10) '  _', byte_x, x, byte_y, y, byte_z, z,&
      & byte_rho, rho1d, byte_p, p1d, byte_T, T1d, byte_M, M1d, byte_v, v1d, lf
    endif
    write(10) '  </AppendedData>'//lf
    write(10) '</VTKFile>'//lf
    close(10)
  end subroutine print_xml

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_vtk_2D(step,nx,ny,x,y,Q,T)
    integer, intent(in)                   :: step, nx, ny
    real(4), intent(in)                   :: x(nx), y(ny)
    real(4), intent(in)                   :: Q(nx,ny,4)
    real(4), intent(in), optional         :: T(nx,ny)
    integer i, j, l, m
    real(4), dimension(nx,ny) :: rho, u, v, p
    real(4) :: z(1) = 0.e0
    real(4), dimension(nx*ny)   :: rho1d, p1d, T1d, M1d
    real(4), dimension(2*nx*ny) :: v1d
    character(len=40) filename
    character :: lf*1
    lf = char(10)
    l = 1
    m = 1
    do j = 1, ny
      do i = 1, nx
        rho(i,j) = Q(i,j,1)
        u(i,j) = Q(i,j,2) / rho(i,j)
        v(i,j) = Q(i,j,3) / rho(i,j)
        p(i,j) = (gamma - 1.e0) * (Q(i,j,4) - 0.5e0 * rho(i,j) * (u(i,j)**2 + v(i,j)**2))
        rho1d(l) = rho(i,j)
        p1d(l)   = p(i,j)
        T1d(l)   = p1d(l) / (real(R) * rho1d(l))
        v1d(m)   = u(i,j)
        v1d(m+1) = v(i,j)
        M1d(l)   = sqrt(v1d(m)**2 + v1d(m+1)**2) / sqrt(gamma * p1d(l) / rho1d(l))
        l = l + 1
        m = m + 2
    enddo;enddo

    write(filename, "(a, i5.5,a)") "data/Q",int(step),".vtr"
    open(10,file=filename,status="replace",action="write",form="unformatted",access="stream",convert="Little_ENDIAN")
    call print_xml(nx,ny,1,2,x,y,z,rho1d,p1d,T1d,M1d,v1d)
  end subroutine print_vtk_2D
  
  subroutine print_vtk_3D(step,nx,ny,nz,x,y,z,Jacobian,QJ,ke0,entropy0,myrank)
    integer, intent(in)           :: step, nx, ny, nz
    real(4), intent(in)           :: x(nx), y(ny), z(nz), Jacobian(nx,ny,nz)
    real(4), intent(in)           :: QJ(nx,ny,nz,5) ! Q / Jacobian
    real(4), intent(inout)        :: ke0, entropy0
    integer, intent(in), optional :: myrank
    integer i, j, k, l, m, len
    real(4) dy
    real(4), allocatable :: rho(:,:,:), u(:,:,:), v(:,:,:), w(:,:,:), p(:,:,:), div(:,:,:), omega(:,:,:,:), Qcriterion(:,:,:), Tw(:,:)
    real(4), allocatable :: rho1d(:), p1d(:), T1d(:), M1d(:), v1d(:), div1d(:), omega1d(:), Qcriterion1d(:)
    character(len=40) filename
    character :: lf*1
    lf = char(10)
    len = nx * ny * nz
    allocate(rho(nx,ny,nz),u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),p(nx,ny,nz),Tw(nx,nz),rho1d(len),p1d(len),T1d(len),M1d(len),v1d(3*len))
    allocate(div(nx,ny,nz),omega(nx,ny,nz,3),Qcriterion(nx,ny,nz),div1d(len),omega1d(3*len),Qcriterion1d(len))
    l = 1
    m = 1
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rho(i,j,k) = Jacobian(i,j,k) * QJ(i,j,k,1)
          u(i,j,k) = Jacobian(i,j,k) * QJ(i,j,k,2) / rho(i,j,k)
          v(i,j,k) = Jacobian(i,j,k) * QJ(i,j,k,3) / rho(i,j,k)
          w(i,j,k) = Jacobian(i,j,k) * QJ(i,j,k,4) / rho(i,j,k)
          p(i,j,k) = (gamma - 1.e0) * (Jacobian(i,j,k) * QJ(i,j,k,5) - 0.e0 * rho(i,j,k) * (u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2))
          rho1d(l) = rho(i,j,k)
          p1d(l)   = p(i,j,k)
          T1d(l)   = p1d(l) / (real(R) * rho1d(l))
          v1d(m)   = u(i,j,k)
          v1d(m+1) = v(i,j,k)
          v1d(m+2) = w(i,j,k)
          M1d(l)   = sqrt(v1d(m)**2 + v1d(m+1)**2 + v1d(m+2)**2) / sqrt(gamma * p1d(l) / rho1d(l))
          l = l + 1
          m = m + 3
    enddo;enddo;enddo

    ! strain tensor
    call calc_strain_tensor(nx,ny,nz,x,y,z,u,v,w,div,omega,Qcriterion)
    l = 1
    m = 1
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          div1d(l)        = div(i,j,k)
          omega1d(m)      = omega(i,j,k,1)
          omega1d(m+1)    = omega(i,j,k,2)
          omega1d(m+2)    = omega(i,j,k,3)
          Qcriterion1d(l) = Qcriterion(i,j,k)
          l = l + 1
          m = m + 3
    enddo;enddo;enddo
    
    if (present(myrank)) then
      call print_entropy(step,nx,ny,nz,rho,p,entropy0,myrank)
      call print_KE(step,nx,ny,nz,rho,u,v,w,ke0,myrank)
      call print_enstrophy(step,nx,ny,nz,real(x),real(y),real(z),rho,omega,myrank)
      call print_rms(step,nx,ny,nz,u,v,w,p,myrank)
      if (myrank == 3) then
        call print_boundary_layer(nx,ny,nz,real(y),u)
        dy = 1.e0 / (-y(1) + y(2))
        Tw(:,:) = p(:,1,:) / (R * rho(:,1,:))
        call print_turbulent_boundary_layer(step,nx,ny,nz,dy,real(y),Tw,u,rho)
      endif
      write(filename, "(a, i1.1, a, i5.5, a)") "data/",int(myrank),"/Q",int(step),".vtr"
    else
      call print_entropy(step,nx,ny,nz,rho,p,entropy0)
      call print_KE(step,nx,ny,nz,rho,u,v,w,ke0)
      call print_enstrophy(step,nx,ny,nz,real(x),real(y),real(z),rho,omega)
      call print_rms(step,nx,ny,nz,u,v,w,p)
      call print_boundary_layer(nx,ny,nz,real(y),u)
      dy = 1.e0 / (-y(1) + y(2))
      Tw(:,:) = p(:,1,:) / (R * rho(:,1,:))
      call print_turbulent_boundary_layer(step,nx,ny,nz,dy,real(y),Tw,u,rho)
      write(filename, "(a, i5.5, a)") "data/Q",int(step),".vtr"
    endif
    open(10,file=filename,status="replace",action="write",form="unformatted",access="stream",convert="Little_ENDIAN")
    call print_xml(nx,ny,nz,3,real(x),real(y),real(z),rho1d,p1d,T1d,M1d,v1d,Qcriterion1d)
    
    deallocate(rho,u,v,w,p,div,omega,Qcriterion,Tw,rho1d,p1d,T1d,M1d,v1d,div1d,omega1d,Qcriterion1d)
  end subroutine print_vtk_3D
end module print

