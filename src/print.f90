module print
  use mod_globals, only : id_accuracy, nt, np, dt, step_offset, gamma, R, Lx
  implicit none
  
  interface
    subroutine print_mass(step,nx,ny,nz,rho,u,v,w,mass0,myrank)
      integer, intent(in)                       :: step, nx, ny, nz
      real(8), intent(in), dimension(nx,ny,nz)  :: rho, u, v, w
      real(8), intent(inout)                    :: mass0
      integer, intent(in), optional             :: myrank
    end subroutine print_mass
  end interface

  interface
    subroutine print_entropy(step,nx,ny,nz,rho,p,entropy0,myrank)
      integer, intent(in)                       :: step, nx, ny, nz
      real(8), intent(in), dimension(nx,ny,nz)  :: rho, p
      real(8), intent(inout)                    :: entropy0
      integer, intent(in), optional             :: myrank
    end subroutine print_entropy
  end interface
  
  interface
    subroutine print_KE(step,nx,ny,nz,rho,u,v,w,ke0,myrank)
      integer, intent(in)                       :: step, nx, ny, nz
      real(8), intent(in), dimension(nx,ny,nz)  :: rho, u, v, w
      real(8), intent(inout)                    :: ke0
      integer, intent(in), optional             :: myrank
    end subroutine print_KE
  end interface

  interface
    subroutine print_enstrophy(step,nx,ny,nz,x,y,z,rho,u,v,w,myrank)
      integer, intent(in)                       :: step, nx, ny, nz
      real(8), intent(in)                       :: x(nx), y(ny), z(nz)
      real(8), intent(in), dimension(nx,ny,nz)  :: rho, u, v, w
      integer, intent(in), optional             :: myrank
    end subroutine print_enstrophy
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
      real(8), intent(in)           :: x(nx), y(ny)
      real(8), intent(in)           :: Q(nx,ny,4)
      real(8), intent(in), optional :: T(nx,ny)
    end subroutine print_vtk_2D

    subroutine print_vtk_3D(step,nx,ny,nz,x,y,z,Jacobian,QJ,ke0,entropy0,myrank)
      integer, intent(in)           :: step, nx, ny, nz
      real(8), intent(in)           :: x(nx), y(ny), z(nz), Jacobian(nx,ny,nz)
      real(8), intent(in)           :: QJ(nx,ny,nz,5)
      real(8), intent(inout)        :: ke0, entropy0
      integer, intent(in), optional :: myrank
    end subroutine print_vtk_3D
  end interface

  interface mean
    module procedure mean1D, mean2D, mean3D
  end interface
contains
  function mean1D(a) result(ans)
    real(8), intent(in) :: a(:)
    real(8) ans
    ans = sum(a) / dble(size(a))
  end function mean1D

  function mean2D(a) result(ans)
    real(8), intent(in) :: a(:,:)
    real(8) ans
    ans = sum(a) / dble(size(a))
  end function mean2D

  function mean3D(a) result(ans)
    real(8), intent(in) :: a(:,:,:)
    real(8) ans
    ans = sum(a) / dble(size(a))
  end function mean3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine simple_bc(nx,ny,nz,A)
    integer, intent(in)    :: nx, ny, nz
    real(8), intent(inout) :: A(nx,ny,nz)
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
    real(8), intent(in) :: dudx, dvdy, dwdz
    real(8) ans
    ans = dudx + dvdy + dwdz
  end function divergence

  function vorticity(dudy,dudz,dvdx,dvdz,dwdx,dwdy) result(ans)
    real(8), intent(in) :: dudy, dudz, dvdx, dvdz, dwdx, dwdy
    real(8) ans(3)
    ans(:) = (/dwdy - dvdz, dudz - dwdx, dvdx - dudy/)
  end function vorticity

  function Qcriterion(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz) result(ans)
    real(8), intent(in) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    real(8) ans
    ans = (dudx * dvdy + dvdy * dwdz + dwdz * dudx) &
        - (dudy * dvdx + dvdz * dwdy + dwdx * dudz)
  end function Qcriterion
  
  subroutine calc_strain_tensor(nx,ny,nz,x,y,z,u,v,w,div,omega,Q)
    integer, intent(in)                         :: nx, ny, nz
    real(8), intent(in)                         :: x(nx), y(ny), z(nz)
    real(8), intent(in), dimension(nx,ny,nz)    :: u, v, w
    real(8), intent(out), dimension(nx,ny,nz)   :: div, Q
    real(8), intent(out), dimension(nx,ny,nz,3) :: omega
    integer i, j, k
    real(8) dudx, dvdx, dwdx, dudy, dvdy, dwdy, dudz, dvdz, dwdz
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

  subroutine print_mass(step,nx,ny,nz,rho,u,v,w,mass0,myrank)
    integer, intent(in)                       :: step, nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz)  :: rho, u, v, w
    real(8), intent(inout)                    :: mass0
    integer, intent(in), optional             :: myrank
    real(8) mass, t
    character(len=40) filename
    mass = sum(rho(3:nx-2,3:ny-2,3:nz-2))
    if (step == 0) then
      mass0 = mass
    endif
    t = nt * step * dt
    if (present(myrank)) then
      write(filename, "(a, i1.1, a)") "data/",int(myrank),"/mass.d"
      open(10,file=filename, position="append")
    else
      open(10,file="data/mass.d", position="append")
    endif
    write(10,"(2e12.4)") t, (mass0 - mass) / mass0
    close(10)
  end subroutine print_mass

  subroutine print_entropy(step,nx,ny,nz,rho,p,entropy0,myrank)
    integer, intent(in)                       :: step, nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz)  :: rho, p
    real(8), intent(inout)                    :: entropy0
    integer, intent(in), optional             :: myrank
    real(8) entropy, t
    character(len=40) filename
    integer i, j, k, accuracy, offset
    if (kind(id_accuracy) == 8) then
      accuracy = 6
      offset   = accuracy / 2
    elseif (kind(id_accuracy) == 4) then
      accuracy = 4
      offset   = accuracy / 2
    else
      accuracy = 2
      offset   = accuracy / 2
    endif
    entropy = 0.d0
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          entropy = entropy + rho(i,j,k) * log(p(i,j,k) * (rho(i,j,k)**(-gamma)))
    enddo;enddo;enddo
    entropy = entropy / dble((nx-accuracy) * (ny-accuracy) * (nz-accuracy))

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
    write(10,"(2e12.4)") t, (entropy0 - entropy) / entropy0
    close(10)
  end subroutine print_entropy

  subroutine print_KE(step,nx,ny,nz,rho,u,v,w,ke0,myrank)
    integer, intent(in)                       :: step, nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz)  :: rho, u, v, w
    real(8), intent(inout)                    :: ke0
    integer, intent(in), optional             :: myrank
    real(8) ke, t
    character(len=40) filename
    integer i, j, k, accuracy, offset
    if (kind(id_accuracy) == 8) then
      accuracy = 6
      offset   = accuracy / 2
    elseif (kind(id_accuracy) == 4) then
      accuracy = 4
      offset   = accuracy / 2
    else
      accuracy = 2
      offset   = accuracy / 2
    endif
    ke = 0.d0
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          ke = ke + 0.5d0 * rho(i,j,k) * (u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2)
    enddo;enddo;enddo
    ke = ke / dble((nx-accuracy) * (ny-accuracy) * (nz-accuracy))

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
    write(10,"(3e12.4)") t, ke, ke / ke0
    close(10)
  end subroutine print_KE

  subroutine print_enstrophy(step,nx,ny,nz,x,y,z,rho,omega,myrank)
    integer, intent(in)                        :: step, nx, ny, nz
    real(8), intent(in)                        :: x(nx), y(ny), z(nz)
    real(8), intent(in), dimension(nx,ny,nz)   :: rho
    real(8), intent(in), dimension(nx,ny,nz,3) :: omega
    integer, intent(in), optional              :: myrank
    real(8) enstrophy, t
    character(len=40) filename
    integer i, j, k, accuracy, offset
    if (kind(id_accuracy) == 8) then
      accuracy = 6
      offset   = accuracy / 2
    elseif (kind(id_accuracy) == 4) then
      accuracy = 4
      offset   = accuracy / 2
    else
      accuracy = 2
      offset   = accuracy / 2
    endif
    t = nt * step * dt
    enstrophy = 0.d0
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          enstrophy = enstrophy + 0.5d0 * rho(i,j,k) * (omega(i,j,k,1)**2 + omega(i,j,k,2)**2 + omega(i,j,k,3)**2)
    enddo;enddo;enddo
    enstrophy = enstrophy / dble((nx-accuracy) * (ny-accuracy) * (nz-accuracy))

    if (present(myrank)) then
      write(filename, "(a, i1.1, a)") "data/",int(myrank),"/enstrophy.d"
      open(10,file=filename, position="append")
    else
      open(10,file="data/enstrophy.d", position="append")
    endif
    write(10,"(2e12.4)") t, enstrophy
    close(10)
  end subroutine print_enstrophy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_1d(step,nx,ny,nz,x,y,z,rho,p,u,sensor)
    integer, intent(in), value               :: step, nx, ny, nz
    real(8), intent(in)                      :: x(nx), y(ny), z(nz)
    real(8), intent(in), dimension(nx,ny,nz) :: rho, p, u, sensor
    real(8) T
    integer i, nyh, nzh
    character(len=40) filename
    nyh = int(0.5 * ny)
    nzh = int(0.5 * nz)
    write(filename, "(a, i5.5, a)") "data/1d/Q", int(step), ".d"
    open(10,file=filename)
    do i = 1, nx
      T = p(i,nyh,nzh) / (rho(i,nyh,nzh) * R)
      write(10,"(7e12.4)") x(i) / Lx, rho(i,nyh,nzh), u(i,nyh,nzh), p(i,nyh,nzh), T, sensor(i,nyh,nzh)
    enddo
    close(10)
  end subroutine print_1d

  function mu(T) result(ans)
    real(8), intent(in), value :: T
    real(8) :: ans, mu0 = 1.716d-5, T0 = 273.2d0, S = 111.d0
    ans = mu0 * ((T0 + S) / (T + S)) * (T / T0)**1.5
  end function mu
  
  subroutine print_turbulent_boundary_layer(step,nx,ny,nz,dy,y,T,u,rho)
    integer, intent(in), value                :: step, nx, ny, nz
    real(8), intent(in), value                :: dy
    real(8), intent(in), dimension(ny)        :: y
    real(8), intent(in), dimension(nx,nz)     :: T
    real(8), intent(in), dimension(nx,ny,nz)  :: u, rho
    integer j
    real(8) rhow, nuw, dudy, tw, ut, uvd
    real(8), dimension(ny) :: yplus, uplus
    rhow = mean(rho(:,1,:))
    nuw  = mu(mean(T(:,:))) / rhow
    dudy = dy * mean(-u(:,2,:) + u(:,3,:))
    tw   = rhow * nuw * dudy
    ut   = sqrt(tw / rhow)
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
    write(10,"(3(f9.4,1x))") nt * dt * step, tw, ut
    close(10)
  end subroutine print_turbulent_boundary_layer

  subroutine print_boundary_layer(nx,ny,nz,y,u)
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), dimension(ny)        :: y
    real(8), intent(in), dimension(nx,ny,nz)  :: u
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
    real(8), intent(in)                   :: x(nx), y(ny)
    real(8), intent(in)                   :: Q(nx,ny,4)
    real(8), intent(in), optional         :: T(nx,ny)
    integer i, j, l, m
    real(8), dimension(nx,ny) :: rho, u, v, p
    real(8) :: z(1) = 0.d0
    real(8), dimension(nx*ny)   :: rho1d, p1d, T1d, M1d
    real(8), dimension(2*nx*ny) :: v1d
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
        p(i,j) = (gamma - 1.d0) * (Q(i,j,4) - 0.5d0 * rho(i,j) * (u(i,j)**2 + v(i,j)**2))
        rho1d(l) = rho(i,j)
        p1d(l)   = p(i,j)
        T1d(l)   = p1d(l) / (R * rho1d(l))
        v1d(m)   = u(i,j)
        v1d(m+1) = v(i,j)
        M1d(l)   = sqrt(v1d(m)**2 + v1d(m+1)**2) / sqrt(gamma * p1d(l) / rho1d(l))
        l = l + 1
        m = m + 2
    enddo;enddo

    write(filename, "(a, i5.5,a)") "data/Q",int(step),".vtr"
    open(10,file=filename,status="replace",action="write",form="unformatted",access="stream",convert="Little_ENDIAN")
    call print_xml(nx,ny,1,2,real(x),real(y),real(z),real(rho1d),real(p1d),real(T1d),real(M1d),real(v1d))
  end subroutine print_vtk_2D
  
  subroutine print_vtk_3D(step,nx,ny,nz,x,y,z,Jacobian,QJ,rhom,pm,Tm,Mm,vm,mass0,ke0,entropy0,myrank)
    integer, intent(in)           :: step, nx, ny, nz
    real(8), intent(in)           :: x(nx), y(ny), z(nz), Jacobian(nx,ny,nz)
    real(8), intent(in)           :: QJ(nx,ny,nz,5) ! Q / Jacobian
    real(8), intent(inout)        :: rhom(nx*ny*nz), pm(nx*ny*nz), Tm(nx*ny*nz), Mm(nx*ny*nz), vm(3*nx*ny*nz)
    real(8), intent(inout)        :: mass0, ke0, entropy0
    integer, intent(in), optional :: myrank
    integer i, j, k, l, m, len
    real(8) dy
    real(8), allocatable :: rho(:,:,:), u(:,:,:), v(:,:,:), w(:,:,:), p(:,:,:), div(:,:,:), omega(:,:,:,:), Qcriterion(:,:,:), Tw(:,:)
    real(8), allocatable :: rho1d(:), p1d(:), T1d(:), M1d(:), v1d(:), div1d(:), omega1d(:), Qcriterion1d(:)
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
          u(i,j,k)   = QJ(i,j,k,2) / QJ(i,j,k,1)
          v(i,j,k)   = QJ(i,j,k,3) / QJ(i,j,k,1)
          w(i,j,k)   = QJ(i,j,k,4) / QJ(i,j,k,1)
          p(i,j,k)   = (gamma - 1.d0) * (Jacobian(i,j,k) * QJ(i,j,k,5) - 0.5d0 * rho(i,j,k) * (u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2))
          rho1d(l)   = rho(i,j,k)
          p1d(l)     = p(i,j,k)
          T1d(l)     = p1d(l) / (R * rho1d(l))
          v1d(m)     = u(i,j,k)
          v1d(m+1)   = v(i,j,k)
          v1d(m+2)   = w(i,j,k)
          M1d(l)     = sqrt(v1d(m)**2 + v1d(m+1)**2 + v1d(m+2)**2) / sqrt(gamma * p1d(l) / rho1d(l))
          rhom(l)    = rhom(l) + rho1d(l)
          pm(l)      = pm(l)   + p1d(l)
          Tm(l)      = Tm(l)   + T1d(l)
          Mm(l)      = Mm(l)   + M1d(l)
          vm(m)      = vm(m)   + v1d(m)
          vm(m+1)    = vm(m+1) + v1d(m+1)
          vm(m+2)    = vm(m+2) + v1d(m+2)
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
      !call print_mass(step,nx,ny,nz,rho,u,v,w,mass0,myrank)
      call print_entropy(step,nx,ny,nz,rho,p,entropy0,myrank)
      call print_KE(step,nx,ny,nz,rho,u,v,w,ke0,myrank)
      call print_enstrophy(step,nx,ny,nz,x,y,z,rho,omega,myrank)
      !call print_1d(step,nx,ny,nz,x,y,z,rho,p,u,sensor)
      if (myrank == 3) then
        !call print_boundary_layer(nx,ny,nz,y,u)
        dy = 1.d0 / (-y(1) + y(2))
        Tw(:,:) = p(:,1,:) / (R * rho(:,1,:))
        call print_turbulent_boundary_layer(step,nx,ny,nz,dy,y,Tw,u,rho)
      endif
      write(filename, "(a, i1.1, a, i5.5, a)") "data/",int(myrank),"/Q",int(step+step_offset),".vtr"
    else
      !call print_mass(step,nx,ny,nz,rho,u,v,w,mass0)
      call print_entropy(step,nx,ny,nz,rho,p,entropy0)
      call print_KE(step,nx,ny,nz,rho,u,v,w,ke0)
      call print_enstrophy(step,nx,ny,nz,x,y,z,rho,omega)
      !call print_1d(step,nx,ny,nz,x,y,z,rho,p,u,sensor)
      !call print_boundary_layer(nx,ny,nz,y,u)
      dy = 1.d0 / (-y(2) + y(3))
      Tw(:,:) = p(:,1,:) / (R * rho(:,1,:))
      call print_turbulent_boundary_layer(step,nx,ny,nz,dy,y,Tw,u,rho)
      write(filename, "(a, i5.5, a)") "data/Q",int(step+step_offset),".vtr"
    endif
    open(10,file=filename,status="replace",action="write",form="unformatted",access="stream",convert="Little_ENDIAN")
    call print_xml(nx,ny,nz,3,real(x),real(y),real(z),real(rho1d),real(p1d),real(T1d),real(M1d),real(v1d),real(Qcriterion1d))
   
    if (step == np) then
      rhom(:) = rhom(:) / dble(np)
      pm(:)   = pm(:)   / dble(np)
      Tm(:)   = Tm(:)   / dble(np)
      Mm(:)   = Mm(:)   / dble(np)
      vm(:)   = vm(:)   / dble(np)
      write(filename, "(a)") "data/Qmean.vtr"
      open(10,file=filename,status="replace",action="write",form="unformatted",access="stream",convert="Little_ENDIAN")
      call print_xml(nx,ny,nz,3,real(x),real(y),real(z),real(rhom),real(pm),real(Tm),real(Mm),real(vm))
    endif

    deallocate(rho,u,v,w,p,div,omega,Qcriterion,Tw,rho1d,p1d,T1d,M1d,v1d,div1d,omega1d,Qcriterion1d)
  end subroutine print_vtk_3D
end module print

