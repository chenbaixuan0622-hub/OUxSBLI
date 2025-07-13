module print
  use mpi
  use mod_globals, only : id_accuracy, nt, np, dt, step_offset, gamma, R, Lx
  implicit none
  interface
    subroutine print_entropy(step, nx, ny, nz, rho1d, p1d, entropy0, myrank)
      integer, intent(in)                      :: step, nx, ny, nz
      real(4), intent(in), dimension(nx*ny*nz) :: rho1d, p1d
      real(4), intent(inout)                   :: entropy0
      integer, intent(in), optional            :: myrank
    end subroutine print_entropy
  end interface
  
  interface
    subroutine print_KE(step, nx, ny, nz, rho1d, v1d, ke0, myrank)
      integer, intent(in)           :: step, nx, ny, nz
      real(4), intent(in)           :: rho1d(nx*ny*nz), v1d(nx*ny*nz*3)
      real(4), intent(inout)        :: ke0
      integer, intent(in), optional :: myrank
    end subroutine print_KE
  end interface

  interface print_vtk
    module procedure print_vtk_2D, print_vtk_3D
  end interface
contains
  subroutine print_entropy(step, nx, ny, nz, rho1d, p1d, entropy0, myrank)
    integer, intent(in)                      :: step, nx, ny, nz
    real(4), intent(in), dimension(nx*ny*nz) :: rho1d, p1d
    real(4), intent(inout)                   :: entropy0
    integer, intent(in), optional            :: myrank
    real(4) entropy, t
    character(len=40) filename
    integer i, j, k, l, accuracy, offset
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
    entropy = 0.e0
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          l = i + (j-1) * nx + (k-1) * nx * ny
          entropy = entropy + rho1d(l) * log(p1d(l) * (rho1d(l)**(-gamma)))
    enddo;enddo;enddo
    entropy = entropy / dble((nx-accuracy) * (ny-accuracy) * (nz-accuracy))
    if (step == 0) then
      entropy0 = entropy
    endif
    t = nt * step * dt
    if (present(myrank)) then
      write(filename, "(a, i0, a)") "data/",int(myrank),"/entropy.d"
      open(10,file=filename, position="append")
    else
      open(10,file="data/entropy.d", position="append")
    endif
    write(10,"(2e12.4)") t, (entropy0 - entropy) / entropy0
    close(10)
  end subroutine print_entropy

  subroutine print_KE(step, nx, ny, nz, rho1d, v1d, ke0, myrank)
    integer, intent(in)           :: step, nx, ny, nz
    real(4), intent(in)           :: rho1d(nx*ny*nz), v1d(nx*ny*nz*3)
    real(4), intent(inout)        :: ke0
    integer, intent(in), optional :: myrank
    real(4) ke, t
    character(len=40) filename
    integer i, j, k, l, m, accuracy, offset
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
    ke = 0.e0
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          l  = i + (j-1) * nx + (k-1) * nx * ny
          m  = 1 + (i-1) * 3 + (j-1) * 3 * nx + (k-1) * 3 * nx * ny 
          ke = ke + 0.5d0 * rho1d(l) * (v1d(m)**2 + v1d(m+1)**2 + v1d(m+2)**2)
    enddo;enddo;enddo
    ke = ke / dble((nx-accuracy) * (ny-accuracy) * (nz-accuracy))
    if (step == 0) then
      ke0 = ke
    endif
    t = nt * step * dt
    if (present(myrank)) then
      write(filename, "(a, i0, a)") "data/",int(myrank),"/kinetic_energy.d"
      open(10,file=filename, position="append")
    else
      open(10,file="data/kinetic_energy.d", position="append")
    endif
    write(10,"(3e12.4)") t, ke, ke / ke0
    close(10)
  end subroutine print_KE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_1d(step,nx,ny,nz,x,y,z,Jacobian,QJ)
    integer, intent(in) :: step, nx, ny, nz
    real(4), intent(in) :: x(nx), y(ny), z(nz)
    real(8), intent(in) :: Jacobian(ny), QJ(nx,ny,nz,5)
    real(8) rho, u, v, w, p, Lx
    integer i, nyh, nzh
    character(len=40) filename
    nyh = int(0.5 * ny)
    nzh = int(0.5 * nz)
    Lx  = x(nx)
    write(filename, "(a, i5.5, a)") "data/1d/Q", int(step), ".d"
    open(10,file=filename)
    do i = 1, nx
      rho = Jacobian(nyh) * QJ(i,nyh,nzh,1)
      u   = QJ(i,nyh,nzh,2) / QJ(i,nyh,nzh,1)
      v   = QJ(i,nyh,nzh,3) / QJ(i,nyh,nzh,1)
      w   = QJ(i,nyh,nzh,4) / QJ(i,nyh,nzh,1)
      p   = (gamma - 1.d0) * (Jacobian(nyh) * QJ(i,nyh,nzh,5) - 0.5d0 * rho * (u**2 + v**2 + w**2))
      write(10,"(7e12.4)") x(i) / Lx, rho, u, p
    enddo
    close(10)
  end subroutine print_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_xml(ni, nj, nk, dimension, x, y, z, rho1d, p1d, v1d)
    integer, intent(in)                                :: ni, nj, nk, dimension
    real(4), intent(in)                                :: x(ni), y(nj), z(nk)
    real(4), intent(in), dimension(ni*nj*nk)           :: rho1d, p1d
    real(4), intent(in), dimension(dimension*ni*nj*nk) :: v1d
    integer(4) byte_x, byte_y, byte_z, byte_rho, byte_p, byte_v
    character :: lf*1, str1*4, str2*4, str3*4, str4*1
    character :: offset1*12, offset2*12, offset3*12, offset4*12, offset5*12
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
    byte_v   = 4 + 4 * (dimension * ni * nj * nk)
    write(offset1(1:12),'(i12)') int(byte_x, kind=8)
    write(offset2(1:12),'(i12)') int(byte_x, kind=8) + int(byte_y, kind=8)
    write(offset3(1:12),'(i12)') int(byte_x, kind=8) + int(byte_y, kind=8) + int(byte_z, kind=8)
    write(offset4(1:12),'(i12)') int(byte_x, kind=8) + int(byte_y, kind=8) + int(byte_z, kind=8) + &
                                 int(byte_rho, kind=8)
    write(offset5(1:12),'(i12)') int(byte_x, kind=8) + int(byte_y, kind=8) + int(byte_z, kind=8) + &
                                 int(byte_rho, kind=8) + int(byte_p, kind=8)
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
    & Name="velocity" NumberOfComponents="'//str4//'"/>'//lf
    write(10) '      </PointData>'//lf
    write(10) '    </Piece>'//lf
    write(10) '  </RectilinearGrid>'//lf
    write(10) '  <AppendedData encoding="raw">'//lf
    write(10) '  _', byte_x, x, byte_y, y, byte_z, z, byte_rho, rho1d, byte_p, p1d, byte_v, v1d, lf
    write(10) '  </AppendedData>'//lf
    write(10) '</VTKFile>'//lf
    close(10)
  end subroutine print_xml

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_1d_for_print(nx, ny, nz, Jacobian, QJ, rho1d, p1d, v1d)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: Jacobian(ny), QJ(nx,ny,nz,5) ! Q / Jacobian
    real(4), intent(out), dimension(nx*ny*nz)   :: rho1d, p1d
    real(4), intent(out), dimension(nx*ny*nz*3) :: v1d
    real(8) rho, u, v, w, p
    integer i, j, k, l, m
    l = 1
    m = 1
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rho      = Jacobian(j) * QJ(i,j,k,1)
          u        = QJ(i,j,k,2) / QJ(i,j,k,1)
          v        = QJ(i,j,k,3) / QJ(i,j,k,1)
          w        = QJ(i,j,k,4) / QJ(i,j,k,1)
          p        = (gamma - 1.d0) * (Jacobian(j) * QJ(i,j,k,5) - 0.5d0 * rho * (u**2 + v**2 + w**2))
          rho1d(l) = real(rho)
          p1d(l)   = real(p)
          v1d(m)   = real(u)
          v1d(m+1) = real(v)
          v1d(m+2) = real(w)
          l = l + 1
          m = m + 3
    enddo;enddo;enddo
  end subroutine make_1d_for_print

  subroutine print_vtk_2D(step, nx, ny, x, y, rho1d, p1d, v1d)
    integer, intent(in) :: step, nx, ny
    real(8), intent(in) :: x(nx), y(ny)
    real(4), intent(in) :: rho1d(nx*ny), p1d(nx*ny), v1d(nx*ny*2)
    real(8) :: z(1) = 0.d0
    character(len=40) filename
    write(filename, "(a, i5.5,a)") "data/Q",int(step),".vtr"
    open(10,file=filename,status="replace",action="write",form="unformatted",access="stream",convert="Little_ENDIAN")
    call print_xml(nx, ny, 1, 2, real(x), real(y), real(z), rho1d, p1d, v1d)
  end subroutine print_vtk_2D
  
  subroutine print_vtk_3D(step, nx, ny, nz, myrank, nranks, x, y, z, rho1d, p1d, v1d, ke0, entropy0)
    integer, intent(in)    :: step, nx, ny, nz, myrank, nranks
    real(8), intent(in)    :: x(nx), y(ny), z(nz)
    real(4), intent(in)    :: rho1d(nx*ny*nz), p1d(nx*ny*nz), v1d(nx*ny*nz*3)
    real(4), intent(inout) :: ke0, entropy0
    character(len=40) filename
    if (nranks >= 4) then
      call print_entropy(step, nx, ny, nz, rho1d, p1d, entropy0, myrank)
      call print_KE(step, nx, ny, nz, rho1d, v1d, ke0, myrank)
      write(filename, "(a, i0, a, i5.5, a)") "data/",int(myrank),"/Q",int(step+step_offset),".vtr"
    else
      call print_entropy(step, nx, ny, nz, rho1d, p1d, entropy0)
      call print_KE(step, nx, ny, nz, rho1d, v1d, ke0)
      write(filename, "(a, i5.5, a)") "data/Q",int(step+step_offset),".vtr"
    endif
    open(10,file=filename,status="replace",action="write",form="unformatted",access="stream",convert="Little_ENDIAN")
    call print_xml(nx, ny, nz, 3, real(x), real(y), real(z), rho1d, p1d, v1d)
  end subroutine print_vtk_3D
end module print

