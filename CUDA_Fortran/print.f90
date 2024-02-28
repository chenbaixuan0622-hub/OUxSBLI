module print
  implicit none
contains
  subroutine print_vtk(step,nx,ny,nz,dx,dy,dz,gamma,Q)
    integer, intent(in) :: step, nx, ny, nz
    real(8), intent(in) :: dx, dy, dz, gamma, Q(nx,ny,nz,5)
    integer i, j, k
    real(8), dimension(nx,ny,nz) :: rho, u, v, w!, p
    character(len=40) filename
    rho = Q(:,:,:,1)
    u = Q(:,:,:,2) / rho
    v = Q(:,:,:,3) / rho
    w = Q(:,:,:,4) / rho
    !p = (gamma - 1.d0) * (Q(:,:,:,5) - 0.5d0 * rho * (u ** 2 + v ** 2 + w ** 2))
    
    write(filename, "(a, i5.5,a)") "data/u",int(step),".vtk"
    open(10,file=filename)
    write(10,"('# vtk DataFile Version 3.0')")
    write(10,"('u')")
    write(10,"('ASCII')")
    write(10,"('DATASET STRUCTURED_GRID')")
    write(10,"('DIMENSIONS',3(1x,i4))") nx, ny, nz
    write(10,"('POINTS',i9,' float')") nx * ny * nz
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          write(10,"(3(f9.4,1x))") (i-1)*dx, (j-1)*dy, (k-1)*dz
        enddo
      enddo
    enddo

    write(10,"('POINT_DATA',i9)") nx * ny * nz
    write(10,"('SCALARS u float')")
    write(10,"('LOOKUP_TABLE default')")
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          write(10,"(3(f9.4,1x))") u(i,j,k)
        enddo
      enddo
    enddo
    close(10)

    write(filename, "(a, i5.5,a)") "data/v",int(step),".vtk"
    open(11,file=filename)
    write(11,"('# vtk DataFile Version 3.0')")
    write(11,"('v')")
    write(11,"('ASCII')")
    write(11,"('DATASET STRUCTURED_GRID')")
    write(11,"('DIMENSIONS',3(1x,i4))") nx, ny, nz
    write(11,"('POINTS',i9,' float')") nx * ny * nz
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          write(11,"(3(f9.4,1x))") (i-1)*dx, (j-1)*dy, (k-1)*dz
        enddo
      enddo
    enddo

    write(11,"('POINT_DATA',i9)") nx * ny * nz
    write(11,"('SCALARS v float')")
    write(11,"('LOOKUP_TABLE default')")
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          write(11,"(3(f9.4,1x))") v(i,j,k)
        enddo
      enddo
    enddo
    close(11)

    write(filename, "(a, i5.5,a)") "data/w",int(step),".vtk"
    open(12,file=filename)
    write(12,"('# vtk DataFile Version 3.0')")
    write(12,"('w')")
    write(12,"('ASCII')")
    write(12,"('DATASET STRUCTURED_GRID')")
    write(12,"('DIMENSIONS',3(1x,i4))") nx, ny, nz
    write(12,"('POINTS',i9,' float')") nx * ny * nz
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          write(12,"(3(f9.4,1x))") (i-1)*dx, (j-1)*dy, (k-1)*dz
        enddo
      enddo
    enddo

    write(12,"('POINT_DATA',i9)") nx * ny * nz
    write(12,"('SCALARS w float')")
    write(12,"('LOOKUP_TABLE default')")
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          write(12,"(3(f9.4,1x))") w(i,j,k)
        enddo
      enddo
    enddo
    close(12)
  end subroutine print_vtk

  subroutine print_message(nx, ny, nz, dx, dy, dz, dt, U, mu, gamma, kappa, Cp, elapsedTime)
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in) :: dt, dx, dy, dz, mu, U, gamma, kappa, Cp, elapsedTime
    real(8) Lx, Ly, Lz, endTime
    Lx = nx * dx
    Ly = ny * dy
    Lz = nz * dz
    open(13,file="output_message.txt")
    write(13,"('mesh information')")
    write(13,"('Lx =',i9,'was devided by',i4,'dx =',i9)") Lx, nx, dx
    write(13,"('Ly =',i9,'was devided by',i4,'dy =',i9)") Ly, ny, dy
    write(13,"('Lz =',i9,'was devided by',i4,'dz =',i9)") Lz, nz, dz
    write(13,"('\n')")

    write(13,"('simulation time')")
    write(13,"('start time     =',i9)") 0.d0
    write(13,"('end time       =',i9)") endTime
    write(13,"('dt             =',i9)") dt
    write(13,"('Courant number =',i9)") U * dt / dx
    write(13,"('\n')")

    write(13,"('Physical properties')")
    write(13,"('gamma =',i9,' [Pa s]')") gamma
    write(13,"('mu    =',i9,' [Pa s]')") mu
    write(13,"('kappa =',i9,' [W/(m K)]')") kappa
    write(13,"('Cp    =',i9,' [J/(kg K)]')") Cp
    write(13,"('\n')")

    write(13,"('elapsed time =',i9)") elapsedTime
    close(13)
  end subroutine print_message
end module print
