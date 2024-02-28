module print
  implicit none
contains
  subroutine print_vtk(step,nx,ny,nz,dx,dy,dz,gamma,Q)
    integer, intent(in) :: step, nx, ny, nz
    real(8), intent(in) :: dx, dy, dz, gamma, Q(nx,ny,nz,5)
    integer i, j, k
    real(8), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8), dimension(nx,ny-1,nz-1) :: wx
    character(len=40) filename
    rho = Q(:,:,:,1)
    !u = Q(:,:,:,2) / rho
    v = Q(:,:,:,3) / rho
    w = Q(:,:,:,4) / rho
    !p = (gamma - 1.d0) * (Q(:,:,:,5) - 0.5d0 * rho * (u ** 2 + v ** 2 + w ** 2))
    
    write(filename, "(a, i5.5,a)") "data/u",int(step),".vtk"
    open(10,file=filename)
    write(10,"('# vtk DataFile Version 3.0')")
    write(10,"('u')")
    write(10,"('ASCII')")
    write(10,"('DATASET STRUCTURED_GRID')")
    write(10,"('DIMENSIONS',3(1x,i3))") nx, ny, nz
    write(10,"('POINTS',i9,' float')") nx * ny * nz
    do k = 1, nz-1
      do j = 1, ny-1
        do i = 1, nx
          write(10,"(3(f9.4,1x))") (i-1)*dx, (j-1)*dy, (k-1)*dz
        enddo
      enddo
    enddo

    write(10,"('POINT_DATA',i9)") nx * ny * nz
    write(10,"('SCALARS u float')")
    write(10,"('LOOKUP_TABLE default')")
    do k = 1, nz-1
      do j = 1, ny-1
        do i = 1, nx
          write(10,"(3(f9.4,1x))") (-w(i,j,k) + w(i,j+1,k)) / dy - (-v(i,j,k) + v(i,j,k+1)) / dz
        enddo
      enddo
    enddo
    close(10)
  end subroutine print_vtk
end module print
