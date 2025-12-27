module preprocess
  use cudafor
  use mpi
  use print
  implicit none
contains
  subroutine check_gpu(mygpu)
    integer, intent(in) :: mygpu
    integer ilen, stat
    type(cudaDeviceProp) prop
    stat = cudaSetDevice(mygpu)
    stat = cudaGetDeviceProperties(prop, mygpu)
    ilen = verify(prop%name, ' ', .true.)
    print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", mygpu, ") is available"
  end subroutine check_gpu


  subroutine allocate_device_mem(myrank, nx, ny, dx, dy, xix, etay, Jacobian, ruvp, T, mu, E, F)
    use mod_globals, only : id_visc
    integer, intent(in)                       :: myrank, nx, ny
    real(8), intent(out), allocatable, device :: dx(:), dy(:), xix(:), etay(:), Jacobian(:,:)
    real(8), intent(out), allocatable, device :: ruvp(:,:,:), T(:,:), mu(:,:)
    real(8), intent(out), allocatable, device :: E(:,:,:), F(:,:,:)
    integer ierr
    allocate(ruvp(4,nx,ny), E(4,nx-1,ny-2), F(4,nx-2,ny-1), stat=ierr)
    allocate(dx(nx-1), dy(ny-1), xix(nx-1), etay(ny-1), Jacobian(nx,ny), stat=ierr)
    if (kind(id_visc) == 2) then
      allocate(T(nx,ny), mu(1,1), stat=ierr)
    elseif (kind(id_visc) == 4) then
      allocate(T(nx,ny), mu(nx,ny), stat=ierr)
    endif
    if (ierr /= 0) then
      print *, "myrank is ", myrank, " memory allocation failed", ierr
    else
      print *, "myrank is ", myrank, " memory allocation has completed"
    endif
  end subroutine allocate_device_mem


  subroutine pre_calc(nx, ny, myrank, nranks, x, dx_cpu, y, dy_cpu, Jacobian_cpu, Q, overlap, &
                      dx, dy, xix, etay, Jacobian, QJ)
    integer, intent(in)    :: nx, ny, myrank, nranks
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1), y(ny), dy_cpu(ny-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(4,nx,ny)
    integer, intent(out)   :: overlap
    real(8), intent(out), device :: dx(nx-1), dy(ny-1), xix(nx-1), etay(ny-1), Jacobian(nx,ny)
    real(8), intent(out), device :: QJ(4,nx,ny)
    real(8) xix_cpu(nx-1), etay_cpu(ny-1)
    real(4) rho1d(nx*ny), p1d(nx*ny), v1d(nx*ny*3)
    integer i, j, l, ierr
    ! set Q / Jacobian
    do j = 1, ny
      do i = 1, nx
        do l = 1, 4
          Q(l,i,j) = Q(l,i,j) / Jacobian_cpu(i,j)
    enddo;enddo;enddo
    ! copy on GPU
    xix_cpu  = 1.d0 / dx_cpu
    etay_cpu = 1.d0 / dy_cpu
    dx       = dx_cpu
    dy       = dy_cpu
    xix      = xix_cpu
    etay     = etay_cpu
    Jacobian = Jacobian_cpu
    QJ = Q
    ! for multi GPU
    if (kind(id_accuracy) == 8) then
      overlap = 3
    elseif (kind(id_accuracy) == 4) then
      overlap = 2
    else
      overlap = 1
    endif
    call make_1d_for_print(nx, ny, Jacobian_cpu, Q, rho1d, p1d, v1d)
    call print_vtk(0, nx, ny, x, y, rho1d, p1d, v1d)
  end subroutine pre_calc
end module preprocess

