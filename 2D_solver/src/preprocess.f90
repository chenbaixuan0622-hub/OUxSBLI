module preprocess
  use cudafor
  use cufft
  use mpi
  use print
  use calc_forcing, only : plan_fwd, plan_inv
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


  subroutine allocate_device_mem(myrank, nx, ny, x, y, dx, dy, xix, etay, Jacobian, ruvp, T, mu, E, F, Fout)
    use mod_globals, only : id_accuracy, id_visc, id_force, kmax
    integer, intent(in)                       :: myrank, nx, ny
    real(8), intent(out), allocatable, device :: x(:), y(:), dx(:), dy(:), xix(:), etay(:), Jacobian(:,:)
    real(8), intent(out), allocatable, device :: ruvp(:,:,:), T(:,:), mu(:,:)
    real(8), intent(out), allocatable, device :: E(:,:,:), F(:,:,:), Fout(:,:,:)
    integer ierr, istat, accuracy
    allocate(ruvp(4,nx,ny), E(4,nx-1,ny-2), F(4,nx-2,ny-1), stat=ierr)
    allocate(x(nx), y(ny), dx(nx-1), dy(ny-1), xix(nx-1), etay(ny-1), Jacobian(nx,ny), stat=ierr)
    if (kind(id_visc) == 2) then
      allocate(T(nx,ny), mu(1,1), stat=ierr)
    elseif (kind(id_visc) == 4) then
      allocate(T(nx,ny), mu(nx,ny), stat=ierr)
    endif
    if (kind(id_force) == 4) then
      allocate(Fout(3,nx-2,ny-2), stat=ierr)
      if (kind(id_accuracy) == 8) then
        accuracy = 6
      elseif (kind(id_accuracy) == 4) then
        accuracy = 4
      elseif (kind(id_accuracy) == 2) then
        accuracy = 2
      endif
      istat = cufftPlan2d(plan_fwd, nx-accuracy, ny-accuracy, CUFFT_D2Z)
      istat = cufftPlan2d(plan_inv, nx-accuracy, ny-accuracy, CUFFT_Z2D)
    else
      allocate(Fout(1,1,1), stat=ierr)
    endif
    if (ierr /= 0) then
      print *, "myrank is ", myrank, " memory allocation failed", ierr
    else
      print *, "myrank is ", myrank, " memory allocation has completed"
    endif
  end subroutine allocate_device_mem


  subroutine pre_calc(nx, ny, myrank, nranks, x_cpu, dx_cpu, y_cpu, dy_cpu, Jacobian_cpu, Q, overlap, &
                      x, y, dx, dy, xix, etay, Jacobian, QJ)
    integer, intent(in)    :: nx, ny, myrank, nranks
    real(8), intent(in)    :: x_cpu(nx), dx_cpu(nx-1), y_cpu(ny), dy_cpu(ny-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(4,nx,ny)
    integer, intent(out)   :: overlap
    real(8), intent(out), device :: x(nx), y(ny), dx(nx-1), dy(ny-1), xix(nx-1), etay(ny-1), Jacobian(nx,ny)
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
    x        = x_cpu
    y        = y_cpu
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
    call print_vtk(0, nx, ny, x_cpu, y_cpu, rho1d, p1d, v1d)
  end subroutine pre_calc
end module preprocess

