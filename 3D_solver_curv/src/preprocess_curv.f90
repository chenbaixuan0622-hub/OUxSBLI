!> Preprocessing and GPU memory management for curvilinear case.
!> Contains only the curvilinear routines; rectangular-grid routines
!> (allocate_device_mem, pre_calc, pre_rescale) are omitted to avoid
!> pulling in calc_flux_base.mod which is not compiled for this case.
module preprocess_curv
  use cudafor
  use mpi
  use print
  use calc_flux_base_curv, only : init_sensor
  implicit none
contains
  !> Check GPU device properties and availability
  !> Prints device name and capability information
  subroutine check_gpu(mygpu)
    integer, intent(in) :: mygpu                         !< GPU device ID to check
    integer ilen, stat
    type(cudaDeviceProp) prop
    stat = cudaSetDevice(mygpu)
    stat = cudaGetDeviceProperties(prop, mygpu)
    ilen = verify(prop%name, ' ', .true.)
    print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", mygpu, ") is available"
  end subroutine check_gpu


  !> Allocate GPU device memory for curvilinear simulation
  !> New arrays: 2D inverse metrics, 2D face normals, scalar dt_Szeta for z-direction
  subroutine allocate_device_mem_curv(myrank, nx, ny, nz, &
      dt_Szeta, n_xi_x, n_xi_y, n_eta_x, n_eta_y, &
      xi_x, xi_y, eta_x, eta_y, Jacobian, &
      ruvwp, T, mu, E, F, G)
    use mod_constant, only : id_visc
    use calc_flux_base_curv, only : init_sensor_curv
    integer, intent(in)                       :: myrank
    integer, intent(in)                       :: nx
    integer, intent(in)                       :: ny
    integer, intent(in)                       :: nz
    real(8), intent(out), allocatable, device :: dt_Szeta(:,:)
    real(8), intent(out), allocatable, device :: n_xi_x(:,:), n_xi_y(:,:)
    real(8), intent(out), allocatable, device :: n_eta_x(:,:), n_eta_y(:,:)
    real(8), intent(out), allocatable, device :: xi_x(:,:), xi_y(:,:)
    real(8), intent(out), allocatable, device :: eta_x(:,:), eta_y(:,:)
    real(8), intent(out), allocatable, device :: Jacobian(:,:)
    real(8), intent(out), allocatable, device :: ruvwp(:,:,:,:)
    real(8), intent(out), allocatable, device :: T(:,:,:)
    real(8), intent(out), allocatable, device :: mu(:,:,:)
    real(8), intent(out), allocatable, device :: E(:,:,:,:)
    real(8), intent(out), allocatable, device :: F(:,:,:,:)
    real(8), intent(out), allocatable, device :: G(:,:,:,:)
    integer ierr
    allocate(ruvwp(nx,5,ny,nz), E(5,nx-1,ny-2,nz-2), F(5,nx-2,ny-1,nz-2), G(5,nx-2,ny-2,nz-1), stat=ierr)
    allocate(dt_Szeta(nx-2,ny-2), stat=ierr)
    allocate(n_xi_x(nx-1,ny-2), n_xi_y(nx-1,ny-2), stat=ierr)
    allocate(n_eta_x(nx-2,ny-1), n_eta_y(nx-2,ny-1), stat=ierr)
    allocate(xi_x(nx,ny), xi_y(nx,ny), stat=ierr)
    allocate(eta_x(nx,ny), eta_y(nx,ny), stat=ierr)
    allocate(Jacobian(nx,ny), stat=ierr)

    if (kind(id_visc) == 2) then
      allocate(T(nx,ny,nz), mu(1,1,1), stat=ierr)
    elseif (kind(id_visc) == 4) then
      allocate(T(nx,ny,nz), mu(nx,ny,nz), stat=ierr)
    endif

    call init_sensor_curv(nx, ny, nz)
    if (ierr /= 0) then
      print *, "myrank is ", myrank, " memory allocation failed", ierr
    else
      print *, "myrank is ", myrank, " memory allocation has completed (curvilinear)"
    endif
  end subroutine allocate_device_mem_curv


  !> Preprocessing for curvilinear grids: compute QJ, dt_Szeta, copy metrics
  !> QJ = Q / Jacobian_code where Jacobian_code = 1/(J_2D*dz)
  subroutine pre_calc_curv(nx, ny, nz, myrank, nranks, dz, &
      Jac_cpu, n_xi_x_cpu, n_xi_y_cpu, n_eta_x_cpu, n_eta_y_cpu, &
      xi_x_cpu, xi_y_cpu, eta_x_cpu, eta_y_cpu, &
      x_phys, y_phys, z, Q, overlap, &
      dt_Szeta, n_xi_x, n_xi_y, n_eta_x, n_eta_y, &
      xi_x, xi_y, eta_x, eta_y, Jacobian, QJ, ke0, entropy0)
    use mod_globals, only : dt
    use mod_constant, only : id_accuracy
    use print_curv, only : print_vtk_curv
    integer, intent(in)                      :: nx
    integer, intent(in)                      :: ny
    integer, intent(in)                      :: nz
    integer, intent(in)                      :: myrank
    integer, intent(in)                      :: nranks
    real(8), intent(in)                      :: dz
    real(8), intent(in)                      :: Jac_cpu(nx,ny)
    real(8), intent(in)                      :: n_xi_x_cpu(nx-1,ny-2), n_xi_y_cpu(nx-1,ny-2)
    real(8), intent(in)                      :: n_eta_x_cpu(nx-2,ny-1), n_eta_y_cpu(nx-2,ny-1)
    real(8), intent(in)                      :: xi_x_cpu(nx,ny), xi_y_cpu(nx,ny)
    real(8), intent(in)                      :: eta_x_cpu(nx,ny), eta_y_cpu(nx,ny)
    real(8), intent(in)                      :: x_phys(nx,ny), y_phys(nx,ny), z(nz)
    real(8), intent(inout)                   :: Q(nx,5,ny,nz)
    integer, intent(out)                     :: overlap
    real(8), intent(out), device, contiguous :: dt_Szeta(nx-2,ny-2)
    real(8), intent(out), device, contiguous :: n_xi_x(nx-1,ny-2), n_xi_y(nx-1,ny-2)
    real(8), intent(out), device, contiguous :: n_eta_x(nx-2,ny-1), n_eta_y(nx-2,ny-1)
    real(8), intent(out), device, contiguous :: xi_x(nx,ny), xi_y(nx,ny)
    real(8), intent(out), device, contiguous :: eta_x(nx,ny), eta_y(nx,ny)
    real(8), intent(out), device, contiguous :: Jacobian(nx,ny)
    real(8), intent(out), device, contiguous :: QJ(nx,5,ny,nz)
    real(4), intent(inout)                   :: ke0, entropy0
    real(8) dt_Szeta_cpu(nx-2,ny-2)
    real(4) rho1d(nx*ny*nz), p1d(nx*ny*nz), v1d(nx*ny*nz*3)
    integer i, j, k, l, ierr
    ! 1. Compute QJ = Q / Jacobian_code and dt_Szeta = dt / (Jacobian_code * dz) = dt * J_2D
    do k = 1, nz
      do j = 1, ny
        do l = 1, 5
          do i = 1, nx
            Q(i,l,j,k) = Q(i,l,j,k) / Jac_cpu(i,j)
    enddo;enddo;enddo;enddo
    ! dt_Szeta(i,j) = dt * J_2D(i+1,j+1) where Jac_cpu(i,j) = 1/(J_2D(i,j)*dz)
    do j = 1, ny-2
      do i = 1, nx-2
        dt_Szeta_cpu(i,j) = dt / (Jac_cpu(i+1,j+1) * dz)
    enddo;enddo
    ! 2. Copy metrics to device
    dt_Szeta = dt_Szeta_cpu
    n_xi_x   = n_xi_x_cpu
    n_xi_y   = n_xi_y_cpu
    n_eta_x  = n_eta_x_cpu
    n_eta_y  = n_eta_y_cpu
    xi_x     = xi_x_cpu
    xi_y     = xi_y_cpu
    eta_x    = eta_x_cpu
    eta_y    = eta_y_cpu
    Jacobian = Jac_cpu
    QJ       = Q
    ! 3. Set overlap based on accuracy
    if (kind(id_accuracy) == 8) then
      overlap = 3
    elseif (kind(id_accuracy) == 4) then
      overlap = 2
    else
      overlap = 1
    endif
    ! 4. Output step 0
    call make_1d_for_print(nx, ny, nz, Jac_cpu, Q, rho1d, p1d, v1d)
    call print_vtk_curv(0, nx, ny, nz, myrank+1, nranks, x_phys, y_phys, z, rho1d, p1d, v1d, ke0, entropy0)
    call MPI_SEND(ke0,      1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
    call MPI_SEND(entropy0, 1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
  end subroutine pre_calc_curv
end module preprocess_curv

