module calc_igr
  use mod_globals, only : blocks, threads
  implicit none
contains
  subroutine set_init_sigma(nx, ny, nz, dx, dy, dz, Jacobian, Q, sigma_gpu)
    integer, intent(in)          :: nx, ny, nz
    real(8), intent(in)          :: dx(nx-1), dy(ny-1), dz(nz-1), Jacobian(nx,ny), Q(5,nx,ny,nz)
    real(4), intent(out), device :: sigma_gpu(nx,ny,nz)
    real(4), allocatable :: Ab(:,:,:,:), sigma(:,:,:), over_rho(:,:,:), u(:,:,:), v(:,:,:), w(:,:,:)
    real(4) :: alpha, dudx, dvdy, dwdz, RHS, a0, alpha_a0, omega = 1.6e0
    real(4) over_rhox1, over_rhox2, over_rhoy1, over_rhoy2, over_rhoz1, over_rhoz2, sigma_new
    integer :: i, j, k, itr, maxitr = 100
    ! set coef
    allocate(Ab(7,nx-2,ny-2,nz-2), sigma(nx,ny,nz), over_rho(nx,ny,nz), u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz))
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          over_rho(i,j,k) = 1.e0 / real(Q(1,i,j,k) * Jacobian(i,j))
          u(i,j,k) = real(Q(2,i,j,k) / Q(1,i,j,k))
          v(i,j,k) = real(Q(3,i,j,k) / Q(1,i,j,k))
          w(i,j,k) = real(Q(4,i,j,k) / Q(1,i,j,k))
    enddo;enddo;enddo
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          alpha = 0.5e0 * dx(i)**2
          dudx = 0.5e0 * (-u(i-1,j,k) + u(i+1,j,k)) / dx(i)
          dvdy = 0.5e0 * (-v(i,j-1,k) + v(i,j+1,k)) / dy(j)
          dwdz = 0.5e0 * (-w(i,j,k-1) + w(i,j,k+1)) / dz(k)
          RHS = (dudx + dvdy + dwdz)**2 + dudx**2 + dvdy**2 + dwdz**2
          over_rhox1 = 0.5e0 * (over_rho(i-1,j,k) + over_rho(i,j,k))
          over_rhox2 = 0.5e0 * (over_rho(i,j,k)   + over_rho(i+1,j,k))
          over_rhoy1 = 0.5e0 * (over_rho(i,j-1,k) + over_rho(i,j,k))
          over_rhoy2 = 0.5e0 * (over_rho(i,j,k)   + over_rho(i,j+1,k))
          over_rhoz1 = 0.5e0 * (over_rho(i,j,k-1) + over_rho(i,j,k))
          over_rhoz2 = 0.5e0 * (over_rho(i,j,k)   + over_rho(i,j,k+1))
          a0 = over_rho(i,j,k) + alpha * ((over_rhox1 + over_rhox2) / dx(i)**2 &
                                        + (over_rhoy1 + over_rhoy2) / dy(j)**2 &
                                        + (over_rhoz1 + over_rhoz2) / dz(k)**2)
          alpha_a0 = alpha / a0
          Ab(1,i-1,j-1,k-1) = over_rhox1 / dx(i)**2 * alpha_a0
          Ab(2,i-1,j-1,k-1) = over_rhox2 / dx(i)**2 * alpha_a0
          Ab(3,i-1,j-1,k-1) = over_rhoy1 / dy(j)**2 * alpha_a0
          Ab(4,i-1,j-1,k-1) = over_rhoy2 / dy(j)**2 * alpha_a0
          Ab(5,i-1,j-1,k-1) = over_rhoz1 / dz(k)**2 * alpha_a0
          Ab(6,i-1,j-1,k-1) = over_rhoz2 / dz(k)**2 * alpha_a0
          Ab(7,i-1,j-1,k-1) = RHS * alpha_a0
    enddo;enddo;enddo
    deallocate(over_rho, u, v, w)
    do itr = 1, maxitr
      do k = 2, nz-1
        do j = 2, ny-1
          do i = 2, nx-1
            sigma_new = Ab(1,i-1,j-1,k-1) * sigma(i-1,j,k) + Ab(2,i-1,j-1,k-1) * sigma(i+1,j,k) &
                      + Ab(3,i-1,j-1,k-1) * sigma(i,j-1,k) + Ab(4,i-1,j-1,k-1) * sigma(i,j+1,k) &
                      + Ab(5,i-1,j-1,k-1) * sigma(i,j,k-1) + Ab(6,i-1,j-1,k-1) * sigma(i,j,k+1) + Ab(7,i-1,j-1,k-1)
            sigma(i,j,k) = sigma(i,j,k) + omega * (sigma_new - sigma(i,j,k))
      enddo;enddo;enddo
      ! set bc sigma
      sigma(1,:,:)  = sigma(2,:,:)
      sigma(nx,:,:) = sigma(nx-1,:,:)
      sigma(:,1,:)  = sigma(:,2,:)
      sigma(:,ny,:) = sigma(:,ny-1,:)
      sigma(:,:,1)  = sigma(:,:,2)
      sigma(:,:,nz) = sigma(:,:,nz-1)
    enddo
    sigma_gpu(:,:,:) = sigma(:,:,:)
    deallocate(Ab, sigma)
  end subroutine set_init_sigma


  subroutine set_bc_sigma(nx, ny, nz, sigma)
    integer, intent(in), value     :: nx, ny, nz
    real(4), intent(inout), device :: sigma(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        sigma(1,j,k)  = sigma(2,j,k)
        sigma(nx,j,k) = sigma(nx-1,j,k)
    enddo;enddo
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        sigma(i,1,k)  = sigma(i,2,k)
        sigma(i,ny,k) = sigma(i,ny-1,k)
    enddo;enddo
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        sigma(i,j,1)  = sigma(i,j,2)
        sigma(i,j,nz) = sigma(i,j,nz-1)
    enddo;enddo
  end subroutine set_bc_sigma


  attributes(global) subroutine set_a(nx, ny, nz, over_dx, over_dy, over_dz, over_rho, ruvwp, sensor, Ab)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: over_dx(nx-1), over_dy(ny-1), over_dz(nz-1)
    real(4), intent(in), device  :: over_rho(nx,ny,nz)
    real(8), intent(in), device  :: ruvwp(5,nx,ny,nz), sensor(nx,ny,nz)
    real(4), intent(out), device :: Ab(7,nx-2,ny-2,nz-2)
    real(4) dudx, dvdy, dwdz, RHS, a0, alpha, alpha_a0
    real(4) over_rhox1, over_rhox2, over_rhoy1, over_rhoy2, over_rhoz1, over_rhoz2
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    alpha = real(sensor(i,j,k)) * (1.e0 / over_dx(i))**2
    dudx = 0.5e0 * real(-ruvwp(2,i-1,j,k) + ruvwp(2,i+1,j,k)) * over_dx(i)
    dvdy = 0.5e0 * real(-ruvwp(3,i,j-1,k) + ruvwp(3,i,j+1,k)) * over_dy(j)
    dwdz = 0.5e0 * real(-ruvwp(4,i,j,k-1) + ruvwp(4,i,j,k+1)) * over_dz(k)
    RHS = (dudx + dvdy + dwdz)**2 + dudx**2 + dvdy**2 + dwdz**2
    over_rhox1 = 0.5e0 * (over_rho(i-1,j,k) + over_rho(i,j,k))
    over_rhox2 = 0.5e0 * (over_rho(i,j,k)   + over_rho(i+1,j,k))
    over_rhoy1 = 0.5e0 * (over_rho(i,j-1,k) + over_rho(i,j,k))
    over_rhoy2 = 0.5e0 * (over_rho(i,j,k)   + over_rho(i,j+1,k))
    over_rhoz1 = 0.5e0 * (over_rho(i,j,k-1) + over_rho(i,j,k))
    over_rhoz2 = 0.5e0 * (over_rho(i,j,k)   + over_rho(i,j,k+1))
    a0 = over_rho(i,j,k) + alpha * ((over_rhox1 + over_rhox2) * over_dx(i)**2 &
                                  + (over_rhoy1 + over_rhoy2) * over_dy(j)**2 &
                                  + (over_rhoz1 + over_rhoz2) * over_dz(k)**2)
    alpha_a0 = alpha / a0
    Ab(1,i-1,j-1,k-1) = over_rhox1 * over_dx(i)**2 * alpha_a0
    Ab(2,i-1,j-1,k-1) = over_rhox2 * over_dx(i)**2 * alpha_a0
    Ab(3,i-1,j-1,k-1) = over_rhoy1 * over_dy(j)**2 * alpha_a0
    Ab(4,i-1,j-1,k-1) = over_rhoy2 * over_dy(j)**2 * alpha_a0
    Ab(5,i-1,j-1,k-1) = over_rhoz1 * over_dz(k)**2 * alpha_a0
    Ab(6,i-1,j-1,k-1) = over_rhoz2 * over_dz(k)**2 * alpha_a0
    Ab(7,i-1,j-1,k-1) = RHS * alpha_a0
  end subroutine set_a


  subroutine set_coef(nx, ny, nz, over_dx, over_dy, over_dz, ruvwp, sensor, Ab)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: over_dx(nx-1), over_dy(ny-1), over_dz(nz-1)
    real(8), intent(in), device  :: ruvwp(5,nx,ny,nz), sensor(nx,ny,nz)
    real(4), intent(out), device :: Ab(7,nx-2,ny-2,nz-2)
    real(4), device :: over_rho(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(3)<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          over_rho(i,j,k) = 1.e0 / real(ruvwp(1,i,j,k))
    enddo;enddo;enddo
    call set_a<<<blocks,threads>>>(nx, ny, nz, over_dx, over_dy, over_dz, over_rho, ruvwp, sensor, Ab)
  end subroutine set_coef


  attributes(device) function update(nx, ny, nz, i, j, k, Ab, sigma) result(sigma_new)
    integer, intent(in), value  :: nx, ny, nz, i, j, k
    real(4), intent(in), device :: Ab(7,nx-2,ny-2,nz-2)
    real(4), intent(in), device :: sigma(nx,ny,nz)
    real(4) sigma_new
    sigma_new = Ab(1,i-1,j-1,k-1) * sigma(i-1,j,k) + Ab(2,i-1,j-1,k-1) * sigma(i+1,j,k) &
              + Ab(3,i-1,j-1,k-1) * sigma(i,j-1,k) + Ab(4,i-1,j-1,k-1) * sigma(i,j+1,k) &
              + Ab(5,i-1,j-1,k-1) * sigma(i,j,k-1) + Ab(6,i-1,j-1,k-1) * sigma(i,j,k+1) + Ab(7,i-1,j-1,k-1)
  end function update


  attributes(global) subroutine point_Jacobi_base(nx, ny, nz, Ab, sigma, sigma_new)
    integer, intent(in), value   :: nx, ny, nz
    real(4), intent(in), device  :: Ab(7,nx-2,ny-2,nz-2)
    real(4), intent(in), device  :: sigma(nx,ny,nz)
    real(4), intent(out), device :: sigma_new(nx,ny,nz)
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1 
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    sigma_new(i,j,k) = update(nx, ny, nz, i, j, k, Ab, sigma)
  end subroutine point_Jacobi_base


  subroutine point_Jacobi(nx, ny, nz, maxitr, Ab, sigma)
    integer, intent(in), value     :: nx, ny, nz, maxitr
    real(4), intent(in), device    :: Ab(7,nx-2,ny-2,nz-2)
    real(4), intent(inout), device :: sigma(nx,ny,nz)
    real(4), device :: sigma_new(nx,ny,nz)
    integer itr
    do itr = 1, maxitr
      call point_Jacobi_base<<<blocks,threads>>>(nx, ny, nz, Ab, sigma, sigma_new)
      call set_bc_sigma(nx, ny, nz, sigma_new)
      sigma(:,:,:) = sigma_new(:,:,:)
    enddo
  end subroutine point_Jacobi

 
  attributes(global) subroutine RB_SOR_base(nx, ny, nz, rb, omega, Ab, sigma)
    integer, intent(in), value     :: nx, ny, nz, rb
    real(4), intent(in), value     :: omega
    real(4), intent(in), device    :: Ab(7,nx-2,ny-2,nz-2)
    real(4), intent(inout), device :: sigma(nx,ny,nz)
    real(4) sigma_new
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1 
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k .or. mod(i+j+k,2) /= rb) return
    sigma_new = update(nx, ny, nz, i, j, k, Ab, sigma)
    sigma(i,j,k) = sigma(i,j,k) + omega * (sigma_new - sigma(i,j,k))
  end subroutine RB_SOR_base


  subroutine RB_SOR(nx, ny, nz, maxitr, Ab, sigma)
    integer, intent(in), value     :: nx, ny, nz, maxitr
    real(4), intent(in), device    :: Ab(7,nx-2,ny-2,nz-2)
    real(4), intent(inout), device :: sigma(nx,ny,nz)
    real(4) :: omega = 1.6e0
    integer itr
    do itr = 1, maxitr
      call RB_SOR_base<<<blocks,threads>>>(nx, ny, nz, 0, omega, Ab, sigma)
      call RB_SOR_base<<<blocks,threads>>>(nx, ny, nz, 1, omega, Ab, sigma)
      call set_bc_sigma(nx, ny, nz, sigma)
    enddo
  end subroutine RB_SOR


  attributes(global) subroutine calc_residual(nx, ny, nz, Ab, sigma, residual)
    integer, intent(in), value   :: nx, ny, nz
    real(4), intent(in), device  :: Ab(7,nx-2,ny-2,nz-2)
    real(4), intent(in), device  :: sigma(nx,ny,nz)
    real(4), intent(out), device :: residual(nx,ny,nz)
    real(4) sigma_new
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1 
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    sigma_new = update(nx, ny, nz, i, j, k, Ab, sigma)
    residual(i,j,k) = sigma_new - sigma(i,j,k)
  end subroutine calc_residual


  subroutine MG_RBSOR(nx, ny, nz, maxitr1, maxitr2, Ab, sigma)
    integer, intent(in), value     :: nx, ny, nz, maxitr1, maxitr2
    real(4), intent(in), device    :: Ab(7,nx-2,ny-2,nz-2)
    real(4), intent(inout), device :: sigma(nx,ny,nz)
    real(4), device :: residual(nx,ny,nz)
    ! pre-smoothing (on fine grid)
    call RB_SOR(nx, ny, nz, maxitr1, Ab, sigma)
    ! compute residual
    call calc_residual<<<blocks,threads>>>(nx, ny, nz, Ab, sigma, residual)
    ! restrict residual
    ! recursive coarse solver
    ! prolongate correction
    ! post-smoothing (on fine grid)
    call RB_SOR(nx, ny, nz, maxitr1, Ab, sigma)
  end subroutine MG_RBSOR


  subroutine calc_sigma(nx, ny, nz, dx, dy, dz, sensor, ruvwp, sigma)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1), dy(ny-1), dz(nz-1), sensor(nx,ny,nz)
    real(8), intent(inout), device :: ruvwp(5,nx,ny,nz)
    real(4), intent(inout), device :: sigma(nx,ny,nz)
    real(4), allocatable, device :: Ab(:,:,:,:)
    integer :: i, j, k, maxitr = 100
    allocate(Ab(7,nx-2,ny-2,nz-2))
    call set_coef(nx, ny, nz, dx, dy, dz, ruvwp, sensor, Ab)
    !call point_Jacobi(nx, ny, nz, maxitr, Ab, sigma)
    call RB_SOR(nx, ny, nz, maxitr, Ab, sigma)
    !call MG_RBSOR(nx, ny, nz, maxitr1, maxitr2, Ab, sigma)
    deallocate(Ab)
    !$cuf kernel do(3)<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ruvwp(5,i,j,k) = ruvwp(5,i,j,k) + dble(sigma(i,j,k))
    enddo;enddo;enddo
  end subroutine calc_sigma
end module calc_igr

