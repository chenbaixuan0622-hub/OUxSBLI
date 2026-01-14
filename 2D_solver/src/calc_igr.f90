module calc_igr
  use mod_globals, only : blocks, threads
  implicit none
contains
  subroutine set_init_sigma(nx, ny, dx, dy, Jacobian, Q, sigma_gpu)
    integer, intent(in)          :: nx, ny
    real(8), intent(in)          :: dx(nx-1), dy(ny-1), Jacobian(nx,ny), Q(4,nx,ny)
    real(4), intent(out), device :: sigma_gpu(nx,ny)
    real(4), allocatable :: Ab(:,:,:), sigma(:,:), over_rho(:,:), u(:,:), v(:,:)
    real(4) :: alpha, dudx, dvdy, RHS, a0, alpha_a0, omega = 1.6e0
    real(4) over_rhox1, over_rhox2, over_rhoy1, over_rhoy2, sigma_new
    integer :: i, j, itr, maxitr = 100
    ! set coef
    allocate(Ab(5,nx-2,ny-2), sigma(nx,ny), over_rho(nx,ny), u(nx,ny), v(nx,ny))
    do j = 1, ny
      do i = 1, nx
        over_rho(i,j) = 1.e0 / real(Q(1,i,j) * Jacobian(i,j))
        u(i,j) = real(Q(2,i,j) / Q(1,i,j))
        v(i,j) = real(Q(3,i,j) / Q(1,i,j))
    enddo;enddo
    do j = 2, ny-1
      do i = 2, nx-1
        alpha = 0.5e0 * dx(i)**2
        dudx = 0.5e0 * (-u(i-1,j) + u(i+1,j)) / dx(i)
        dvdy = 0.5e0 * (-v(i,j-1) + v(i,j+1)) / dy(j)
        RHS = (dudx + dvdy)**2 + dudx**2 + dvdy**2
        over_rhox1 = 0.5e0 * (over_rho(i-1,j) + over_rho(i,j))
        over_rhox2 = 0.5e0 * (over_rho(i,j)   + over_rho(i+1,j))
        over_rhoy1 = 0.5e0 * (over_rho(i,j-1) + over_rho(i,j))
        over_rhoy2 = 0.5e0 * (over_rho(i,j)   + over_rho(i,j+1))
        a0 = over_rho(i,j) + alpha * ((over_rhox1 + over_rhox2) / dx(i)**2 &
                                    + (over_rhoy1 + over_rhoy2) / dy(j)**2)
        alpha_a0 = alpha / a0
        Ab(1,i-1,j-1) = over_rhox1 / dx(i)**2 * alpha_a0
        Ab(2,i-1,j-1) = over_rhox2 / dx(i)**2 * alpha_a0
        Ab(3,i-1,j-1) = over_rhoy1 / dy(j)**2 * alpha_a0
        Ab(4,i-1,j-1) = over_rhoy2 / dy(j)**2 * alpha_a0
        Ab(5,i-1,j-1) = RHS * alpha_a0
    enddo;enddo
    deallocate(over_rho, u, v)
    do itr = 1, maxitr
      do j = 2, ny-1
        do i = 2, nx-1
          sigma_new = Ab(1,i-1,j-1) * sigma(i-1,j) + Ab(2,i-1,j-1) * sigma(i+1,j) &
                    + Ab(3,i-1,j-1) * sigma(i,j-1) + Ab(4,i-1,j-1) * sigma(i,j+1) + Ab(5,i-1,j-1)
          sigma(i,j) = sigma(i,j) + omega * (sigma_new - sigma(i,j))
      enddo;enddo
      ! set bc sigma
      sigma(1,:)  = sigma(2,:)
      sigma(nx,:) = sigma(nx-1,:)
      sigma(:,1)  = sigma(:,2)
      sigma(:,ny) = sigma(:,ny-1)
    enddo
    sigma_gpu(:,:) = sigma(:,:)
    deallocate(Ab, sigma)
  end subroutine set_init_sigma


  subroutine set_bc_sigma(nx, ny, sigma)
    integer, intent(in), value     :: nx, ny
    real(4), intent(inout), device :: sigma(nx,ny)
    integer i, j
    !$cuf kernel do(1)<<<*,*>>>
    do j = 1, ny
      sigma(1,j)  = sigma(2,j)
      sigma(nx,j) = sigma(nx-1,j)
    enddo
    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      sigma(i,1)  = sigma(i,2)
      sigma(i,ny) = sigma(i,ny-1)
    enddo
  end subroutine set_bc_sigma


  attributes(global) subroutine set_a(nx, ny, over_dx, over_dy, over_rho, ruvp, sensor, Ab)
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: over_dx(nx-1), over_dy(ny-1)
    real(4), intent(in), device  :: over_rho(nx,ny)
    real(8), intent(in), device  :: ruvp(4,nx,ny), sensor(nx,ny)
    real(4), intent(out), device :: Ab(5,nx-2,ny-2)
    real(4) dudx, dvdy, RHS, a0, alpha, alpha_a0
    real(4) over_rhox1, over_rhox2, over_rhoy1, over_rhoy2
    integer i, j
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    if (nx-1 < i .or. ny-1 < j) return
    !alpha = real(sensor(i,j)) * (1.e0 / over_dx(i))**2
    alpha = (1.e0 / over_dx(i))**2
    dudx = 0.5e0 * real(-ruvp(2,i-1,j) + ruvp(2,i+1,j)) * over_dx(i)
    dvdy = 0.5e0 * real(-ruvp(3,i,j-1) + ruvp(3,i,j+1)) * over_dy(j)
    RHS = (dudx + dvdy)**2 + dudx**2 + dvdy**2
    over_rhox1 = 0.5e0 * (over_rho(i-1,j) + over_rho(i,j))
    over_rhox2 = 0.5e0 * (over_rho(i,j)   + over_rho(i+1,j))
    over_rhoy1 = 0.5e0 * (over_rho(i,j-1) + over_rho(i,j))
    over_rhoy2 = 0.5e0 * (over_rho(i,j)   + over_rho(i,j+1))
    a0 = over_rho(i,j) + alpha * ((over_rhox1 + over_rhox2) * over_dx(i)**2 &
                                + (over_rhoy1 + over_rhoy2) * over_dy(j)**2)
    alpha_a0 = alpha / a0
    Ab(1,i-1,j-1) = over_rhox1 * over_dx(i)**2 * alpha_a0
    Ab(2,i-1,j-1) = over_rhox2 * over_dx(i)**2 * alpha_a0
    Ab(3,i-1,j-1) = over_rhoy1 * over_dy(j)**2 * alpha_a0
    Ab(4,i-1,j-1) = over_rhoy2 * over_dy(j)**2 * alpha_a0
    Ab(5,i-1,j-1) = RHS * alpha_a0
  end subroutine set_a


  subroutine set_coef(nx, ny, over_dx, over_dy, ruvp, sensor, Ab)
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: over_dx(nx-1), over_dy(ny-1)
    real(8), intent(in), device  :: ruvp(4,nx,ny), sensor(nx,ny)
    real(4), intent(out), device :: Ab(5,nx-2,ny-2)
    real(4), device :: over_rho(nx,ny)
    integer i, j
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        over_rho(i,j) = 1.e0 / real(ruvp(1,i,j))
    enddo;enddo
    call set_a<<<blocks,threads>>>(nx, ny, over_dx, over_dy, over_rho, ruvp, sensor, Ab)
  end subroutine set_coef


  attributes(device) function update(nx, ny, i, j, Ab, sigma) result(sigma_new)
    integer, intent(in), value  :: nx, ny, i, j
    real(4), intent(in), device :: Ab(5,nx-2,ny-2)
    real(4), intent(in), device :: sigma(nx,ny)
    real(4) sigma_new
    sigma_new = Ab(1,i-1,j-1) * sigma(i-1,j) + Ab(2,i-1,j-1) * sigma(i+1,j) &
              + Ab(3,i-1,j-1) * sigma(i,j-1) + Ab(4,i-1,j-1) * sigma(i,j+1) + Ab(5,i-1,j-1)
  end function update

 
  attributes(global) subroutine RB_SOR_base(nx, ny, rb, omega, Ab, sigma)
    integer, intent(in), value     :: nx, ny, rb
    real(4), intent(in), value     :: omega
    real(4), intent(in), device    :: Ab(5,nx-2,ny-2)
    real(4), intent(inout), device :: sigma(nx,ny)
    real(4) sigma_new
    integer i, j
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1 
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    if (nx-1 < i .or. ny-1 < j .or. mod(i+j,2) /= rb) return
    sigma_new = update(nx, ny, i, j, Ab, sigma)
    sigma(i,j) = sigma(i,j) + omega * (sigma_new - sigma(i,j))
  end subroutine RB_SOR_base


  subroutine RB_SOR(nx, ny, maxitr, Ab, sigma)
    integer, intent(in), value     :: nx, ny, maxitr
    real(4), intent(in), device    :: Ab(5,nx-2,ny-2)
    real(4), intent(inout), device :: sigma(nx,ny)
    real(4) :: omega = 1.6e0
    integer itr
    do itr = 1, maxitr
      call RB_SOR_base<<<blocks,threads>>>(nx, ny, 0, omega, Ab, sigma)
      call RB_SOR_base<<<blocks,threads>>>(nx, ny, 1, omega, Ab, sigma)
      call set_bc_sigma(nx, ny, sigma)
    enddo
  end subroutine RB_SOR


  attributes(global) subroutine calc_residual(nx, ny, Ab, sigma, residual)
    integer, intent(in), value   :: nx, ny
    real(4), intent(in), device  :: Ab(5,nx-2,ny-2)
    real(4), intent(in), device  :: sigma(nx,ny)
    real(4), intent(out), device :: residual(nx,ny)
    real(4) sigma_new
    integer i, j
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1 
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    if (nx-1 < i .or. ny-1 < j) return
    sigma_new = update(nx, ny, i, j, Ab, sigma)
    residual(i,j) = sigma_new - sigma(i,j)
  end subroutine calc_residual


  subroutine calc_sigma(nx, ny, dx, dy, sensor, ruvp, sigma)
    integer, intent(in), value     :: nx, ny
    real(8), intent(in), device    :: dx(nx-1), dy(ny-1), sensor(nx,ny)
    real(8), intent(in), device    :: ruvp(4,nx,ny)
    real(4), intent(inout), device :: sigma(nx,ny)
    real(4), allocatable, device :: Ab(:,:,:)
    integer :: maxitr = 10
    allocate(Ab(5,nx-2,ny-2))
    call set_coef(nx, ny, dx, dy, ruvp, sensor, Ab)
    call RB_SOR(nx, ny, maxitr, Ab, sigma)
    deallocate(Ab)
  end subroutine calc_sigma
end module calc_igr

