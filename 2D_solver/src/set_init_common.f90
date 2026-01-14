module set_init_common
  use cufft
  use curand
  use mod_globals, only : pi
  implicit none
contains
  subroutine calc_init_pressure(accuracy, nx_a, ny_a, nx_h, plan_fwd, plan_inv, rho0, p0, norm, u_h, v_h, p)
    integer, intent(in)            :: accuracy, nx_a, ny_a, nx_h, plan_fwd, plan_inv
    real(8), intent(in)            :: rho0, p0, norm
    complex(8), intent(in), device :: u_h(nx_h,ny_a), v_h(nx_h,ny_a)
    real(8), intent(out)           :: p(nx_a,ny_a)
    complex(8), allocatable, device :: temp_h(:,:), s_h(:,:)
    real(8), allocatable, device    :: ux(:,:), uy(:,:), vx(:,:), vy(:,:), source(:,:), p_d(:,:)
    integer :: i, j, istat
    real(8) :: kx, ky, k2
    complex(8) :: eye
    eye = (0.d0, 1.d0)
    allocate(temp_h(nx_h,ny_a),ux(nx_a,ny_a), uy(nx_a,ny_a), vx(nx_a,ny_a), vy(nx_a,ny_a))
    ! ux = ifft(ikx * u_h)
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny_a
      do i = 1, nx_h
        kx = dble(i-1)
        temp_h(i,j) = eye * kx * u_h(i,j)
    enddo;enddo
    istat = cufftExecZ2D(plan_inv, temp_h, ux)
    ! uy = ifft(iky * u_h)
    !$cuf kernel do(1) <<<*,*>>>
    do j = 1, ny_a
      ky = merge(dble(j-1), dble(j-1-ny_a), j <= ny_a/2+1)
      temp_h(:,j) = eye * ky * u_h(:,j)
    enddo
    istat = cufftExecZ2D(plan_inv, temp_h, uy)
    ! vx = ifft(ikx * v_h)
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny_a
      do i = 1, nx_h
        kx = dble(i-1)
        temp_h(i,j) = eye * kx * v_h(i,j)
    enddo;enddo
    istat = cufftExecZ2D(plan_inv, temp_h, vx)
    ! vy = ifft(iky * v_h)
    !$cuf kernel do(1) <<<*,*>>>
    do j = 1, ny_a
      ky = merge(dble(j-1), dble(j-1-ny_a), j <= ny_a/2+1)
      temp_h(:,j) = eye * ky * v_h(:,j)
    enddo
    istat = cufftExecZ2D(plan_inv, temp_h, vy)
    deallocate(temp_h)
    allocate(source(nx_a,ny_a))
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny_a
      do i = 1, nx_a
        ! source term: S = - rho (ux^2 + 2*uy*vx + vy^2)
        source(i,j) = -rho0 * ( (ux(i,j)*norm)**2 + 2.d0*(uy(i,j)*norm)*(vx(i,j)*norm) + (vy(i,j)*norm)**2 )
    enddo;enddo
    deallocate(uy, vx, vy)
    allocate(s_h(nx_h,ny_a))
    istat = cufftExecD2Z(plan_fwd, source, s_h)
    deallocate(source)
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny_a
      do i = 1, nx_h
        ky = merge(dble(j-1), dble(j-1-ny_a), j <= ny_a/2+1)
        kx = dble(i-1)
        k2 = kx*kx + ky*ky
        if (k2 > 0.d0) then
          s_h(i,j) = s_h(i,j) / (-k2)
        else
          s_h(i,j) = (0.d0, 0.d0)
        endif
    enddo;enddo
    istat = cufftExecZ2D(plan_inv, s_h, ux) ! ux is a temporary matrix
    deallocate(s_h)
    allocate(p_d(nx_a,ny_a))
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny_a
      do i = 1, nx_a
        p_d(i,j) = (ux(i,j) * norm)! + p0
    enddo;enddo
    p = p_d
    deallocate(ux, p_d)
  end subroutine calc_init_pressure


  subroutine calc_spectra(accuracy, nx, ny, k0, energy_scale, rho0, p0, u, v, p)
    integer, intent(in)  :: accuracy, nx, ny, k0
    real(8), intent(in)  :: energy_scale, rho0, p0
    real(8), intent(out) :: u(nx-accuracy,ny-accuracy), v(nx-accuracy,ny-accuracy)
    real(8), intent(out) :: p(nx-accuracy,ny-accuracy)
    complex(8), allocatable, device :: u_h(:,:), v_h(:,:)
    real(8), allocatable, device    :: u_d(:,:), v_d(:,:)
    real(8) kx, ky, k2, kabs, E_k, amp, dot_prod, norm
    real(4), allocatable, device :: phi(:)
    type(curandGenerator) gen
    integer i, j, istat, nx_a, ny_a, nx_h, plan_fwd, plan_inv
    nx_a = nx - accuracy
    ny_a = ny - accuracy
    nx_h = nx_a / 2 + 1
    norm = 1.d0 / dble((nx_a) * (ny_a))
    allocate(u_h(nx_h,ny_a), v_h(nx_h,ny_a))
    allocate(u_d(nx_a,ny_a), v_d(nx_a,ny_a))
    allocate(phi(nx_h*ny_a*2))
    istat = cufftPlan2d(plan_fwd, nx_a, ny_a, CUFFT_D2Z)
    istat = cufftPlan2d(plan_inv, nx_a, ny_a, CUFFT_Z2D)
    istat = curandCreateGenerator(gen, CURAND_RNG_PSEUDO_DEFAULT)
    istat = curandSetPseudoRandomGeneratorSeed(gen, 12345_8)
    istat = curandGenerateUniform(gen, phi, 2*nx_h*ny_a)
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny_a
      do i = 1, nx_h
        kx = dble(i-1)
        if (j <= ny_a/2+1) then
          ky = dble(j-1)
        else
          ky = dble(j-1-ny_a)
        endif
        k2 = kx*kx + ky*ky
        kabs = sqrt(k2)
        if (k2 > 0.d0) then
          E_k = energy_scale * (kabs**4) / ((kabs + dble(k0))**7)
          amp = sqrt(E_k / (pi * kabs))
          u_h(i,j) = amp * exp((0.d0, 1.d0) * 2.d0 * pi * phi(i+(j-1)*nx_h))
          v_h(i,j) = amp * exp((0.d0, 1.d0) * 2.d0 * pi * phi(i+(j-1)*nx_h+nx_h*ny_a))
          dot_prod = (kx * u_h(i,j) + ky * v_h(i,j)) / k2
          u_h(i,j) = u_h(i,j) - dot_prod * kx
          v_h(i,j) = v_h(i,j) - dot_prod * ky
        else
          u_h(i,j) = (0.d0, 0.d0)
          v_h(i,j) = (0.d0, 0.d0)
        endif
    enddo;enddo
    istat = cufftExecZ2D(plan_inv, u_h, u_d)
    istat = cufftExecZ2D(plan_inv, v_h, v_d)
    call calc_init_pressure(accuracy, nx_a, ny_a, nx_h, plan_fwd, plan_inv, rho0, p0, norm, u_h, v_h, p)
    istat = cufftDestroy(plan_fwd)
    istat = cufftDestroy(plan_inv)
    istat = curandDestroyGenerator(gen)
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny_a
      do i = 1, nx_a
        u_d(i,j) = u_d(i,j) * norm
        v_d(i,j) = v_d(i,j) * norm
    enddo;enddo
    u = u_d
    v = v_d
    deallocate(u_h, v_h, u_d, v_d, phi)
  end subroutine calc_spectra
end module set_init_common

