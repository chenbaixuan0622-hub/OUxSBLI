module set_init_common
  use mod_globals, only : gamma, R, Taw, rf
  use mod_constant, only : Cp, gamma_1, over_gamma_1
  use set_compressible_bl
  implicit none
contains
  subroutine calc_Gaussian_filter_x(nx, ny, n, x, phi)
    integer, intent(in)    :: nx, ny, n
    real(8), intent(in)    :: x(nx)
    real(8), intent(inout) :: phi(nx,ny)
    real(8) sigma, d, tmp
    integer :: i, j, ix
    real(8), allocatable :: phi_tmp(:,:), wg(:)
    allocate(phi_tmp(nx,ny))
    phi_tmp = phi
    d = -x(1) + x(2)
    sigma = dble(n) * d / 3.d0
    allocate(wg(2*n+1))
    do ix = -n, n
      wg(ix+n+1) = exp(-((ix*d)**2)/(2.d0*sigma**2))
    enddo
    wg = wg / sum(wg)
    do j = 1, ny
        do i = 1, nx
          tmp = 0.d0
          do ix = -n, n
            if (i+ix >= 1 .and. i+ix <= nx) then
              tmp = tmp + phi_tmp(i+ix,j)*wg(ix+n+1)
            else
              tmp = tmp + phi_tmp(i,j)
            endif
          enddo
          phi(i,j) = tmp
    enddo;enddo
    deallocate(phi_tmp, wg)
  end subroutine calc_Gaussian_filter_x


  subroutine calc_Gaussian_filter_y(nx, ny, n, y, phi)
    integer, intent(in)    :: nx, ny, n
    real(8), intent(in)    :: y(ny)
    real(8), intent(inout) :: phi(nx,ny)
    real(8) sigma, d, tmp
    integer :: i, j, iy
    real(8), allocatable :: phi_tmp(:,:), wg(:)
    allocate(phi_tmp(nx,ny))
    phi_tmp = phi
    d = -y(ny-1) + y(ny)
    sigma = dble(n) * d / 3.d0
    allocate(wg(2*n+1))
    do iy = -n, n
      wg(iy+n+1) = exp(-((iy*d)**2)/(2.d0*sigma**2))
    enddo
    wg = wg / sum(wg)
    do j = 1, ny
        do i = 1, nx
          tmp = 0.d0
          do iy = -n, n
            if (j+iy >= 1 .and. j+iy <= ny) then
              tmp = tmp + phi_tmp(i,j+iy)*wg(iy+n+1)
            else
              tmp = tmp + phi_tmp(i,j)
            endif
          enddo
          phi(i,j) = tmp
    enddo;enddo
    deallocate(phi_tmp, wg)
  end subroutine calc_Gaussian_filter_y


  ! subroutine calc_Gaussian_filter_z(nx, ny, nz, n, z, phi)
  !   integer, intent(in)    :: nx, ny, nz, n
  !   real(8), intent(in)    :: z(nz)
  !   real(8), intent(inout) :: phi(nx,ny,nz)
  !   real(8) sigma, d, tmp
  !   integer :: i, j, k, iz
  !   real(8), allocatable :: phi_tmp(:,:,:), wg(:)
  !   allocate(phi_tmp(nx,ny,nz))
  !   phi_tmp = phi
  !   d = -z(1) + z(2)
  !   sigma = dble(n) * d / 3.d0
  !   allocate(wg(2*n+1))
  !   do iz = -n, n
  !     wg(iz+n+1) = exp(-((iz*d)**2)/(2.d0*sigma**2))
  !   enddo
  !   wg = wg / sum(wg)
  !   do k = 1, nz
  !     do j = 1, ny
  !       do i = 1, nx
  !         tmp = 0.d0
  !         do iz = -n, n
  !           if (k+iz >= 1 .and. k+iz <= nz) then
  !             tmp = tmp + phi_tmp(i,j,k+iz)*wg(iz+n+1)
  !           else
  !             tmp = tmp + phi_tmp(i,j,k)
  !           endif
  !         enddo
  !         phi(i,j,k) = tmp
  !   enddo;enddo;enddo
  !   deallocate(phi_tmp, wg)
  ! end subroutine calc_Gaussian_filter_z


  subroutine set_bc_cyclic_x_cpu(nx, ny, ustd, vstd, Tstd)
    integer, intent(in)    :: nx, ny
    real(8), intent(inout) :: ustd(nx,ny), vstd(nx,ny), Tstd(nx,ny)
    integer j !, k
    !do k = 1, nz
      do j = 1, ny
        ustd(1,j)    = ustd(nx-5,j)
        ustd(2,j)    = ustd(nx-4,j)
        ustd(3,j)    = ustd(nx-3,j)
        ustd(nx-2,j) = ustd(4,j)
        ustd(nx-1,j) = ustd(5,j)
        ustd(nx,j)   = ustd(6,j)
        vstd(1,j)    = vstd(nx-5,j)
        vstd(2,j)    = vstd(nx-4,j)
        vstd(3,j)    = vstd(nx-3,j)
        vstd(nx-2,j) = vstd(4,j)
        vstd(nx-1,j) = vstd(5,j)
        vstd(nx,j)   = vstd(6,j)
        ! wstd(1,j,k)    = wstd(nx-5,j,k)
        ! wstd(2,j,k)    = wstd(nx-4,j,k)
        ! wstd(3,j,k)    = wstd(nx-3,j,k)
        ! wstd(nx-2,j,k) = wstd(4,j,k)
        ! wstd(nx-1,j,k) = wstd(5,j,k)
        ! wstd(nx,j,k)   = wstd(6,j,k)
        Tstd(1,j)    = Tstd(nx-5,j)
        Tstd(2,j)    = Tstd(nx-4,j)
        Tstd(3,j)    = Tstd(nx-3,j)
        Tstd(nx-2,j) = Tstd(4,j)
        Tstd(nx-1,j) = Tstd(5,j)
        Tstd(nx,j)   = Tstd(6,j)
      enddo !enddo
  end subroutine set_bc_cyclic_x_cpu


  ! subroutine set_bc_cyclic_z_cpu(nx, ny, nz, ustd, vstd, wstd, Tstd)
  !   integer, intent(in)    :: nx, ny, nz
  !   real(8), intent(inout) :: ustd(nx,ny,nz), vstd(nx,ny,nz), wstd(nx,ny,nz), Tstd(nx,ny,nz)
  !   integer i, j
  !   do j = 1, ny
  !     do i = 1, nx
  !       ustd(i,j,1)    = ustd(i,j,nz-5)
  !       ustd(i,j,2)    = ustd(i,j,nz-4)
  !       ustd(i,j,3)    = ustd(i,j,nz-3)
  !       ustd(i,j,nz-2) = ustd(i,j,4)
  !       ustd(i,j,nz-1) = ustd(i,j,5)
  !       ustd(i,j,nz)   = ustd(i,j,6)
  !       vstd(i,j,1)    = vstd(i,j,nz-5)
  !       vstd(i,j,2)    = vstd(i,j,nz-4)
  !       vstd(i,j,3)    = vstd(i,j,nz-3)
  !       vstd(i,j,nz-2) = vstd(i,j,4)
  !       vstd(i,j,nz-1) = vstd(i,j,5)
  !       vstd(i,j,nz)   = vstd(i,j,6)
  !       wstd(i,j,1)    = wstd(i,j,nz-5)
  !       wstd(i,j,2)    = wstd(i,j,nz-4)
  !       wstd(i,j,3)    = wstd(i,j,nz-3)
  !       wstd(i,j,nz-2) = wstd(i,j,4)
  !       wstd(i,j,nz-1) = wstd(i,j,5)
  !       wstd(i,j,nz)   = wstd(i,j,6)
  !       Tstd(i,j,1)    = Tstd(i,j,nz-5)
  !       Tstd(i,j,2)    = Tstd(i,j,nz-4)
  !       Tstd(i,j,3)    = Tstd(i,j,nz-3)
  !       Tstd(i,j,nz-2) = Tstd(i,j,4)
  !       Tstd(i,j,nz-1) = Tstd(i,j,5)
  !       Tstd(i,j,nz)   = Tstd(i,j,6)
  !   enddo;enddo
  ! end subroutine set_bc_cyclic_z_cpu


  subroutine calc_rms(nx, ny, phi, rms)
    integer, intent(in)  :: nx, ny
    real(8), intent(in)  :: phi(nx,ny)
    real(8), intent(out) :: rms
    real(8), allocatable :: phi2(:,:)
    real(8) phi_mean, phi2_mean
    allocate(phi2(nx,ny))
    phi2 = phi * phi
    phi2_mean = sum(phi2) / size(phi2)
    phi_mean  = sum(phi)  / size(phi)
    rms = sqrt(phi2_mean - phi_mean * phi_mean)
    deallocate(phi2)
  end subroutine calc_rms


  subroutine set_init_tbl(nx, ny, x, y, rand, blt0, blt, u0, p0, T0, M0, Q)
    integer, intent(in)  :: nx, ny
    real(8), intent(in)  :: x(nx), y(ny)
    real(8), intent(in)  :: rand, blt0, blt, u0, p0, T0, M0
    real(8), intent(out) :: Q(nx,4,ny)
    integer i, j
    integer(4) sz
    integer(4), allocatable :: seed(:)
    real(8) :: p_wall
    real(8) :: fd, pi = acos(-1.d0)
    ! random
    real(8), allocatable :: rho(:), u(:), v(:), T(:), randum(:,:,:), ustd(:,:), vstd(:,:), Tstd(:,:)
    integer ir, jr, nxr, nyr, n
    ! generate randum
    nxr = (nx+9) / 10
    nyr = (ny+1) / 2
    !nzr = (nz+1) / 2
    allocate(rho(ny), u(ny), v(ny), T(ny))
    call calc_HD_Blasius(ny, y, blt0, u0, T0, p0, M0, rho, u, v, T)
    allocate(randum(3,nxr,nyr), ustd(nx,ny), vstd(nx,ny), Tstd(nx,ny))
    call random_seed(size=sz)
    allocate(seed(sz))
    call random_seed(get=seed)
    seed(1) = 48
    call random_seed(put=seed)
    ! generate random number
    !do k = 1, nzr
      do j = 1, nyr
        do i = 1, nxr
          call random_number(randum(1,i,j))
          call random_number(randum(2,i,j))
          call random_number(randum(3,i,j))
          !call random_number(randum(4,i,j))
          randum(1,i,j) = 2.d0 * randum(1,i,j) - 1.d0
          randum(2,i,j) = 2.d0 * randum(2,i,j) - 1.d0
          randum(3,i,j) = 2.d0 * randum(3,i,j) - 1.d0
          !randum(4,i,j) = 2.d0 * randum(4,i,j) - 1.d0
    enddo;enddo
    !enddo
    ! generate fluctuation
    !do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (y(j) <= blt) then
            ir = i / 10 + 1
            jr = j / 2  + 1
            !kr = k / 2  + 1
            ustd(i,j) = 2.d0 * rand * u0 * randum(1,ir,jr)
            vstd(i,j) =        rand * u0 * randum(2,ir,jr)
            !wstd(i,j,k) =        rand * u0 * randum(3,ir,jr,kr)
            !Tstd(i,j,k) = T0 * gamma_1 * M0**2 * randum(4,ir,jr,kr) * rand ! old version. This is not wrong.
            Tstd(i,j) = T0 * gamma_1 * M0**2 * randum(3,ir,jr) * rand * u(j) / u0 ! fixed 2026/03/30
          else
            ustd(i,j) = 0.d0
            vstd(i,j) = 0.d0
            !wstd(i,j,k) = 0.d0
            Tstd(i,j) = 0.d0
          endif
    enddo;enddo
    !enddo
    ! damping
    !do k = 1, nz
      do j = 1, ny
        if (j >= 10) then
          fd = sin(pi * y(j) / (2.d0 * blt))
        else
          fd = 0.d0
        endif
        do i = 1, nx
          ustd(i,j) = fd * ustd(i,j)
          vstd(i,j) = fd * vstd(i,j)
          !wstd(i,j,k) = fd * wstd(i,j,k)
          Tstd(i,j) = fd * Tstd(i,j)
    enddo;enddo;
    !enddo
    call calc_Gaussian_filter_x(nx, ny, 3, x, ustd)
    call calc_Gaussian_filter_y(nx, ny, 3, y, ustd)
    !call calc_Gaussian_filter_z(nx, ny, nz, 3, z, ustd)
    call calc_Gaussian_filter_x(nx, ny, 3, x, vstd)
    call calc_Gaussian_filter_y(nx, ny, 3, y, vstd)
    !call calc_Gaussian_filter_z(nx, ny, nz, 3, z, vstd)
    !call calc_Gaussian_filter_x(nx, ny, 3, x, wstd)
    !call calc_Gaussian_filter_y(nx, ny, 3, y, wstd)
    !call calc_Gaussian_filter_z(nx, ny, nz, 3, z, wstd)
    call calc_Gaussian_filter_x(nx, ny, 3, x, Tstd)
    call calc_Gaussian_filter_y(nx, ny, 3, y, Tstd)
    !call calc_Gaussian_filter_z(nx, ny, nz, 3, z, Tstd)
    call set_bc_cyclic_x_cpu(nx, ny, ustd, vstd, Tstd)
    !call set_bc_cyclic_z_cpu(nx, ny, nz, ustd, vstd, wstd, Tstd)
    ! add fluctuation
    !do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          Q(i,1,j) = p0 / (R * (T(j) + Tstd(i,j)))          !rho(j)
          Q(i,2,j) = Q(i,1,j) * (u(j) + ustd(i,j))
          Q(i,3,j) = Q(i,1,j) * (v(j) + vstd(i,j))
          !Q(i,4,j,k) = Q(i,1,j,k) * wstd(i,j,k)
          Q(i,4,j) = p0 * over_gamma_1 + 0.5d0 * (Q(i,2,j)**2 + Q(i,3,j)**2) / Q(i,1,j)
      enddo;enddo
    !enddo
    deallocate(rho, u, v, T, randum, ustd, vstd, Tstd, seed)
    ! bottom
    Q(:,1,1) = Q(:,1,2)
    Q(:,2,1) = 0.d0
    Q(:,3,1) = 0.d0
    !Q(:,4,1,:) = 0.d0
    p_wall = gamma_1 * (Q(2,4,2) - 0.5d0 * (Q(2,2,2)**2 + Q(2,3,2)**2) / Q(2,1,2))
    Q(:,4,1) = p_wall * over_gamma_1
  end subroutine set_init_tbl
end module set_init_common

