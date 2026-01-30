module set
  use cudafor
  use mpi
  use mod_globals, only : ny1, gamma, R, Cp, Pr, u0, p0, T0, M0, blt, beta, &
                          ny2, rho2, p2, ux, uy, rf, Taw, rho3, p3, ux3, uy3
  use mod_constant, only : Cp, gamma_1, over_gamma, over_gamma_1
  use set_bc_common
  use set_compressible_bl
  implicit none
  integer No
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(1)
    integer i, j, ny_b, nx1, nx2
    real(8) dx1, dy1, ximp, Lx_s
    dx1  = 20.d0 * blt / dble(512)
    dy1  = dx1
    ximp = 100.d0 * blt

    x(1) = 0.d0
    ! buffer region
    nx1 = 512
    nx2 = nx - 256
    do i = 1, nx1
      dx(i) = dx1 * (1.d0 + 3.d0 * dble(nx1-i) / dble(nx1))
      x(i+1) = x(i) + dx(i)
    enddo
    ! computational region
    do i = nx1, nx2-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo
    ! buffer region
    do i = nx2, nx-1
      dx(i) = dx1 * (1.d0 + 3.d0 * dble(i-nx2) / dble(nx-nx2))
      x(i+1) = x(i) + dx(i)
    enddo
    x(:) = x(:) - ximp

    y(1) = 0.d0
    do j = 1, ny-1
      if (y(j) <= 5.d0 * blt) then
        dy(j) = min(1.d0, max(0.14d0, dble(j)/dble(128))) * dy1
        ny_b  = j
      else
        dy(j) = dy1 * (1.d0 + 1.5d0 * dble(j-ny_b) / dble(ny-ny_b))
      endif
      y(j+1) = y(j) + dy(j)
    enddo

    Lx_s = y(ny) / dble(beta) / blt
    do i = 1, nx
      if (x(i) / blt + Lx_s >= 0.d0) then
        No = i
        exit
      endif
    enddo
    print *, "tan(beta) = ", tan(beta)
    print *, "length of shock = ", Lx_s
    print *, "height of shock = ", y(ny) / blt
    print *, "No    = ", No
    print *, "x(No) = ", x(No) / blt
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(4,nx,ny)
    real(8), allocatable :: rho(:), u(:), v(:), T(:)
    integer i, j
    allocate(rho(ny), u(ny), v(ny), T(ny))
    call calc_HD_Blasius(ny, y, blt, u0, T0, p0, M0, rho, u, v, T)
    do j = 1, ny
      do i = 1, nx
        Q(1,i,j) = rho(j)
        Q(2,i,j) = rho(j) * u(j)
        Q(3,i,j) = rho(j) * v(j)
        Q(4,i,j) = p0 * over_gamma_1 + 0.5d0 * rho(j) * (u(j)**2 + v(j)**2)
    enddo;enddo
    deallocate(rho, u, v, T)
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, x, y, Jacobian, QJ)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: x(nx), y(ny)
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ(4,nx,ny) ! Q / Jacobian
    integer i, j, l
    real(8) :: p_wall
    ! Blasius
    real(8) eta, rho, u, v, T, Tw
    ! Riemann invariants
    real(8) :: rhoin, pin, cin, vin, Rp, Rm, rhob, ub, vb, cb, pb, v0 = 0.d0
    ! cache
    real(8) Jacobian_tmp
    ! temperature and density at top
    real(8), parameter :: Ttop = Taw - rf * u0**2 / (2.d0 * Cp)
    real(8), parameter :: rho0 = p0 / (R * Ttop)
    real(8), parameter :: c0   = sqrt(gamma * p0 / rho0)
    real(8), parameter :: c3   = sqrt(gamma * p3 / rho3)
    !$cuf kernel do(1)<<<*,*>>>
    do j = 2, ny-1
      do l = 1, 4
        ! outlet
        QJ(l,nx,j)   = QJ(l,nx-1,j)
      enddo
    enddo

    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, No-1
      ! Riemann invariants
      Jacobian_tmp = 1.d0 / Jacobian(i,ny)
      pin   = gamma_1 * (QJ(4,i,ny-1) - 0.5d0 * (QJ(2,i,ny-1)**2 + QJ(3,i,ny-1)**2) &
              / QJ(1,i,ny-1)) * Jacobian(i,ny-1)
      rhoin = QJ(1,i,ny-1) * Jacobian(i,ny-1)
      cin   = sqrt(gamma * pin / rhoin)
      vin   = QJ(3,i,ny-1) / QJ(1,i,ny-1)
      Rp    = vin + 2.d0 * cin * over_gamma_1
      Rm    = v0  - 2.d0 * c0  * over_gamma_1
      vb    = 0.5d0 * (Rp + Rm)
      cb    = 0.25d0 * gamma_1 * (Rp - Rm)
      rhob  = (cb / c0)**(2.d0 * over_gamma_1) * rho0
      pb    = (rhob * cb**2) * over_gamma
      QJ(1,i,ny) = rhob * Jacobian_tmp
      QJ(2,i,ny) = rhob * u0 * Jacobian_tmp
      QJ(3,i,ny) = rhob * vb * Jacobian_tmp
      QJ(4,i,ny) = (pb * over_gamma_1 + 0.5d0 * rhob * (u0**2 + vb**2)) * Jacobian_tmp
      ! NoSlip
      QJ(1,i,1) = QJ(1,i,2)
      QJ(2,i,1) = 0.d0
      QJ(3,i,1) = 0.d0
      p_wall = gamma_1 * (QJ(4,i,2) - 0.5d0 * (QJ(2,i,2)**2 + QJ(3,i,2)**2) / QJ(1,i,2))
      QJ(4,i,1) = p_wall * over_gamma_1
    enddo

    !$cuf kernel do(1)<<<*,*>>>
    do i = No, nx
      Jacobian_tmp = 1.d0 / Jacobian(i,ny)
      vin   = QJ(3,i,ny-1) / QJ(1,i,ny-1)
      if (0.5d0 * uy > vin) then 
        QJ(1,i,ny) = rho2 * Jacobian_tmp
        QJ(2,i,ny) = rho2 * ux * Jacobian_tmp
        QJ(3,i,ny) = rho2 * uy * Jacobian_tmp
        QJ(4,i,ny) = (p2 * over_gamma_1 + 0.5d0 * rho2 * (ux**2 + uy**2)) * Jacobian_tmp
      else
        pin   = gamma_1 * (QJ(4,i,ny-1) - 0.5d0 * (QJ(2,i,ny-1)**2 + QJ(3,i,ny-1)**2) &
                / QJ(1,i,ny-1)) * Jacobian(i,ny-1)
        rhoin = QJ(1,i,ny-1) * Jacobian(i,ny-1)
        cin   = sqrt(gamma * pin / rhoin)
        Rp    = vin + 2.d0 * cin * over_gamma_1
        Rm    = uy3 - 2.d0 * c3  * over_gamma_1
        vb    = 0.5d0 * (Rp + Rm)
        cb    = 0.25d0 * gamma_1 * (Rp - Rm)
        rhob  = cin * rhoin / cb
        pb    = (rhob * cb**2) * over_gamma
        ub    = sqrt(2.d0 * gamma * (p3 / rho3 - pb / rhob) * over_gamma_1 + ux3**2 + uy3**2 - vb**2)
        QJ(1,i,ny) = rhob * Jacobian_tmp
        QJ(2,i,ny) = rhob * ub * Jacobian_tmp
        QJ(3,i,ny) = rhob * vb * Jacobian_tmp
        QJ(4,i,ny) = (pb * over_gamma_1 + 0.5d0 * rhob * (ub**2 + vb**2)) * Jacobian_tmp
      endif
      ! NoSlip
      QJ(1,i,1) = QJ(1,i,2)
      QJ(2,i,1) = 0.d0
      QJ(3,i,1) = 0.d0
      p_wall = gamma_1 * (QJ(4,i,2) - 0.5d0 * (QJ(2,i,2)**2 + QJ(3,i,2)**2) / QJ(1,i,2))
      QJ(4,i,1) = p_wall * over_gamma_1
    enddo
  end subroutine set_bc


  subroutine calc_force(nx, ny, x, y, dx, dy, Q, Fout)
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: x(nx), y(ny), dx(nx-1), dy(ny-1)
    real(8), intent(in), device  :: Q(4,nx,ny)
    real(8), intent(out), device :: Fout(3,nx-2,ny-2)
  end subroutine calc_force
end module set

