module set
  use cudafor
  use mpi
  use mod_globals, only : gamma, R, Cp, Pr, u0, rho0, p0, T0, M0, blt, beta, &
                          rho2, p2, ux, uy, rf, Taw, rho3, p3, ux3, uy3
  use mod_constant, only : Cp, gamma_1, over_gamma, over_gamma_1
  use set_bc_common
  use set_init_common
  implicit none
  integer No
contains
  subroutine set_grid(myrank, nx, ny, Lx, Ly, x, y, dx, dy)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: Lx, Ly
    real(8), intent(out) :: x(nx), y(ny), dx(nx-1), dy(ny-1)
    integer i, j, ny_b
    real(8) dx1, dy1, ximp, Lx_s
    dx1  = Lx / dble(nx-1)
    dy1  = dx1
    ximp = 70.d0 * blt

    x(1) = 0.d0
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo
    x(:) = x(:) - ximp

    y(1) = 0.d0
    do j = 1, ny-1
      if (y(j) <= 3.d0 * blt) then
        dy(j) = min(1.d0, max(0.07d0, dble(j)/dble(128))) * dy1
        ny_b  = j
      else
        dy(j) = dy1 * (1.d0 + 0.75d0 * dble(j-ny_b) / dble(ny-ny_b))
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
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, xs, ys, Q)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: xs(nx), ys(ny)
    real(8), intent(out) :: Q(nx,4,ny)
    call set_init_tbl(nx, ny, xs, ys, blt, blt, u0, p0, T0, M0, Q)
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, Jacobian, QJ)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ(nx,4,ny) ! Q / Jacobian
    integer i, j, k
    real(8) :: p_wall, pre, rho, rhou, rhov, rhow, p, e
    ! Riemann invariants
    real(8) :: rhoin, pin, cin, vin, Rp, Rm, rhob, ub, vb, cb, pb, v0 = 0.d0
    ! cache
    real(8) Jacobian_tmp
    ! temperature and density at top
    real(8), parameter :: T    = Taw - rf * u0**2 / (2.d0 * Cp)
    real(8), parameter :: rho0 = p0 / (R * T)
    real(8), parameter :: c0   = sqrt(gamma * p0 / rho0)
    real(8), parameter :: c3   = sqrt(gamma * p3 / rho3)
    !$cuf kernel do(1)<<<*,*>>>
    do j = 2, ny-1
      do k = 1, 4
        ! cyclic
        QJ(1,k,j)    = QJ(nx-5,k,j)
        QJ(2,k,j)    = QJ(nx-4,k,j)
        QJ(3,k,j)    = QJ(nx-3,k,j)
        QJ(nx-2,k,j) = QJ(4,k,j)
        QJ(nx-1,k,j) = QJ(5,k,j)
        QJ(nx,k,j)   = QJ(6,k,j)
        ! outlet
        !QJ(nx,k,j) = QJ(nx-1,k,j)
    enddo;enddo
 
    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      ! top
      Jacobian_tmp = 1.d0 / Jacobian(i,ny)
      QJ(i,1,ny) = rho0 * Jacobian_tmp
      QJ(i,2,ny) = rho0 * u0 * Jacobian_tmp
      QJ(i,3,ny) = 0.d0
      QJ(i,4,ny) = (p0 * over_gamma_1 + 0.5d0 * rho0 * u0**2) * Jacobian_tmp
      ! NoSlip
      QJ(i,1,1) = QJ(i,1,2)
      QJ(i,2,1) = 0.d0
      QJ(i,3,1) = 0.d0
      p_wall = gamma_1 * (QJ(i,4,2) - 0.5d0 * (QJ(i,2,2)**2 + QJ(i,3,2)**2) / QJ(i,1,2))
      QJ(i,4,1) = p_wall * over_gamma_1
    enddo

    !!$cuf kernel do(1)<<<*,*>>>
    !do i = No, nx
    !  Jacobian_tmp = 1.d0 / Jacobian(i,ny)
    !  vin = QJ(i,3,ny-1) / QJ(i,1,ny-1)
    !  if (0.5d0 * uy > vin) then 
    !    QJ(i,1,ny) = rho2 * Jacobian_tmp
    !    QJ(i,2,ny) = rho2 * ux * Jacobian_tmp
    !    QJ(i,3,ny) = rho2 * uy * Jacobian_tmp
    !    QJ(i,4,ny) = (p2 * over_gamma_1 + 0.5d0 * rho2 * (ux**2 + uy**2)) * Jacobian_tmp
    !  else
    !    pin   = gamma_1 * (QJ(i,4,ny-1) - 0.5d0 * (QJ(i,2,ny-1)**2 + QJ(i,3,ny-1)**2) &
    !            / QJ(i,1,ny-1)) * Jacobian(i,ny-1)
    !    rhoin = QJ(i,1,ny-1) * Jacobian(i,ny-1)
    !    cin   = sqrt(gamma * pin / rhoin)
    !    Rp    = vin + 2.d0 * cin * over_gamma_1
    !    Rm    = uy3 - 2.d0 * c3  * over_gamma_1
    !    vb    = 0.5d0 * (Rp + Rm)
    !    cb    = 0.25d0 * gamma_1 * (Rp - Rm)
    !    rhob  = cin * rhoin / cb
    !    pb    = (rhob * cb**2) * over_gamma
    !    ub    = sqrt(2.d0 * gamma * (p3 / rho3 - pb / rhob) * over_gamma_1 + ux3**2 + uy3**2 - vb**2)
    !    QJ(i,1,ny) = rhob * Jacobian_tmp
    !    QJ(i,2,ny) = rhob * ub * Jacobian_tmp
    !    QJ(i,3,ny) = rhob * vb * Jacobian_tmp
    !    QJ(i,4,ny) = (pb * over_gamma_1 + 0.5d0 * rhob * (ub**2 + vb**2)) * Jacobian_tmp
    !  endif
    !enddo
  end subroutine set_bc
end module set

