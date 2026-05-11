module set_bc_tbl_sbli
  use mod_globals, only : R, gamma
  use mod_constant, only : gamma_1, over_gamma, over_gamma_1, Cp
  implicit none
contains
  subroutine set_bc_Neumann_tbl_top_down(nx, ny, offset, istart, iend, Jacobian, QJ)
    use mod_globals, only : Taw, rf, u0, p0
    integer, intent(in), value     :: nx, ny, offset, istart, iend
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ(nx,4,ny)
    real(8) :: Jacobian_tmp, over_QJ1, p_wall, rhob, ub, vb, pb
    integer i   !, k
    !$cuf kernel do(1)<<<*,*>>>
    !do k = 1+offset, nz-offset
      do i = istart, iend
        Jacobian_tmp = 1.d0 / Jacobian(i,ny)
        over_QJ1 = 1.d0 / QJ(i,1,ny-1)         !J/rho?
        rhob = QJ(i,1,ny-1) * Jacobian(i,ny-1)  !(rho/j)*j?
        ub   = QJ(i,2,ny-1) * over_QJ1          !(rho*u/j)*(j/rho)?
        vb   = QJ(i,3,ny-1) * over_QJ1          !(rho*v/j)*(j/rho)?
        !wb   = QJ(i,4,ny-1,k) * over_QJ1
        pb   = gamma_1 * (QJ(i,4,ny-1) * Jacobian(i,ny-1) - 0.5d0 * rhob * (ub**2 + vb**2))
        QJ(i,1,ny) = rhob * Jacobian_tmp
        QJ(i,2,ny) = rhob * ub * Jacobian_tmp
        QJ(i,3,ny) = rhob * vb * Jacobian_tmp
        !QJ(i,4,ny,k) = rhob * wb * Jacobian_tmp
        QJ(i,4,ny) = (pb * over_gamma_1 + 0.5d0 * rhob * (ub**2 + vb**2)) * Jacobian_tmp
        ! NoSlip
        QJ(i,1,1) = QJ(i,1,2)  !rho/j
        QJ(i,2,1) = 0.d0      !(rho*u)/j
        QJ(i,3,1) = 0.d0      !(rho*v)/j
        !QJ(i,4,1,k) = 0.d0
        p_wall = gamma_1 * (QJ(i,4,2) - 0.5d0 * (QJ(i,2,2)**2 + QJ(i,3,2)**2) / QJ(i,1,2))
        QJ(i,4,1) = p_wall * over_gamma_1
      enddo
    !enddo
  end subroutine set_bc_Neumann_tbl_top_down


  subroutine set_bc_Riemann_tbl_top_down(nx, ny, offset, istart, iend, Jacobian, QJ)
    use mod_globals, only : Taw, rf, u0, p0
    integer, intent(in), value     :: nx, ny, offset, istart, iend
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ(nx,4,ny)
    real(8), parameter :: T       = Taw - rf * u0**2 / (2.d0 * Cp)
    real(8), parameter :: rho0    = p0 / (R * T)
    real(8), parameter :: c0      = sqrt(gamma * p0 / rho0)
    real(8), parameter :: over_c0 = 1.d0 / c0
    real(8) :: Jacobian_tmp, p_wall, rhoin, pin, cin, vin, Rp, Rm, rhob, vb, cb, pb, v0 = 0.d0
    integer i   !, k
    !$cuf kernel do(1)<<<*,*>>>
    !do k = 1+offset, nz-offset
      do i = istart, iend
        Jacobian_tmp = 1.d0 / Jacobian(i,ny)
        pin   = gamma_1 * (QJ(i,4,ny-1) - 0.5d0 * (QJ(i,2,ny-1)**2 + QJ(i,3,ny-1)**2) &
                / QJ(i,1,ny-1)) * Jacobian(i,ny-1)
        rhoin = QJ(i,1,ny-1) * Jacobian(i,ny-1)
        cin   = sqrt(gamma * pin / rhoin)
        vin   = QJ(i,3,ny-1) / QJ(i,1,ny-1)
        Rp   = vin + 2.d0 * cin * over_gamma_1
        Rm   = v0  - 2.d0 * c0  * over_gamma_1
        vb   = 0.5d0 * (Rp + Rm)
        cb   = 0.25d0 * gamma_1 * (Rp - Rm)
        rhob = (cb * over_c0)**(2.d0 * over_gamma_1) * rho0
        pb   = (rhob * cb**2) * over_gamma
        QJ(i,1,ny) = rhob * Jacobian_tmp
        QJ(i,2,ny) = rhob * u0 * Jacobian_tmp
        QJ(i,3,ny) = rhob * vb * Jacobian_tmp
        !QJ(i,4,ny,k) = 0.d0
        QJ(i,4,ny) = (pb * over_gamma_1 + 0.5d0 * rhob * (u0**2 + vb**2)) * Jacobian_tmp
        ! NoSlip
        QJ(i,1,1) = QJ(i,1,2)
        QJ(i,2,1) = 0.d0
        QJ(i,3,1) = 0.d0
        !QJ(i,4,1,k) = 0.d0
        p_wall = gamma_1 * (QJ(i,4,2) - 0.5d0 * (QJ(i,2,2)**2 + QJ(i,3,2)**2) / QJ(i,1,2))
        QJ(i,4,1) = p_wall * over_gamma_1
      enddo
    !enddo
  end subroutine set_bc_Riemann_tbl_top_down
end module set_bc_tbl_sbli
