module calc_visc_base
  use mod_globals, only  : Pr
  use mod_constant, only : Cp
  implicit none
  interface tauxx
    module procedure tauxx_visc, tauxx_sgs
  end interface tauxx

  interface tauxy
    module procedure tauxy_visc, tauxy_sgs
  end interface tauxy
contains
  !dir$ inline
  attributes(device) function interpolation6(a) result(ans)
    real(8), intent(in), device :: a(6)
    real(8) ans(3)
    ans(:) = 0.0625d0 * (-a(1:3) + 9.d0 * (a(2:4) + a(3:5)) -a(4:6))
  end function interpolation6

  !dir$ inline
  attributes(device) function dx6(a, dx) result(ans)
    real(8), intent(in), device :: a(6), dx
    real(8) ans(3)
    ans(:) = 0.125d0 * (9.d0 * (-a(2:4) + a(3:5)) - (-a(1:3) + a(4:6)) / 3.d0) * dx
  end function dx6

  !dir$ inline
  attributes(device) function flux4(a) result(ans)
    real(8), intent(in), device :: a(3)
    real(8) ans
    ans = 0.125d0 * ((9.d0 - 1.d0 / 3.d0) * a(2) - (a(1) + a(3)) / 3.d0)
  end function flux4

  attributes(device) subroutine tauxx_visc(mu, ux, vy, wz, u, txx, utxx)
    real(8), intent(in), dimension(3), device :: mu, ux, vy, wz, u
    real(8), intent(out)  :: txx, utxx
    real(8), dimension(3) :: tau, utau
    tau(:)  = 2.d0 * mu(:) * (2.d0 * ux(:) - vy(:) - wz(:)) / 3.d0
    utau(:) = u(:) * tau(:)
    txx     = flux4(tau(:))
    utxx    = flux4(utau(:))
  end subroutine tauxx_visc
  
  attributes(device) subroutine tauxx_sgs(mu, ux, vy, wz, u, txx)
    real(8), intent(in), dimension(3), device :: mu, ux, vy, wz, u
    real(8), intent(out)  :: txx
    real(8), dimension(3) :: tau
    tau(:) = 2.d0 * mu(:) * (2.d0 * ux(:) - vy(:) - wz(:)) / 3.d0
    txx    = flux4(tau(:))
  end subroutine tauxx_sgs

  attributes(device) subroutine tauxy_visc(mu, uy, vx, v, txy, vtxy)
    real(8), intent(in), dimension(3), device :: mu, uy, vx, v
    real(8), intent(out)  :: txy, vtxy
    real(8), dimension(3) :: tau, vtau
    tau(:)  = mu(:) * (uy(:) + vx(:))
    vtau(:) = v(:) * tau(:)
    txy     = flux4(tau(:))
    vtxy    = flux4(vtau(:))
  end subroutine tauxy_visc

  attributes(device) subroutine tauxy_sgs(mu, uy, vx, v, txy)
    real(8), intent(in), dimension(3), device :: mu, uy, vx, v
    real(8), intent(out)  :: txy
    real(8), dimension(3) :: tau
    tau(:) = mu(:) * (uy(:) + vx(:))
    txy    = flux4(tau(:))
  end subroutine tauxy_sgs

  attributes(device) function heat_conduction6(mu, T, dx) result(ans)
    real(8), intent(in), device :: mu(3), T(6)
    real(8), intent(in), value  :: dx
    real(8), dimension(3) ::  kTx
    real(8) :: ans
    kTx(:) = Cp * mu(:) * dx6(T, dx) / Pr
    ans = flux4(kTx(:))
  end function heat_conduction6
end module calc_visc_base

