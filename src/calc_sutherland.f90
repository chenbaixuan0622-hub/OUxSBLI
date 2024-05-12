module calc_sutherland
  use mod_globals, only : gamma, R, Pr
  implicit none
  interface calc_mu
    module procedure calc_mu2, calc_mu4
  end interface
contains
  attributes(device) function mu(T) result(ans)
    real(8), intent(in), value :: T
    real(8) :: ans
    real(8) :: mu0 = 1.716d-5
    real(8) :: T0 = 273.2d0
    real(8) :: S = 111.d0
    ans = mu0 * ((T0 + S) / (T + S)) * (T / T0) ** 1.5d0
  end function mu

  attributes(device) subroutine calc_mu2(T1,T2,mu_mean)
    real(8), intent(in), value  :: T1, T2
    real(8), intent(out)        :: mu_mean
    mu_mean = 0.5d0 * (mu(T1) + mu(T2)) 
  end subroutine calc_mu2

  attributes(device) subroutine calc_mu4(T1,T2,T3,T4,mu_mean)
    real(8), intent(in), value  :: T1, T2, T3, T4
    real(8), intent(out)        :: mu_mean
    mu_mean = 0.25d0 * (mu(T1) + mu(T2) + mu(T3) + mu(T4)) 
  end subroutine calc_mu4

  attributes(device) subroutine calc_kappa(T1,T2,kappa)
    real(8), intent(in), value  :: T1, T2
    real(8), intent(out)        :: kappa
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    kappa =  0.5d0 * (mu(T1) + mu(T2)) * Cp / Pr
  end subroutine calc_kappa
end module calc_sutherland

