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

  attributes(device) function mu2(T) result(ans)
    real(8), intent(in), device :: T(2)
    real(8) ans
    ans = 0.5d0 * (mu(T(1)) + mu(T(2)))
  end function mu2

  attributes(device) function mu6(T) result(ans)
    real(8), intent(in), device :: T(6)
    real(8) ans(3), mus(6)
    mus(1) = mu(T(1)); mus(2) = mu(T(2)); mus(3) = mu(T(3))
    mus(4) = mu(T(4)); mus(5) = mu(T(5)); mus(6) = mu(T(6))
    ans(:) = 0.0625d0 * (9.d0 * (mus(2:4) + mus(3:5)) - (mus(1:3) + mus(4:6)))
  end function mu6

  attributes(device) function mu23(T) result(ans)
    real(8), intent(in), device :: T(2,3)
    real(8) ans(2)
    ans(1) = 0.25d0 * (mu(T(1,1)) + mu(T(1,2)) + mu(T(2,1)) + mu(T(2,2)))
    ans(2) = 0.25d0 * (mu(T(1,2)) + mu(T(1,3)) + mu(T(2,2)) + mu(T(2,3)))
  end function mu23

  attributes(device) function mu32(T) result(ans)
    real(8), intent(in), device :: T(3,2)
    real(8) ans(2)
    ans(1) = 0.25d0 * (mu(T(1,1)) + mu(T(2,1)) + mu(T(1,2)) + mu(T(2,2)))
    ans(2) = 0.25d0 * (mu(T(2,1)) + mu(T(3,1)) + mu(T(2,2)) + mu(T(3,2)))
  end function mu32

  attributes(device) subroutine calc_kappa(T1,T2,kappa)
    real(8), intent(in), value  :: T1, T2
    real(8), intent(out)        :: kappa
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    kappa =  0.5d0 * (mu(T1) + mu(T2)) * Cp / Pr
  end subroutine calc_kappa
end module calc_sutherland

