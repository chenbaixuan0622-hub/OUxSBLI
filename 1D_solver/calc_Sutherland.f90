module calc_Sutherland
  implicit none
contains
  attributes(device) function mu(T) result(ans)
    real(8), intent(in), value :: T
    real(8) :: ans
    ans = (1.4592d-6 * T ** 1.5d0) / (109.1d0 + T)
  end function mu

  attributes(device) subroutine calc_mu(T1,T2,mu_mean)
    real(8), intent(in), value :: T1, T2
    real(8), intent(out) :: mu_mean
    real(8) mu1, mu2
    mu_mean = 0.5d0 * (mu(T1) + mu(T2)) 
  end subroutine calc_mu

  attributes(device) subroutine calc_kappa(T1,T2,kappa)
    real(8), intent(in), value :: T1, T2
    real(8), intent(out) :: kappa
    real(8) kappa1, kappa2
    kappa1 = (2.334d-3 * T1 ** 1.5d0) / (164.54d0 + T1)
    kappa2 = (2.334d-3 * T2 ** 1.5d0) / (164.54d0 + T2)
    kappa =  0.5d0 * (kappa1 + kappa2)
  end subroutine calc_kappa
end module calc_Sutherland

