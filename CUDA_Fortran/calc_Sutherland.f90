module calc_Sutherland
  implicit none
  interface calc_mu
    module procedure calc_mu2, calc_mu4
  end interface
contains
  attributes(device) subroutine calc_mu2(T1,T2,mu)
    real(8), intent(in), value :: T1, T2
    real(8), intent(out) :: mu
    real(8) mu1, mu2
    mu1 = (1.4592d-6 * T1 ** (1.5d0)) / (109.1d0 + T1) 
    mu2 = (1.4592d-6 * T2 ** (1.5d0)) / (109.1d0 + T2)
    mu = 0.5d0 * (mu1 + mu2) 
  end subroutine calc_mu2

  attributes(device) subroutine calc_mu4(T1,T2,T3,T4,mu)
    real(8), intent(in), value :: T1, T2, T3, T4
    real(8), intent(out) :: mu
    real(8) mu1, mu2, mu3, mu4
    mu1 = (1.4592d-6 * T1 ** (1.5d0)) / (109.1d0 + T1) 
    mu2 = (1.4592d-6 * T2 ** (1.5d0)) / (109.1d0 + T2) 
    mu3 = (1.4592d-6 * T3 ** (1.5d0)) / (109.1d0 + T3) 
    mu4 = (1.4592d-6 * T4 ** (1.5d0)) / (109.1d0 + T4)
    mu = 0.25d0 * (T1 + T2 + T3 + T4) 
  end subroutine calc_mu4

  attributes(device) subroutine calc_kappa(T1,T2,kappa)
    real(8), intent(in), value :: T1, T2
    real(8), intent(out) :: kappa
    real(8) kappa1, kappa2
    kappa1 = (2.334d-3 * T1 ** 1.5d0) / (164.54d0 + T1)
    kappa2 = (2.334d-3 * T2 ** 1.5d0) / (164.54d0 + T2)
    kappa =  0.5d0 * (kappa1 + kappa2)
  end subroutine calc_kappa
end module calc_Sutherland
