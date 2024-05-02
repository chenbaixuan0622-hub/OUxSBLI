module calc_common
  implicit none
contains
  function energy(gamma,p,rho,u) result(e)
    real(8), intent(in) :: gamma, p, rho, u
    real(8) :: e
    e = p / (gamma - 1.d0) + 0.5d0 * rho * u ** 2
  end function energy

  function enthalpy(e,p,rho) result(h)
    real(8), intent(in) :: e, p, rho
    real(8) :: h
    h = (e + p) / rho
  end function enthalpy

  function speed_of_sound(gamma,p,rho) result(c)
    real(8), intent(in) :: gamma, p, rho
    real(8) :: c
    c = sqrt(gamma * p / rho)
  end function speed_of_sound
end module calc_common

