module calc_common
  use mod_globals, only : dimension, gamma
  implicit none
contains
  attributes(device) function energy(p,rho,V) result(e)
    real(8), intent(in), value :: p, rho
    real(8), intent(in), device :: V(dimension)
    real(8) :: e
    e = p / (gamma - 1.d0) + 0.5d0 * rho * sum(V**2)
  end function energy

  attributes(device) function ENTHALPY(e,p,rho) result(H)
    real(8), intent(in), value :: e, p, rho
    real(8) :: H
    H = (e + p) / rho
  end function ENTHALPY

  attributes(device) function speed_of_sound(p,rho) result(c)
    real(8), intent(in), value :: p, rho
    real(8) :: c
    c = sqrt(gamma * p / rho)
  end function speed_of_sound
end module calc_common

