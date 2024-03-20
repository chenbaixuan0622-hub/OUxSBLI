module calc_common
  use mod_globals, only : gamma
  implicit none
  interface energy
    module procedure energy_2D, energy_3D
  end interface

contains
  attributes(device) function energy_2D(p,rho,u,v) result(e)
    real(8), intent(in), value :: p, rho, u , v
    real(8) :: e
    e = p / (gamma - 1.d0) + 0.5d0 * rho * (u ** 2 + v ** 2)
  end function energy_2D
  
  attributes(device) function energy_3D(p,rho,u,v,w) result(e)
    real(8), intent(in), value :: p, rho, u , v, w
    real(8) :: e
    e = p / (gamma - 1.d0) + 0.5d0 * rho * (u ** 2 + v ** 2 + w ** 2)
  end function energy_3D

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

