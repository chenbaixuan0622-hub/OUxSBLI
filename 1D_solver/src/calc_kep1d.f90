module calc_kep1d
  use mod_globals, only : gamma
  implicit none
contains
  attributes(device) function KEP2(rho,p,u) result(F)
    real(8), intent(in), dimension(2), device :: rho, p, u
    real(8), dimension(3) :: F
    integer i
    real(8) KE, H
    F(1) = 0.25d0 * (rho(1) + rho(2)) * (u(1) + u(2))

    F(2) = 0.5d0 * (F(1) * (u(1) + u(2)) + (p(1) + p(2)))

    ! original kinetic energy term
    KE = 0.25d0 * F(1) * (u(1)**2 + u(2)**2)

    ! original enthalpy term
    H = F(1) * 0.5d0 * gamma / (gamma - 1.d0) * (p(1) / rho(1) + p(2) / rho(2))
    !H = F(1) * 0.5d0 * gamma / (gamma - 1.d0) * (p(1) / rho(1) + p(2) / rho(2) - abs(F(1)) * (p(2) / rho(2) - p(1) / rho(1)))
    ! KEEP PE like enthalpy term
    !H = 0.25d0 * (u(1) + u(2)) * gamma * (p(1) + p(2)) / (gamma - 1.d0)

    ! KEP energy term
    F(3) = KE + H
  end function KEP2

  attributes(device) function KEEP2(rho,p,u) result(F)
    real(8), intent(in), dimension(2), device :: rho, p, u
    real(8), dimension(3) :: F
    integer i
    real(8) KE, IE, PV
    F(1) = 0.25d0 * (rho(1) + rho(2)) * (u(1) + u(2))

    F(2) = 0.5d0 * (F(1) * (u(1) + u(2)) + (p(1) + p(2)))

    ! KEEP kinetic energy term
    KE = 0.5d0 * F(1) * u(1) * u(2)
    
    ! KEEP internal energy term
    IE = F(1) * 0.5d0 * (p(1) / rho(1) + p(2) / rho(2)) / (gamma - 1.d0)
    ! KEEPPE internal energy term
    !IE = 0.25d0 * (u(1) + u(2)) * (p(1) + p(2)) / (gamma - 1.d0)

    ! KEEP pressure diffusion term
    PV = 0.5d0 * (u(1) * p(2) + u(2) * p(1))

    ! KEEP energy term
    F(3) = KE + IE + PV
  end function KEEP2

  attributes(device) function KEP4(rho,p,u) result(F)
    real(8), intent(in), dimension(4), device :: rho, p, u
    real(8), dimension(3) :: F
    real(8) KE, H
    F(1) = 0.25d0 * (rho(2) + rho(3)) * (u(2) + u(3))

    F(2) = 0.5d0 * (F(1) * (u(2) + u(3)) + (p(2) + p(3)))

    KE = 0.25d0 * F(1) * (u(2)**2 + u(3)**2)
    ! original enthalpy term
    H = F(1) * 0.5d0 * gamma * (p(2) / rho(2) + p(3) / rho(3)) / (gamma - 1.d0)
    ! KEEP PE like enthalpy term
    !H = 0.25d0 * (u(2) + u(3)) * gamma * (p(2) + p(3)) / (gamma - 1.d0)

    F(3) = KE + H
  end function KEP4

  attributes(device) function KEEP4(rho,p,u) result(F)
    real(8), intent(in), dimension(4), device :: rho, p, u
    real(8), dimension(3) :: F
    real(8) KE, IE, PV
    F(1) = 0.25d0 * (rho(2) + rho(3)) * (u(2) + u(3))

    F(2) = 0.5d0 * (F(1) * (u(2) + u(3)) + (p(2) + p(3)))

    KE = 0.5d0 * F(1) * u(2) * u(3)
    ! KEEP internal energy term
    IE = F(1) * 0.5d0 * (p(2) / rho(2) + p(3) / rho(3)) / (gamma - 1.d0)
    ! KEEPPE internal energy term
    !IE = 0.25d0 * (u(2) + u(3)) * (p(2) + p(3)) / (gamma - 1.d0)
    PV = 0.5d0 * (u(2) * p(3) + u(3) * p(2))

    F(3) = KE + IE + PV
  end function KEEP4
end module calc_kep1d

