module calc_flux
  use mod_globals, only : id_scheme, id_tvd, id_accuracy, accuracy, offset, gamma
  use calc_kep1d
  use calc_slau1d
  implicit none
contains
  attributes(device) function Albada(rho,p,u) result(phi)
    real(8), intent(in), dimension(4), device :: rho, p, u
    real(8) :: e(4), d1, d2, d3, phip, phim, phi, eps = 1.d-16
    e(:) = p(:) / (gamma - 1.d0) + 0.5d0 * rho(:) * u(:)**2 
    d1 = -e(1) / rho(1) + e(2) / rho(2)
    d2 = -e(2) / rho(2) + e(3) / rho(3)
    d3 = -e(3) / rho(3) + e(4) / rho(4)
    phip = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    phim = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    phi = min(phim, phip)
  end function Albada

  attributes(device) function Jameson(p) result(nu)
    real(8), intent(in), dimension(4), device :: p
    real(8) :: nu1, nu2, nu, eps = 1.d-16
    nu1 = abs(p(1) - 2.d0 * p(2) + p(3)) / (abs(p(1)) + 2.d0 * abs(p(2)) + abs(p(3)) + eps)
    nu2 = abs(p(2) - 2.d0 * p(3) + p(4)) / (abs(p(2)) + 2.d0 * abs(p(3)) + abs(p(4)) + eps)
    nu  = max(nu1, nu2)
  end function Jameson

  attributes(device) function Ducros(u) result(fd)
    real(8), intent(in), dimension(4), device :: u
    real(8) :: fd, eps = 1.d-16
    fd = (-u(2) + u(3))**2 / ((-u(2) + u(3))**2 + eps)
  end function Ducros

  attributes(device) function minmod(x,y) result(ans)
    real(8), intent(in), value :: x, y
    real(8) :: ans, sgn
    sgn = sign(1.d0, x)
    ans = sgn * max(min(abs(x), sgn * y), 0.d0)
  end function minmod

  attributes(device) function MUSCL(eps,a) result(alr)
    real(8), intent(in), value  :: eps
    real(8), intent(in), device :: a(4)
    real(8) :: alr(2), d1, d2, d3, dt1, dt2, dt3, dt4
    real(8) :: k = 1.d0 / 3.d0, b
    b = (3.d0 - k) / (1.d0 - k)
    d1 = -a(1) + a(2)
    d2 = -a(2) + a(3)
    d3 = -a(3) + a(4)
    dt1 = minmod(d1, b * d2)
    dt2 = minmod(d2, b * d1)
    dt3 = minmod(d3, b * d2)
    dt4 = minmod(d2, b * d3)
    ! non TVD
    !alr(1) = a(2) + 0.25d0 * eps * ((1.d0 - k) * d1 + (1.d0 + k) * d2)
    !alr(2) = a(3) - 0.25d0 * eps * ((1.d0 - k) * d3 + (1.d0 + k) * d2)
    ! TVD
    alr(1) = a(2) + 0.25d0 * eps * ((1.d0 - k) * dt1 + (1.d0 + k) * dt2)
    alr(2) = a(3) - 0.25d0 * eps * ((1.d0 - k) * dt3 + (1.d0 + k) * dt4)
  end function MUSCL

  attributes(global) subroutine calc_E(nx, rho, u, p, E, sensor)
    integer, intent(in), value                 :: nx
    real(8), intent(in), dimension(nx), device :: rho, u, p
    real(8), intent(out), device               :: E(nx-1,3), sensor(nx-1)
    integer i
    real(8), dimension(2)   :: rho2, p2, u2
    real(8), dimension(4)   :: rho4, p4, u4
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (2 <= i .and. i <= nx-1) then
      rho4(:) = rho(i-1:i+2)
      p4(:)   = p(i-1:i+2)
      u4(:)   = u(i-1:i+2)
      !sensor(i-1) = (1.d0 - Albada(rho4,p4,u4))
      !sensor(i-1) = Jameson(p4)
      !sensor(i-1) = Ducros(u4)
      !sensor(i-1) = Ducros(u4) * (1.d0 - Albada(rho4,p4,u4))
      sensor(i-1) = (1.d0 - Albada(rho4,p4,u4))
      if (id_tvd /= 0) then
        rho2 = MUSCL(sensor(i-1),rho4)
        p2   = MUSCL(sensor(i-1),p4)
        u2   = MUSCL(sensor(i-1),u4)
        if (id_scheme == 0) then
          E(i,:) = KEEP2(rho2,p2,u2)
        elseif (id_scheme == 1) then
          E(i,:) = KEP2(rho2,p2,u2)
        elseif (id_scheme == 2) then
          E(i,:) = SLAU(rho2,p2,u2)
        endif
      else
        if (id_scheme == 0) then
          E(i,:) = KEEP4(rho4,p4,u4)
        elseif (id_scheme == 1) then
          E(i,:) = KEP4(rho4,p4,u4)
        elseif (id_scheme == 2) then
          rho2   = rho4(2:3)
          p2     = p4(2:3)
          u2     = u4(2:3)
          E(i,:) = SLAU(rho2,p2,u2)
        endif
      endif
    else
      rho2(:) = rho(i:i+1)
      p2(:)   = p(i:i+1)
      u2(:)   = u(i:i+1)
      E(i,:) = SLAU(rho2,p2,u2)
    endif
  end subroutine calc_E
end module calc_flux

