module calc_hybrid_flux
  use mod_globals, only : accuracy
  implicit none
  real(8), parameter :: eps = 1.d-16
contains
  subroutine calc_E_tvd(nx,Energy,E_keep,E_roe,E_tvd)
    integer, intent(in) :: nx
    real(8), intent(in), dimension(nx) :: Energy
    real(8), intent(in), dimension(nx-accuracy+1,3) :: E_keep, E_roe
    real(8), intent(out), dimension(nx-accuracy+1,3) :: E_tvd
    integer i
    real(8) phi, phi_p, phi_m, d1, d2, d3
    do i = 1, nx-1
      if (i == 1) then
        d1 = 0.d0
      else
        d1 = -Energy(i-1) + Energy(i)
      endif
      d2 = -Energy(i) + Energy(i+1)
      if (i == nx-1) then
        d3 = 0.d0
      else
        d3 = -Energy(i+1) + Energy(i+2)
      endif
      phi_p = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
      phi_m = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
      phi = min(phi_p, phi_m)
      E_tvd(i,:) = phi * E_keep(i,:) + (1.d0 - phi) * E_roe(i,:)
    enddo
  end subroutine calc_E_tvd
end module calc_hybrid_flux

