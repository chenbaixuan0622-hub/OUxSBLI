module set_compressible_bl
  use mod_globals, only : R, rf,  Taw
  use mod_constant, only : Cp
  implicit none
  integer, parameter :: Neta = 29
  real(8), parameter :: eta_tab(Neta) = (/ &
  0.d0, 0.2d0, 0.4d0, 0.6d0, 0.8d0, 1.d0, 1.2d0, 1.4d0, 1.6d0, 1.8d0, &
  2.d0, 2.2d0, 2.4d0, 2.6d0, 2.8d0, 3.d0, 3.2d0, 3.4d0, 3.6d0, 3.8d0, &
  4.d0, 4.2d0, 4.4d0, 4.6d0, 4.8d0, 5.d0, 6.0d0, 7.0d0, 8.0d0/)
  real(8), parameter :: fp_tab(Neta) = (/ &
  0.0000d0, 0.06643d0, 0.1327d0, 0.19855d0, 0.2638d0, &
  0.3282d0, 0.39160d0, 0.4538d0, 0.51450d0, 0.5736d0, &
  0.6310d0, 0.68660d0, 0.7404d0, 0.79230d0, 0.8423d0, &
  0.8904d0, 0.93650d0, 0.9606d0, 0.97620d0, 0.9866d0, &
  0.9931d0, 0.99670d0, 0.9987d0, 0.99960d0, 0.9999d0, &
  1.0000d0, 1.00000d0, 1.0000d0, 1.00000d0/)
  real(8), parameter :: f_tab(Neta) = (/ &
  0.00000d0, 0.002214d0, 0.008849d0, 0.019901d0, 0.03536d0, &
  0.05522d0, 0.079470d0, 0.108100d0, 0.141000d0, 0.17820d0, &
  0.21960d0, 0.265200d0, 0.314800d0, 0.368500d0, 0.42610d0, &
  0.48750d0, 0.552700d0, 0.621600d0, 0.694200d0, 0.77030d0, &
  0.84980d0, 0.932600d0, 1.018600d0, 1.107600d0, 1.19950d0, &
  1.29420d0, 1.788100d0, 2.347000d0, 2.964100d0/)
contains
  function lin_interp(n, x, x_tab, y_tab) result(y)
    integer, intent(in) :: n
    real(8), intent(in) :: x
    real(8), intent(in) :: x_tab(n), y_tab(n)
    integer :: i
    real(8) :: y, w
    if (x <= x_tab(1)) then
      lin_interp = y_tab(1)
      return
    endif
    do i = 1, n-1
      if (x_tab(i) <= x .and. x <= x_tab(i+1)) then
        w = (x - x_tab(i)) / (x_tab(i+1) - x_tab(i))
        lin_interp = (1.0d0 - w) * y_tab(i) + w * y_tab(i+1)
        return
      endif
    enddo
    lin_interp = y_tab(n)
  end function lin_interp


  subroutine calc_HD_Blasius(ny, y, blt0, u0, T0, p0, M0, rho_out, u_out, v_out, T_out)
    integer, intent(in) :: ny
    real(8), intent(in) :: y(ny), blt0, u0, T0, p0, M0
    real(8), intent(out) :: rho_out(ny), u_out(ny), v_out(ny), T_out(ny)
    real(8), allocatable :: rho(:), u(:), v(:), T(:), T_r(:), integral_part(:), y_physical(:)
    real(8) :: eta99 = 4.91d0, Tw = Taw
    real(8) coeff, I_99, v_coeff
    integer j
    allocate(rho(Neta), u(Neta), v(Neta), T(Neta), T_r(Neta), integral_part(Neta), y_physical(Neta))
    do j = 1, Neta
      u(j)   = u0 * fp_tab(j) ! Blasius
      T(j)   = Tw + (Taw - Tw) * u(j) / u0 - rf * u(j)**2 / (2.0d0 * Cp) ! Crocco-Busemann
      T_r(j) = T(j) / T0
      rho(j) = p0 / (R * T(j))
    enddo
    integral_part(1) = 0.0d0 ! Howarth-Dorodnitsyn transform
    do j = 2, Neta
      integral_part(j) = integral_part(j-1) + 0.5d0 * (1.d0 / T_r(j-1) + 1.d0 / T_r(j)) * (-eta_tab(j-1) + eta_tab(j))
    enddo
    I_99 = 0.0d0
    do j = 2, Neta
      if (eta_tab(j) <= eta99) then
        I_99 = I_99 + 0.5d0 * (T_r(j-1) + T_r(j)) * (eta_tab(j) - eta_tab(j-1))
      else
        exit
      endif
    enddo
    coeff   = I_99 / blt0
    v_coeff = 0.5d0 * u0 / coeff
    do j = 1, Neta
      y_physical(j) = integral_part(j) / coeff
      v(j) = v_coeff * (integral_part(j) * fp_tab(j) - T_r(j) * f_tab(j))
    enddo
    do j = 1, ny
      rho_out(j) = lin_interp(Neta, y(j), y_physical, rho)
      u_out(j)   = lin_interp(Neta, y(j), y_physical, u)
      v_out(j)   = lin_interp(Neta, y(j), y_physical, v)
      T_out(j)   = lin_interp(Neta, y(j), y_physical, T)
    enddo
    deallocate(rho, u, v, T, T_r, integral_part, y_physical)
  end subroutine calc_HD_Blasius
end module set_compressible_bl

