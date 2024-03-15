module calc_physical_quantities
  implicit none
contains
  subroutine calc_quantities(nx,gamma,Q,rho,u,p)
    integer, intent(in), value :: nx
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,3) :: Q
    real(8), intent(out), dimension(nx) :: rho, u, p
    integer i
    do i = 1, nx
      rho(i) = Q(i,1)
      u(i) = Q(i,2) / rho(i)
      p(i) = (gamma-1.d0)*(Q(i,3)-0.5d0*rho(i)*(u(i)**2))
    enddo
  end subroutine calc_quantities
  
  subroutine calc_quantities_T(nx,gamma,Q,rho,u,p,T)
    integer, intent(in), value :: nx
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,3) :: Q
    real(8), intent(out), dimension(nx) :: rho, u, p, T
    integer i
    real(8) Cp
    do i = 1, nx
      rho(i) = Q(i,1)
      u(i) = Q(i,2) / rho(i)
      p(i) = (gamma-1.d0)*(Q(i,3)-0.5d0*rho(i)*(u(i)**2))
      Cp = 1030.5d0 - 0.19975d0 * T(i) + 3.9734d0 * T(i) ** 2
      T(i) = ((Q(i,3) + p(i)) / rho(i) - 0.5d0 * (u(i) ** 2)) / Cp
    enddo
  end subroutine calc_quantities_T
end module calc_physical_quantities

