module calc_physical_quantities
  implicit none
contains
  subroutine vecadd(nx,ny,Ev,E)
    integer, intent(in), value :: nx, ny
    real(8), intent(in), device :: Ev(nx,ny,4)
    real(8), intent(inout), device :: E(nx,ny,4)
    integer i, j
    !$acc kernels deviceptr(E,Ev)
    !$acc loop collapse(2)
    do j = 1, ny
      do i = 1, nx
        E(i,j,:) = E(i,j,:) - Ev(i,j,:)
      enddo
    enddo
    !$acc end kernels
  end subroutine

  subroutine calc_quantities(nx,ny,gamma,Q,rho,u,v,p)
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny,4), device :: Q
    real(8), intent(out), dimension(nx,ny), device :: rho, u, v, p
    integer i, j
    !$acc kernels deviceptr(Q,rho,u,v,p)
    !$acc loop collapse(2)
    do j = 1, ny
      do i = 1, nx
        rho(i,j) = Q(i,j,1)
        u(i,j) = Q(i,j,2) / rho(i,j)
        v(i,j) = Q(i,j,3) / rho(i,j)
        p(i,j) = (gamma-1.d0)*(Q(i,j,4)-0.5d0*rho(i,j)*(u(i,j)**2+v(i,j)**2))
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_quantities
  
  subroutine calc_quantities_T(nx,ny,gamma,Q,rho,u,v,p,T)
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny,4), device :: Q
    real(8), intent(out), dimension(nx,ny), device :: rho, u, v, p, T
    integer i, j
    real(8) Cp
    !$acc kernels deviceptr(Q,rho,u,v,p,T)
    !$acc loop collapse(2)
    do j = 1, ny
      do i = 1, nx
        rho(i,j) = Q(i,j,1)
        u(i,j) = Q(i,j,2) / rho(i,j)
        v(i,j) = Q(i,j,3) / rho(i,j)
        p(i,j) = (gamma-1.d0)*(Q(i,j,4)-0.5d0*rho(i,j)*(u(i,j)**2+v(i,j)**2))
        Cp = 1030.5d0 - 0.19975d0 * T(i,j) + 3.9734d0 * T(i,j) ** 2
        T(i,j) = ((Q(i,j,4) + p(i,j)) / rho(i,j) - 0.5d0 * (u(i,j) ** 2 + v(i,j) ** 2)) / Cp
      enddo
    enddo
    !$acc end kernels     
  end subroutine calc_quantities_T
end module calc_physical_quantities

