module calc_physical_quantities
  implicit none
  interface calc_quantities
    module procedure calc_quantities_noT, calc_quantities_T
  end interface  
contains
  subroutine calc_quantities_noT(nx,ny,gamma,Q,rho,u,v,p)
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
  end subroutine calc_quantities_noT
  
  subroutine calc_quantities_T(nx,ny,gamma,Cp,Q,rho,u,v,p,T)
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: gamma, Cp
    real(8), intent(in), dimension(nx,ny,4), device :: Q
    real(8), intent(out), dimension(nx,ny), device :: rho, u, v, p, T
    integer i, j
    !$acc kernels deviceptr(Q,rho,u,v,p,T)
    !$acc loop collapse(2)
    do j = 1, ny
      do i = 1, nx
        rho(i,j) = Q(i,j,1)
        u(i,j) = Q(i,j,2) / rho(i,j)
        v(i,j) = Q(i,j,3) / rho(i,j)
        p(i,j) = (gamma-1.d0)*(Q(i,j,4)-0.5d0*rho(i,j)*(u(i,j)**2+v(i,j)**2))
        T(i,j) = ((Q(i,j,4) + p(i,j)) / rho(i,j) - 0.5d0 * (u(i,j) ** 2 + v(i,j) ** 2)) / Cp
      enddo
    enddo
    !$acc end kernels     
  end subroutine calc_quantities_T
end module calc_physical_quantities
  
