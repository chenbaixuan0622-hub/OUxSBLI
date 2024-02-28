module calc_physical_quantities
  implicit none
  interface calc_quantities
    module procedure calc_quantities_noT, calc_quantities_T
  end interface  
contains
  subroutine calc_quantities_noT(nx,ny,nz,gamma,Q,rho,u,v,w,p)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q
    real(8), intent(out), dimension(nx,ny,nz), device :: rho, u, v, w, p
    integer i, j, k
    !$acc kernels deviceptr(Q,rho,u,v,w,p,T)
    !$acc loop collapse(3)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rho(i,j,k) = Q(i,j,k,1)
          u(i,j,k) = Q(i,j,k,2) / rho(i,j,k)
          v(i,j,k) = Q(i,j,k,3) / rho(i,j,k)
          w(i,j,k) = Q(i,j,k,4) / rho(i,j,k)
          p(i,j,k) = (gamma-1.d0)*(Q(i,j,k,5)-0.5d0*rho(i,j,k)*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2))
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_quantities_noT
  
  subroutine calc_quantities_T(nx,ny,nz,gamma,Cp,Q,rho,u,v,w,p,T)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: gamma, Cp
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q
    real(8), intent(out), dimension(nx,ny,nz), device :: rho, u, v, w, p, T
    integer i, j, k
    !$acc kernels deviceptr(Q,rho,u,v,w,p,T)
    !$acc loop collapse(3)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rho(i,j,k) = Q(i,j,k,1)
          u(i,j,k) = Q(i,j,k,2) / rho(i,j,k)
          v(i,j,k) = Q(i,j,k,3) / rho(i,j,k)
          w(i,j,k) = Q(i,j,k,4) / rho(i,j,k)
          p(i,j,k) = (gamma-1.d0)*(Q(i,j,k,5)-0.5d0*rho(i,j,k)*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2))
          T(i,j,k) = ((Q(i,j,k,5) + p(i,j,k)) / rho(i,j,k) &
          & - 0.5d0 * (u(i,j,k) ** 2 + v(i,j,k) ** 2 + w(i,j,k) ** 2)) / Cp
        enddo
      enddo
    enddo
    !$acc end kernels     
  end subroutine calc_quantities_T
end module calc_physical_quantities
  