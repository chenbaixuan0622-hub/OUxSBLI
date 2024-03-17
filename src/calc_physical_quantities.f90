module calc_physical_quantities
  use mod_globals, only : nx, ny, nz, gamma
  implicit none
  interface 
    subroutine calc_quantities_2D(Q,rho,u,v,p,T)  
      real(8), intent(in), dimension(nx,ny,4), device :: Q
      real(8), intent(out), dimension(nx,ny), device :: rho, u, v, p
      real(8), intent(out), dimension(nx,ny), optional, device :: T
      integer i, j
      real(8) Cp
    end subroutine calc_quantities_2D
  end interface

  interface
    subroutine calc_quantities_3D(Q,rho,u,v,w,p,T)
      real(8), intent(in), dimension(nx,ny,nz,5), device :: Q
      real(8), intent(out), dimension(nx,ny,nz), device :: rho, u, v, w, p
      real(8), intent(out), dimension(nx,ny,nz), optional, device :: T
      integer i, j, k
      real(8) Cp
    end subroutine calc_quantities_3D
  end interface

contains
  subroutine calc_quantities_2D(Q,rho,u,v,p,T)
    real(8), intent(in), dimension(nx,ny,4), device :: Q
    real(8), intent(out), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(out), dimension(nx,ny), optional, device :: T
    integer i, j
    real(8) Cp
    !$acc kernels deviceptr(Q,rho,u,v,p,T)
    !$acc loop collapse(2) private(Cp)
    do j = 1, ny
      do i = 1, nx
        rho(i,j) = Q(i,j,1)
        u(i,j) = Q(i,j,2) / rho(i,j)
        v(i,j) = Q(i,j,3) / rho(i,j)
        p(i,j) = (gamma - 1.d0) * (Q(i,j,4)- 0.5d0 * rho(i,j) * (u(i,j)**2 + v(i,j)**2))
        Cp = 1030.5d0 - 0.19975d0 * T(i,j) + 3.9734d0 * T(i,j) ** 2
        T(i,j) = ((Q(i,j,4) + p(i,j)) / rho(i,j) - 0.5d0 * (u(i,j)**2 + v(i,j)**2)) / Cp
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_quantities_2D

  subroutine calc_quantities_3D(Q,rho,u,v,w,p,T)
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q
    real(8), intent(out), dimension(nx,ny,nz), device :: rho, u, v, w, p
    real(8), intent(out), dimension(nx,ny,nz), optional, device :: T
    integer i, j, k
    real(8) Cp
    !$acc kernels deviceptr(Q,rho,u,v,w,p,T)
    !$acc loop collapse(3) private(Cp)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rho(i,j,k) = Q(i,j,k,1)
          u(i,j,k) = Q(i,j,k,2) / rho(i,j,k)
          v(i,j,k) = Q(i,j,k,3) / rho(i,j,k)
          w(i,j,k) = Q(i,j,k,4) / rho(i,j,k)
          p(i,j,k) = (gamma - 1.d0) * (Q(i,j,k,5) - 0.5d0 * rho(i,j,k) * (u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2))
          Cp = 1030.5d0 - 0.19975d0 * T(i,j,k) + 3.9734d0 * T(i,j,k) ** 2
          T(i,j,k) = ((Q(i,j,k,5) + p(i,j,k)) / rho(i,j,k) - 0.5d0 * (u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2)) / Cp
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_quantities_3D
end module calc_physical_quantities

