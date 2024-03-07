module calc_physical_quantities
  implicit none
contains
  subroutine calc_quantities(nx,ny,nz,gamma,Q,rho,u,v,w,p)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q
    real(8), intent(out), dimension(nx,ny,nz), device :: rho, u, v, w, p
    integer i, j, k
    !$acc kernels deviceptr(Q,rho,u,v,w,p)
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
  end subroutine calc_quantities
  
  subroutine calc_quantities_T(nx,ny,nz,gamma,Q,rho,u,v,w,p,T)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q
    real(8), intent(out), dimension(nx,ny,nz), device :: rho, u, v, w, p
    real(8), intent(inout), dimension(nx,ny,nz), device :: T
    integer i, j, k
    real(8) Cp
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
          Cp = 1030.5d0 - 0.19975d0 * T(i,j,k) + 3.9734 * T(i,j,k) ** 2
          T(i,j,k) = ((Q(i,j,k,5) + p(i,j,k)) / rho(i,j,k) &
          & - 0.5d0 * (u(i,j,k) ** 2 + v(i,j,k) ** 2 + w(i,j,k) ** 2)) / Cp
        enddo
      enddo
    enddo
    !$acc end kernels     
  end subroutine calc_quantities_T

  subroutine calc_kinetic_energy(nx,ny,nz,rho,u,v,w,ke)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w
    real(8), intent(out) :: ke
    integer i, j, k
    !$acc kernels deviceptr(rho,u,v,w)
    !$acc loop collapse(3) reduction(+:ke)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ke = rho(i,j,k) * (u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2)
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_kinetic_energy

  subroutine calc_entropy(nx,ny,nz,Cv,Cp,rho,p,s)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: Cv, Cp
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, p
    real(8), intent(out) :: s
    integer i, j, k
    !$acc kernels deviceptr(rho,p)
    !$acc loop collapse(3) reduction(+:s)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          s = Cv * log(p(i,j,k)) - Cp * log(rho(i,j,k))
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_entropy
end module calc_physical_quantities
  
