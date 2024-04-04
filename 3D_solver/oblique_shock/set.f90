module set
  use mod_globals, only : nx, ny, nz, dy, gamma
  implicit none
  real(8) :: M0 = 1.9d0
  real(8) :: u0 = 506.8d0
  real(8) :: p0 = 14.924d0
  real(8) :: rho0 = gamma * p0 / M0**2
contains
  subroutine set_init(Q)
    real(8), intent(out), dimension(nx,ny,nz,5) :: Q
    Q(:,:,:,1) = rho0
    Q(:,:,:,2) = rho0 * u0
    Q(:,:,:,3) = 0.d0
    Q(:,:,:,4) = 0.d0
    Q(:,:,:,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
  end subroutine set_init
  
  subroutine set_bc(id_accuracy,Q,T)
    integer(kind=2), intent(in), value :: id_accuracy
    real(8), intent(inout), device :: Q(nx,ny,nz,5)
    real(8), intent(inout), device :: T(nx,ny,nz)
    integer i, j, k, l
    integer :: No = int(0.4 * nx)
    ! slip wall
    real(8) T_bottom, Cp, p_bottom, e_bottom
    ! oblique shock
    real(8) :: beta = 37.2d0 / 180.d0
    real(8) :: Ms = M0 * sin(beta)
    real(8) :: p1 = p0
    real(8) :: p2 = p1 * (1.d0 + 2.d0 * gamma * (Ms**2 - 1.d0) / (gamma + 1.d0)) 
    real(8) :: rho1 = rho0
    real(8) :: rho2 = rho1 * (gamma + 1.d0) * Ms**2 / ((gamma - 1.d0) * Ms**2 + 2.d0)
    real(8) :: u1 = u0 * sin(beta)
    real(8) :: a1 =
    real(8) :: u2 = u1 - 2.d0  * a1 * (Ms**2 - 1.d0 / Ms**2) / (gamma + 1.d0)
    real(8) :: v1 = u0 * cos(beta)
    real(8) :: v2 = u0 * cos(beta) 
    ! inlet
    !$cuf kernel do<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        Q(1,j,k,1) = rho0
        Q(1,j,k,2) = rho0 * u0
        Q(1,j,k,3) = 0.d0
        Q(1,j,k,4) = 0.d0
        Q(1,j,k,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
      enddo
    enddo

    ! outlet
    !$cuf kernel do<<<*,*>>>
    do l = 1, 5
      do k = 1, nz
        do j = 1, ny
          Q(nx,j,k,l) = Q(nx-1,j,k,l)
        enddo
      enddo
    enddo

    ! top
    !$cuf kernel do<<<*,*>>>
    do k = 1, nz
      do i = 1, No
        Q(i,ny,k,1) = rho1
        Q(i,ny,k,2) = rho1 * u0
        Q(i,ny,k,3) = 0.d0
        Q(i,ny,k,4) = 0.d0
        Q(i,ny,k,5) = p1 / (gamma - 1.d0) + 0.5d0 * (Q(i,ny,k,2)**2 + Q(i,ny,k,3)**2) / Q(i,ny,k,1)
      enddo
    enddo
    !$cuf kernel do<<<*,*>>>
    do k = 1, nz
      do i = No + 1, nx
        Q(i,ny,k,1) = rho2
        Q(i,ny,k,2) = rho2 * (v2 * cos(beta) + u2 * sin(beta))
        Q(i,ny,k,3) = rho2 * (-v2 * sin(beta) + u2 * cos(beta))
        Q(i,ny,k,4) = 0.d0
        Q(i,ny,k,5) = p2 / (gamma - 1.d0) + 0.5d0 * (Q(i,ny,k,2)**2 + Q(i,ny,k,3)**2) / Q(i,ny,k,1)
      enddo
    enddo

    ! bottom
    !$cuf kernel do<<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        ! dT/dy = 0
        T_bottom = T(i,2,k)
        T(i,1,k) = T_bottom
        Cp = 1030.5d0 - 0.19975d0 * T_bottom + 3.9734 * T_bottom**2
        ! dp/dy = 0
        p_bottom = (gamma - 1.d0) * (Q(i,2,k,5) - 0.5d0 * (Q(i,2,k,2)**2 + Q(i,2,k,3)**2 + Q(i,2,k,4)**2) / Q(i,2,k,1))
        ! slip wall
        e_bottom = p_bottom / (gamma - 1.d0) + 0.5d0 * (Q(i,2,k,2)**2) / Q(i,2,k,1) 
        Q(i,1,k,1) = (e_bottom + p_bottom) / (Cp * T_bottom + 0.5d0 * (Q(i,2,k,2)/Q(i,2,k,1))**2)
        Q(i,1,k,2) = Q(i,2,k,2)
        Q(i,1,k,3) = 0.d0
        Q(i,1,k,4) = 0.d0
        Q(i,1,k,5) = e_bottom
      enddo
    enddo
  end subroutine set_bc
end module set

