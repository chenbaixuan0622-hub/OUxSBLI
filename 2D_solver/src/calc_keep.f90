module calc_keep
  use calc_term
  implicit none
  interface calc_E
    module procedure calc_E2, calc_E4
  end interface

  interface calc_F
    module procedure calc_F2, calc_F4
  end interface

contains
  attributes(global) subroutine calc_E2(id, nx, ny, gamma, rho, u, v, p, E)
    integer(kind=2), intent(in), value :: id
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(out), dimension(nx-1,ny-2,4), device :: E
    integer i, j
    real(8) Rho_m, U_m, P_m
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    
    Rho_m = 0.5d0 * (rho(i,j) + rho(i+1,j))
    U_m = 0.5d0 * (u(i,j) + u(i+1,j))
    P_m = 0.5d0 * (p(i,j) + p(i+1,j))
    E(i,j-1,1) = Rho_m * U_m
    E(i,j-1,2) = E(i,j-1,1) * U_m + P_m
    E(i,j-1,3) = E(i,j-1,1) * 0.5d0 * (v(i,j) + v(i+1,j))
    E(i,j-1,4) = E(i,j-1,1) * P_m / ((gamma - 1.d0) * Rho_m) &
    + 0.5d0 * E(i,j-1,1) * (u(i,j) * u(i+1,j) + v(i,j) * v(i+1,j)) &
    + 0.5d0 * (u(i,j) * p(i+1,j) + u(i+1,j) * p(i,j))
  end subroutine calc_E2

  attributes(global) subroutine calc_E4(id, nx, ny, gamma, rho, u, v, p, E)
    integer(kind=4), intent(in), value :: id
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(out), dimension(nx-3,ny-4,4), device :: E
    integer i, j
    real(8), dimension(3) :: RhoU, RhoUU_P, IE, Energy
    real(8), dimension(4) :: rhos, us, vs, ps
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 2
    
    rhos(:) = rho(i:i+3,j)
    us(:) = u(i:i+3,j)
    vs(:) = v(i:i+3,j)
    ps(:) = p(i:i+3,j)
    RhoU(:)  = RhoPhi(rhos(:), us(:))
    RhoUU_P(:) = RhoPhiU(RhoU(:), us(:)) + Phi(ps(:))
    IE(:) = p(i:i+3,j) / ((gamma - 1.d0) * rho(i:i+3,j))
    Energy(:) = IE(:) + RhoUPhiPhi(RhoU(:), us(:), vs(:)) + PhiPsi(us(:), ps(:))
    E(i,j-2,1) = Flux(RhoU(:))
    E(i,j-2,2) = Flux(RhoUU_P(:))
    E(i,j-2,3) = Flux(RhoPhiU(RhoU(:), vs(:)))
    E(i,j-2,4) = Flux(Energy(:))
  end subroutine calc_E4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(global) subroutine calc_F2(id, nx, ny, gamma, rho, u, v, p, F)
    integer(kind=2), intent(in), value :: id
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(out), device :: F(nx-2,ny-1,4)
    integer i, j
    real(8) Rho_m, V_m, P_m
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y

    Rho_m = 0.5d0 * (rho(i,j) + rho(i,j+1))
    V_m = 0.5d0 * (v(i,j) + v(i,j+1))
    P_m = 0.5d0 * (p(i,j) + p(i,j+1))
    F(i-1,j,1) = Rho_m * V_m
    F(i-1,j,2) = F(i-1,j,1) * 0.5d0 * (u(i,j) + u(i,j+1))
    F(i-1,j,3) = F(i-1,j,1) * V_m + P_m
    F(i-1,j,4) = F(i-1,j,1) * P_m / ((gamma - 1.d0) * Rho_m) &
    + 0.5d0 * F(i-1,j,1) * (u(i,j) * u(i,j+1) + v(i,j) * v(i,j+1)) &
    + 0.5d0 * (v(i,j) * p(i,j+1) + v(i,j+1) * p(i,j))
  end subroutine calc_F2

  attributes(global) subroutine calc_F4(id, nx, ny, gamma, rho, u, v, p, F)
    integer(kind=4), intent(in), value :: id
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(out), device :: F(nx-4,ny-3,5)
    integer i, j
    real(8), dimension(3) :: RhoV, RhoVV_P, IE, Energy
    real(8), dimension(4) :: rhos, us, vs, ps
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 2
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y

    rhos(:) = rho(i,j:j+3)
    us(:) = u(i,j:j+3)
    vs(:) = v(i,j:j+3)
    ps(:) = p(i,j:j+3)
    RhoV(:) = RhoPhi(rhos(:), vs(:))
    RhoVV_P(:) = RhoPhiU(RhoV(:), vs(:)) + Phi(ps(:))
    IE(:) = p(i,j:j+3) / ((gamma - 1.0d0) * rho(i,j:j+3))
    Energy(:) = IE(:) + RhoUPhiPhi(RhoV(:), us(:), vs(:)) + PhiPsi(vs(:), ps(:))
    F(i-2,j,1) = Flux(RhoV(:))
    F(i-2,j,2) = Flux(RhoPhiU(RhoV(:), us(:)))
    F(i-2,j,3) = Flux(RhoVV_P(:))
    F(i-2,j,4) = Flux(Energy(:))
  end subroutine calc_F4
end module calc_keep

