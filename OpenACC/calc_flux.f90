module calc_flux
  use calc_term
  implicit none
  interface calc_E
    module procedure calc_E2, calc_E4
  end interface

  interface calc_F
    module procedure calc_F2, calc_F4
  end interface

  interface calc_G
    module procedure calc_G2, calc_G4
  end interface

contains
  subroutine calc_E2(id, nx, ny, nz, gamma, rho, u, v, w, p, E)
    integer(kind=2), intent(in), value :: id
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w, p
    real(8), intent(out), dimension(nx-1,ny-2,nz-2,5), device :: E
    integer i, j, k
    real(8) Rho_m, U_m, P_m
    !$acc kernels deviceptr(rho,u,v,w,p,E)
    !$acc loop collapse(3) private(Rho_m, U_m, P_m)
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 1, nx-1
          Rho_m = 0.5d0 * (rho(i,j,k) + rho(i+1,j,k))
          U_m = 0.5d0 * (u(i,j,k) + u(i+1,j,k))
          P_m = 0.5d0 * (p(i,j,k) + p(i+1,j,k))
          E(i,j-1,k-1,1) = Rho_m * U_m
          E(i,j-1,k-1,2) = E(i,j-1,k-1,1) * U_m + P_m
          E(i,j-1,k-1,3) = E(i,j-1,k-1,1) * 0.5d0 * (v(i,j,k) + v(i+1,j,k))
          E(i,j-1,k-1,4) = E(i,j-1,k-1,1) * 0.5d0 * (w(i,j,k) + w(i+1,j,k))
          E(i,j-1,k-1,5) = E(i,j-1,k-1,1) * P_m / ((gamma - 1.d0) * Rho_m) &
          + 0.5d0 * E(i,j-1,k-1,1) * (u(i,j,k) * u(i+1,j,k) + v(i,j,k) * v(i+1,j,k) + w(i,j,k) * w(i+1,j,k)) &
          + 0.5d0 * (u(i,j,k) * p(i+1,j,k) + u(i+1,j,k) * p(i,j,k))
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_E2

  subroutine calc_E4(id, nx, ny, nz, gamma, rho, u, v, w, p, E)
    integer(kind=4), intent(in), value :: id
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8), intent(out), dimension(nx-3,ny-4,nz-4,5) :: E
    integer i, j, k
    real(8), dimension(3) :: RhoU, RhoUU_P, IE, Energy
    do k = 3, nz-2
      do j = 3, ny-2
        do i = 1, nx-3
          RhoU(:)  = RhoPhi(rho(i:i+3,j,k), u(i:i+3,j,k))
          RhoUU_P(:) = RhoPhiU(RhoU(:), u(i:i+3,j,k)) + Phi(p(i:i+3,j,k))
          IE(:) = p(i:i+3,j,k) / ((gamma - 1.d0) * rho(i:i+3,j,k))
          Energy(:) = IE(:) + RhoUPhiPhi(RhoU(:), u(i:i+3,j,k), v(i:i+3,j,k), w(i:i+3,j,k)) + PhiPsi(u(i:i+3,j,k), p(i:i+3,j,k))
          E(i,j-2,k-2,1) = Flux(RhoU(:))
          E(i,j-2,k-2,2) = Flux(RhoUU_P(:))
          E(i,j-2,k-2,3) = Flux(RhoPhiU(RhoU(:), v(i:i+3,j,k)))
          E(i,j-2,k-2,4) = Flux(RhoPhiU(RhoU(:), w(i:i+3,j,k)))
          E(i,j-2,k-2,5) = Flux(Energy(:))
        enddo
      enddo
    enddo
  end subroutine calc_E4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_F2(id, nx, ny, nz, gamma, rho, u, v, w, p, F)
    integer(kind=2), intent(in), value :: id
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w, p
    real(8), intent(out), device :: F(nx-2,ny-1,nz-2,5)
    integer i, j, k
    real(8) Rho_m, V_m, P_m
    !$acc kernels deviceptr(rho,u,v,w,p,F)
    !$acc loop collapse(3) private(Rho_m, V_m, P_m)
    do k = 2, nz-1
      do j = 1, ny-1
        do i = 2, nx-1
          Rho_m = 0.5d0 * (rho(i,j,k) + rho(i,j+1,k))
          V_m = 0.5d0 * (v(i,j,k) + v(i,j+1,k))
          P_m = 0.5d0 * (p(i,j,k) + p(i,j+1,k))
          F(i-1,j,k-1,1) = Rho_m * V_m
          F(i-1,j,k-1,2) = F(i-1,j,k-1,1) * 0.5d0 * (u(i,j,k) + u(i,j+1,k))
          F(i-1,j,k-1,3) = F(i-1,j,k-1,1) * V_m + P_m
          F(i-1,j,k-1,4) = F(i-1,j,k-1,1) * 0.5d0 * (w(i,j,k) + w(i,j+1,k))
          F(i-1,j,k-1,5) = F(i-1,j,k-1,1) * P_m / ((gamma - 1.d0) * Rho_m) &
          + 0.5d0 * F(i-1,j,k-1,1) * (u(i,j,k) * u(i,j+1,k) + v(i,j,k) * v(i,j+1,k) + w(i,j,k) * w(i,j+1,k)) &
          + 0.5d0 * (v(i,j,k) * p(i,j+1,k) + v(i,j+1,k) * p(i,j,k))
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_F2

  subroutine calc_F4(id, nx, ny, nz, gamma, rho, u, v, w, p, F)
    integer(kind=4), intent(in), value :: id
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8), intent(out) :: F(nx-4,ny-3,nz-4,5)
    integer i, j, k
    real(8), dimension(3) :: RhoV, RhoVV_P, IE, Energy
    do k = 3, nz-2
      do j = 1, ny-3
        do i = 3, nx-2
          RhoV(:) = RhoPhi(rho(i,j:j+3,k), v(i,j:j+3,k))
          RhoVV_P(:) = RhoPhiU(RhoV(:), v(i,j:j+3,k)) + Phi(p(i,j:j+3,k))
          IE(:) = p(i,j:j+3,k) / ((gamma - 1.0d0) * rho(i,j:j+3,k))
          Energy(:) = IE(:) + RhoUPhiPhi(RhoV(:), u(i,j:j+3,k), v(i,j:j+3,k), w(i,j:j+3,k)) + PhiPsi(v(i,j:j+3,k), p(i,j:j+3,k))
          F(i-2,j,k-2,1) = Flux(RhoV(:))
          F(i-2,j,k-2,2) = Flux(RhoPhiU(RhoV(:), u(i,j:j+3,k)))
          F(i-2,j,k-2,3) = Flux(RhoVV_P(:))
          F(i-2,j,k-2,4) = Flux(RhoPhiU(RhoV(:), w(i,j:j+3,k)))
          F(i-2,j,k-2,5) = Flux(Energy(:))
        enddo
      enddo
    enddo
  end subroutine calc_F4


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_G2(id, nx, ny, nz, gamma, rho, u, v, w, p, G)
    integer(kind=2), intent(in), value :: id
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w, p
    real(8), intent(out), device :: G(nx-2,ny-2,nz-1,5)
    integer i, j, k
    real(8) Rho_m, W_m, P_m
    !$acc kernels deviceptr(rho,u,v,w,p,G)
    !$acc loop collapse(3) private(Rho_m, W_m, P_m)
    do k = 1, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          Rho_m = 0.5d0 * (rho(i,j,k) + rho(i,j,k+1))
          W_m = 0.5d0 * (w(i,j,k) + w(i,j,k+1))
          P_m = 0.5d0 * (p(i,j,k) + p(i,j,k+1))
          G(i-1,j-1,k,1) = Rho_m * W_m
          G(i-1,j-1,k,2) = G(i-1,j-1,k,1) * 0.5d0 * (u(i,j,k) + u(i,j,k+1))
          G(i-1,j-1,k,3) = G(i-1,j-1,k,1) * 0.5d0 * (v(i,j,k) + v(i,j,k+1))
          G(i-1,j-1,k,4) = G(i-1,j-1,k,1) * W_m + P_m
          G(i-1,j-1,k,5) = G(i-1,j-1,k,1) * P_m / ((gamma - 1.d0) * Rho_m) &
          + 0.5d0 * G(i-1,j-1,k,1) * (u(i,j,k) * u(i,j,k+1) + v(i,j,k) * v(i,j,k+1) + w(i,j,k) * w(i,j,k+1)) &
          + 0.5d0 * (w(i,j,k) * p(i,j,k+1) + w(i,j,k+1) * p(i,j,k))
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_G2

  subroutine calc_G4(id, nx, ny, nz, gamma, rho, u, v, w, p, G)
    integer(kind=4), intent(in), value :: id
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8), intent(out) :: G(nx-4,ny-4,nz-3,5)
    integer i, j, k
    real(8), dimension(3) :: RhoW, RhoWW_P, IE, Energy
    do k = 1, nz-3
      do j = 3, ny-2
        do i = 3, nx-2
          RhoW(:) = RhoPhi(rho(i,j,k:k+3), w(i,j,k:k+3))
          RhoWW_P(:) = RhoPhiU(RhoW(:), w(i,j,k:k+3)) + Phi(p(i,j,k:k+3))
          IE(:) = p(i,j,k:k+3) / ((gamma - 1.d0) * rho(i,j,k:k+3))
          Energy(:) = IE(:) + RhoUPhiPhi(RhoW(:), u(i,j,k:k+3), v(i,j,k:k+3), w(i,j,k:k+3)) + PhiPsi(w(i,j,k:k+3), p(i,j,k:k+3))
          G(i-2,j-2,k,1) = Flux(RhoW(:))
          G(i-2,j-2,k,2) = Flux(RhoPhiU(RhoW(:), u(i,j,k:k+3)))
          G(i-2,j-2,k,3) = Flux(RhoPhiU(RhoW(:), v(i,j,k:k+3)))
          G(i-2,j-2,k,4) = Flux(RhoWW_P(:))
          G(i-2,j-2,k,5) = Flux(Energy(:))
        enddo
      enddo
    enddo
  end subroutine calc_G4
end module calc_flux