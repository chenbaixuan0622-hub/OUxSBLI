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
  subroutine calc_E2(na, nx, ny, nz, gamma, rho, u, v, w, p, E)
    integer, intent(in) :: na, nx, ny, nz
    real(8), intent(in) :: gamma
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8), intent(out), dimension(nx-1,ny-2,nz-2,5) :: E
    integer i, j, k
    real(8) :: Rho_m, U_m, P_m
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
  end subroutine calc_E2

  subroutine calc_E4(na, nx, ny, nz, gamma, rho, u, v, w, p, E)
    integer, intent(in) :: na, nx, ny, nz
    real(8), intent(in) :: gamma
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8), intent(out), dimension(nx-3,ny-4,nz-4,5) :: E
    integer i, j, k, offset
    real(8), dimension(3) :: P_keep, RhoU, UP
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    offset = na / 2
    do k = 3, nz-2
      do j = 3, ny-2
        do i = 1, nx-3
          call Phi(p(i:i+3,j,k),P_keep)
          call RhoPhi(rho(i:i+3,j,k), u(i:i+3,j,k), RhoU)
          call PhiPsi(u(i:i+3,j,k), p(i:i+3,j,k), UP)
          E(i,j-2,k-2,1) = Flux(RhoU)
          E(i,j-2,k-2,2) = Flux(RhoPhiU(RhoU, u(i:i+3,j,k)) + P_keep)
          E(i,j-2,k-2,3) = Flux(RhoPhiU(RhoU, v(i:i+3,j,k)))
          E(i,j-2,k-2,4) = Flux(RhoPhiU(RhoU, w(i:i+3,j,k)))
          E(i,j-2,k-2,5) = Flux(RhoPhiU(RhoU, p(i:i+3,j,k) / rho(i:i+3,j,k)) / (gamma - 1.d0)&
          + RhoUPhiPhi(RhoU, u(i:i+3,j,k), v(i:i+3,j,k), w(i:i+3,j,k)) + UP)
        enddo
      enddo
    enddo
  end subroutine calc_E4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_F2(nx, ny, nz, gamma, rho, u, v, w, p, F)
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in) :: gamma
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8), intent(out) :: F(nx-2,ny-1,nz-2,5)
    integer i, j, k
    real(8) Rho_m, V_m, P_m
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
  end subroutine calc_F2

  subroutine calc_F4(na, nx, ny, nz, gamma, rho, u, v, w, p, F)
    integer, intent(in) :: na, nx, ny, nz
    real(8), intent(in) :: gamma
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8), intent(out) :: F(nx-na,ny-na+1,nz-na,5)
    integer i, j, k, offset
    real(8), dimension(na-1) :: P_keep, RhoV, VP
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    offset = na / 2
    do k = 1+offset, nz-offset
      do j = 1, ny-na+1
        do i = 1+offset, nx-offset
          call Phi(p(i,j:j+na-1,k), P_keep)
          call RhoPhi(rho(i,j:j+na-1,k), v(i,j:j+na-1,k), RhoV)
          call PhiPsi(v(i,j:j+na-1,k), p(i,j:j+na-1,k), VP)
          F(i-offset,j,k-offset,1) = Flux(RhoV)
          F(i-offset,j,k-offset,2) = Flux(RhoPhiU(RhoV, u(i,j:j+na-1,k)))
          F(i-offset,j,k-offset,3) = Flux(RhoPhiU(RhoV, v(i,j:j+na-1,k)) + P_keep)
          F(i-offset,j,k-offset,4) = Flux(RhoPhiU(RhoV, w(i,j:j+na-1,k)))
          F(i-offset,j,k-offset,5) = Flux(RhoPhiU(RhoV, p(i,j:j+na-1,k) / rho(i,j:j+na-1,k)) / (gamma - 1.0d0)&
          + RhoUPhiPhi(RhoV, u(i,j:j+na-1,k), v(i,j:j+na-1,k), w(i,j:j+na-1,k)) + VP)
        enddo
      enddo
    enddo
  end subroutine calc_F4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_G2(nx, ny, nz, gamma, rho, u, v, w, p, G)
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in) :: gamma
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8), intent(out) :: G(nx-2,ny-2,nz-1,5)
    integer i, j, k
    real(8) Rho_m, W_m, P_m
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
  end subroutine calc_G2

  subroutine calc_G4(na, nx, ny, nz, gamma, rho, u, v, w, p, G)
    integer, intent(in) :: na, nx, ny, nz
    real(8), intent(in) :: gamma
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8), intent(out) :: G(nx-na,ny-na,nz-na+1,5)
    integer i, j, k, offset
    real(8), dimension(na-1) :: P_keep, RhoW, WP
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    offset = na / 2
    do k = 1, nz-na+1
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          call Phi(p(i,j,k:k+na-1), P_keep)
          call RhoPhi(rho(i,j,k:k+na-1), w(i,j,k:k+na-1), RhoW)
          call PhiPsi(w(i,j,k:k+na-1), p(i,j,k:k+na-1), WP)
          G(i-offset,j-offset,k,1) = Flux(RhoW)
          G(i-offset,j-offset,k,2) = Flux(RhoPhiU(RhoW, u(i,j,k:k+na-1)))
          G(i-offset,j-offset,k,3) = Flux(RhoPhiU(RhoW, v(i,j,k:k+na-1)))
          G(i-offset,j-offset,k,4) = Flux(RhoPhiU(RhoW, w(i,j,k:k+na-1)) + P_keep)
          G(i-offset,j-offset,k,5) = Flux(RhoPhiU(RhoW, p(i,j,k:k+na-1) / rho(i,j,k:k+na-1)) / (gamma - 1.d0)&
          + RhoUPhiPhi(RhoW, u(i,j,k:k+na-1), v(i,j,k:k+na-1), w(i,j,k:k+na-1)) + WP)
        enddo
      enddo
    enddo
  end subroutine calc_G4
end module calc_flux