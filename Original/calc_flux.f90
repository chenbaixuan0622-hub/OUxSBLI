module calc_flux
  use calc_term
  implicit none
contains
  subroutine calc_E(na, nx, ny, nz, gamma, rho, u, v, w, p, E)
    integer, intent(in) :: na, nx, ny, nz
    real(8), intent(in) :: gamma
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8), intent(out), dimension(nx-na+1,ny-na,nz-na,5) :: E
    integer i, j, k, offset
    real(8), dimension(na-1) :: P_keep, RhoU, UP
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    offset = na / 2
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1, nx-na+1
          call Phi(p(i:i+na-1,j,k),P_keep)
          call RhoPhi(rho(i:i+na-1,j,k), u(i:i+na-1,j,k), RhoU)
          call PhiPsi(u(i:i+na-1,j,k), p(i:i+na-1,j,k), UP)
          E(i,j-offset,k-offset,1) = Flux(RhoU)
          E(i,j-offset,k-offset,2) = Flux(RhoPhiU(RhoU, u(i:i+na-1,j,k)) + P_keep)
          E(i,j-offset,k-offset,3) = Flux(RhoPhiU(RhoU, v(i:i+na-1,j,k)))
          E(i,j-offset,k-offset,4) = Flux(RhoPhiU(RhoU, w(i:i+na-1,j,k)))
          E(i,j-offset,k-offset,5) = Flux(RhoPhiU(RhoU, p(i:i+na-1,j,k) / rho(i:i+na-1,j,k)) / (gamma - 1.d0)&
          + RhoUPhiPhi(RhoU, u(i:i+na-1,j,k), v(i:i+na-1,j,k), w(i:i+na-1,j,k)) + UP)
        enddo
      enddo
    enddo
  end subroutine calc_E

  subroutine calc_F(na, nx, ny, nz, gamma, rho, u, v, w, p, F)
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
  end subroutine calc_F

  subroutine calc_G(na, nx, ny, nz, gamma, rho, u, v, w, p, G)
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
  end subroutine calc_G
end module calc_flux