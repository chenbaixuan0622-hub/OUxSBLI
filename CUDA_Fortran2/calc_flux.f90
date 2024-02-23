module calc_flux
  use calc_term
  implicit none
contains
  attributes(global) subroutine calcE(nx, ny, nz, gamma, rho, u, v, w, p, E)
    integer, value :: nx, ny, nz
    real(8), value :: gamma
    real(8), dimension(nx,ny,nz,4) :: rho, u, v, w, p
    real(8) E(nx-3,ny-4,nz-4,5,4)
    integer :: i, j, k
    k = 3 + nk * mod(int(threadIdx%x/64),8) !k = 3, nz-nk-1, nk :3 <= k <= nz-2
    j = 3 + nj * mod(int(threadIdx%x/8),8)  !j = 3, ny-nj-1, nj :3 <= j <= ny-2
    i = 3 + ni * mod(threadIdx%x,8)         !i = 3, nx-ni-1, ni :3 <= i <= nx-2
    !real(8), dimension(3) :: rho_d, u_d, v_d, w_d, p_d 
    !real(8), dimension(3) :: P_keep, RhoU, RhoUU, RhoUV, RhoUW, IE, RhoUIE, RhoUKE, UP, RhoUU_P, E_total
    !do k = 3, nz-2
    !  do j = 3, ny-2
    !    do i = 1, nx-3
    !      rho_d = rho(i:i+3,j,k)
    !      u_d = u(i:i+3,j,k)
    !      v_d = v(i:i+3,j,k)
    !      w_d = w(i:i+3,j,k)
    !      p_d = p(i:i+3,j,k)
    !      call calc_Phi(p_d, P_keep)
    !      call calc_RhoPhi(rho_d, u_d, RhoU)
    !      call calc_RhoPhiU(RhoU, u_d, RhoUU)
    !      call calc_RhoPhiU(RhoU, v_d, RhoUV)
    !      call calc_RhoPhiU(RhoU, w_d, RhoUW)
    !      IE = p_d / ((gamma - 1.d0) * rho_d)
    !      call calc_RhoPhiU(RhoU, IE, RhoUIE)
    !      call calc_RhoUPhi2(RhoU, u_d, v_d, w_d, RhoUKE)
    !      call calc_PhiPsi(u_d, p_d, UP)
    !      call merge_term(RhoU, E(i,j,k,1))
    !      RhoUU_P = RhoUU + P_keep
    !      call merge_term(RhoUU_P, E(i,j,k,2))
    !      call merge_term(RhoUV, E(i,j,k,3))
    !      call merge_term(RhoUW, E(i,j,k,4))
    !      E_total = RhoUIE + RhoUKE + UP
    !      call merge_term(E_total, E(i,j,k,5))
    !    enddo
    !  enddo
    !enddo
  end subroutine calcE
  
  attributes(global) subroutine calcF(nx, ny, nz, gamma, rho, u, v, w, p, F)
    integer, value :: nx, ny, nz
    real(8), value :: gamma
    real(8), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8) F(nx,ny-3,nz,5)
    integer :: i, j, k
    real(8), dimension(3) :: rho_d, u_d, v_d, w_d, p_d 
    real(8), dimension(3) :: P_keep, RhoV, RhoVU, RhoVV, RhoVW, IE, RhoVIE, RhoVKE, VP, RhoVV_P, E_total
    do k = 3, nz-2
      do j = 1, ny-3
        do i = 3, nx-2
          rho_d = rho(i,j:j+3,k)
          u_d = u(i,j:j+3,k)
          v_d = v(i,j:j+3,k)
          w_d = w(i,j:j+3,k)
          p_d = p(i,j:j+3,k)
          call calc_Phi(p_d, P_keep)
          call calc_RhoPhi(rho_d, v_d, RhoV)
          call calc_RhoPhiU(RhoV, u_d, RhoVU)
          call calc_RhoPhiU(RhoV, v_d, RhoVV)
          call calc_RhoPhiU(RhoV, w_d, RhoVW)
          IE = p_d / ((gamma - 1.0d0) * rho_d)
          call calc_RhoPhiU(RhoV, IE, RhoVIE)
          call calc_RhoUPhi2(RhoV, u_d, v_d, w_d, RhoVKE)
          call calc_PhiPsi(v_d, p_d, VP)
          call merge_term(RhoV, F(i,j,k,1))
          call merge_term(RhoVU, F(i,j,k,2))
          RhoVV_P = RhoVV + P_keep
          call merge_term(RhoVV_P, F(i,j,k,3))
          call merge_term(RHoVW, F(i,j,k,4))
          E_total = RhoVIE + RhoVKE + VP
          call merge_term(E_total, F(i,j,k,5))
        enddo
      enddo
    enddo
  end subroutine calcF

  attributes(global) subroutine calcG(nx, ny, nz, gamma, rho, u, v, w, p, G)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8), intent(out) :: G(nx,ny,nz-3,5)
    integer :: i, j, k
    real(8), dimension(3) :: rho_d, u_d, v_d, w_d, p_d 
    real(8), dimension(3) :: P_keep, RhoW, RhoWU, RhoWV, RhoWW, IE, RhoWIE, RhoWKE, WP, RhoWW_P, E_total
    do k = 1, nz-3
      do j = 3, ny-2
        do i = 3, nx-2
          rho_d = rho(i,j,k:k+3)
          u_d = u(i,j,k:k+3)
          v_d = v(i,j,k:k+3)
          w_d = w(i,j,k:k+3)
          p_d = p(i,j,k:k+3)
          call calc_Phi(p_d, P_keep)
          call calc_RhoPhi(rho_d, w_d, RhoW)
          call calc_RhoPhiU(RhoW, u_d, RhoWU)
          call calc_RhoPhiU(RhoW, v_d, RhoWV)
          call calc_RhoPhiU(RhoW, w_d, RhoWW)
          IE = p_d / ((gamma - 1.d0) * rho_d)
          call calc_RhoPhiU(RhoW, IE, RhoWIE)
          call calc_RhoUPhi2(RhoW, u_d, v_d, w_d, RhoWKE)
          call calc_PhiPsi(w_d, p_d, WP)
          call merge_term(RhoW, G(i,j,k,1))
          call merge_term(RhoWU, G(i,j,k,2))
          call merge_term(RhoWV, G(i,j,k,3))
          RhoWW_P = RhoWW + P_keep
          call merge_term(RhoWW_P, G(i,j,k,4))
          E_total = RhoWIE + RhoWKE + WP
          call merge_term(E_total, G(i,j,k,5))
        enddo
      enddo
    enddo
  end subroutine calcG
end module calc_flux
