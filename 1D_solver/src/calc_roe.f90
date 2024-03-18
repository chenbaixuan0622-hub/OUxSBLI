module calc_roe
  use calc_Qlr
  use calc_common
  use calc_mat
  implicit none
contains
  subroutine calc_E(nx,gamma,k,b,rho,u,p,E)
    integer, intent(in) :: nx
    real(8), intent(in) :: gamma, k, b
    real(8), intent(in), dimension(nx) :: rho, u, p
    real(8), intent(out), dimension(nx-1,3) :: E
    integer i
    real(8) rhol, rhor, ul, ur, pl, pr, energyl, energyr, Hl, Hr, c, cl, cr, rho_ave, u_ave, H_ave, c_ave
    real(8), dimension(3) :: d1, d2, d3, Ql, Qr, El, Er
    real(8), dimension(3,3) :: A
    real(8), dimension(nx,3) :: Q
    Q(:,1) = rho
    Q(:,2) = u
    Q(:,3) = p
    do i = 1, nx-1
      if (2 <= i .and. i <= nx-2) then
        call Qlr(k,b,Q(i-1,:),Q(i,:),Q(i+1,:),Q(i+2,:),rhol,rhor,ul,ur,pl,pr)
      elseif (i == 1) then
        call Qlr(k,b,0.d0,Q(i,:),Q(i+1,:),Q(i+2,:),rhol,rhor,ul,ur,pl,pr)
      else
        call Qlr(k,b,Q(i-1,:),Q(i,:),Q(i+1,:),0.d0,rhol,rhor,ul,ur,pl,pr)
      endif
      
      energyl = energy(gamma,pl,rhol,ul)
      energyr = energy(gamma,pr,rhor,ur)

      Hl = (energyl + pl) / rhol
      Hr = (energyr + pr) / rhor
      rho_ave = sqrt(rhol * rhor)
      u_ave = (sqrt(rhol) * ul + sqrt(rhor) * ur) / (sqrt(rhol) + sqrt(rhor))
      H_ave = (sqrt(rhol) * Hl + sqrt(rhor) * Hr) / (sqrt(rhol) + sqrt(rhor))
      c_ave = sqrt((gamma - 1.d0) * (H_ave - 0.5d0 * u_ave**2))
      
      call calc_A(gamma,rho_ave,u_ave,H_ave,c_ave,A)

      Ql(1) = rhol
      Qr(1) = rhor
      Ql(2) = rhol * ul
      Qr(2) = rhor * ur
      Ql(3) = energyl
      Qr(3) = energyr
      El(1) = Ql(1) 
      Er(1) = Qr(2)
      El(2) = ul * Ql(1) + pl
      Er(2) = ur * Qr(1) + pr
      El(3) = (energyl + pl) * ul
      Er(3) = (energyr + pr) * ur
      ! calc Roe
      E(i,:) = 0.5d0 * (El(:) + Er(:) - matmul(A,Qr(:) - Ql(:)))
    enddo
  end subroutine calc_E
end module calc_roe

