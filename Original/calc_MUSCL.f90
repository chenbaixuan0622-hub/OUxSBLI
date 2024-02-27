module calc_MUSCL
  use calc_physical_quantities
  implicit none
contains
  function minmod(x, y) result(ans)
    real(8), intent(in), value :: x, y
    real(8) sgn
    real(8) ans
    sgn = sign(1.d0,x)
    ans = sgn * max(min(abs(x), sgn * y), 0.d0)
  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine MUSCLx(nx,ny,nz,k_cof,b,Q,Ql,Qr)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: k_cof, b
    real(8), intent(in) :: Q(nx,ny,nz,5)
    real(8), intent(out), dimension(nx-1,ny,nz,5) :: Ql, Qr
    real(8) d1, d2, d3, dm1, dp1, dm2, dp2
    integer i, j, k, l
    do l = 1, 5
      do k = 1, nz
        do j = 1, ny
          do i = 2, nx-1
            d1 = Q(i,j,k,l) - Q(i-1,j,k,l)
            d2 = Q(i+1,j,k,l) - Q(i,j,k,l)
            d3 = Q(i+2,j,k,l) - Q(i+1,j,k,l)
            dp1 = minmod(d2, b * d1)
            dm1 = minmod(d1, b * d2)
            dp2 = minmod(d3, b * d2)
            dm2 = minmod(d2, b * d3)
            Ql(i,j,k,l) = Q(i,j,k,l) + 0.25d0 * ((1.d0 - k_cof) * dm1 + (1.d0 + k_cof) * dp1)
            Qr(i,j,k,l) = Q(i+1,j,k,l) - 0.25d0 * ((1.d0 - k_cof) * dp2 + (1.d0 + k_cof) * dm2)
          enddo
        enddo
      enddo
    enddo
  end subroutine MUSCLx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module
  