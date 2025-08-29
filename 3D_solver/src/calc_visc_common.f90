module calc_visc_common
  use mod_constant, only : Cp, Cp_over_Pr
  implicit none
contains
  !dir$ inline
  attributes(device) function interpolation6(a) result(ans)
    real(8), intent(in), device :: a(6)
    real(8) ans(3)
    ans(:) = 0.0625d0 * (-a(1:3) + 9.d0 * (a(2:4) + a(3:5)) -a(4:6))
  end function interpolation6

  !dir$ inline
  attributes(device) function dx6(a, dx) result(ans)
    real(8), intent(in), device :: a(6)
    real(8), intent(in), value  :: dx
    real(8) ans(3)
    ans(:) = 0.125d0 * (9.d0 * (-a(2:4) + a(3:5)) - (-a(1:3) + a(4:6)) / 3.d0) * dx
  end function dx6

  !dir$ inline
  attributes(device) function dy5(a, dy) result(ans)
    real(8), intent(in), device :: a(5)
    real(8), intent(in), value  :: dy
    real(8) ans
    ans = (2.d0 * (-a(2) + a(4)) - 0.25d0 * (-a(1) + a(5))) * dy / 3.d0
  end function dy5

  !dir$ inline
  attributes(device) function dy23(mu, a, dy) result(ans)
    real(8), intent(in), device :: mu(2), a(2,3)
    real(8), intent(in), value  :: dy
    real(8) ans
    ans = 0.25d0 * (mu(1) * (-a(1,1) + a(1,2) -a(2,1) + a(2,2)) &
                  + mu(2) * (-a(1,2) + a(1,3) -a(2,2) + a(2,3))) * dy
  end function dy23

  attributes(device) function dy65(a, dy) result(ans)
    real(8), intent(in), device :: a(6,5)
    real(8), intent(in), value  :: dy
    real(8) ans(3), a1(5), a2(5), a3(5), a4(5), a5(5), a6(5), ay(6)
    a1 = a(1,:); a2 = a(2,:); a3 = a(3,:); a4 = a(4,:); a5 = a(5,:); a6 = a(6,:)
    ay(1) = dy5(a1(:), dy); ay(2) = dy5(a2(:), dy); ay(3) = dy5(a3(:), dy)
    ay(4) = dy5(a4(:), dy); ay(5) = dy5(a5(:), dy); ay(6) = dy5(a6(:), dy)
    ans = interpolation6(ay(:))
  end function dy65
  
  !dir$ inline
  attributes(device) function dy32(mu, a, dy) result(ans)
    real(8), intent(in), device :: mu(2), a(3,2)
    real(8), intent(in), value  :: dy
    real(8) ans
    ans = 0.25d0 * (mu(1) * (-a(1,1) + a(2,1) - a(1,2) + a(2,2)) &
                  + mu(2) * (-a(2,1) + a(3,1) - a(2,2) + a(3,2))) * dy
  end function dy32

  attributes(device) function dy56(a, dy) result(ans)
    real(8), intent(in), device :: a(5,6)
    real(8), intent(in), value  :: dy
    real(8) ans(3), ay(6)
    ay(1) = dy5(a(:,1), dy); ay(2) = dy5(a(:,2), dy); ay(3) = dy5(a(:,3), dy)
    ay(4) = dy5(a(:,4), dy); ay(5) = dy5(a(:,5), dy); ay(6) = dy5(a(:,6), dy)
    ans = interpolation6(ay(:))
  end function dy56

  !dir$ inline
  attributes(device) function flux4(a) result(ans)
    real(8), intent(in), device :: a(3)
    real(8) ans
    ans = 0.125d0 * ((9.d0 - 1.d0 / 3.d0) * a(2) - (a(1) + a(3)) / 3.d0)
  end function flux4

  attributes(device) subroutine tauxx4(mu, ux, vy, wz, u6, txx, utxx)
    real(8), intent(in), dimension(3), device :: mu, ux, vy, wz
    real(8), intent(in), dimension(6), device :: u6
    real(8), intent(out)  :: txx, utxx
    real(8), dimension(3) :: tmp
    tmp(:) = 2.d0 * mu(:) * (2.d0 * ux(:) - vy(:) - wz(:)) / 3.d0 ! tau
    txx    = flux4(tmp(:))
    tmp(:) = interpolation6(u6(:)) * tmp(:) ! u tau
    utxx   = flux4(tmp(:))
  end subroutine tauxx4

  attributes(device) subroutine tauxy4(mu, uy, vx, v6, txy, vtxy)
    real(8), intent(in), dimension(3), device :: mu, uy, vx
    real(8), intent(in), dimension(6), device :: v6
    real(8), intent(out)  :: txy, vtxy
    real(8), dimension(3) :: tmp
    tmp(:) = mu(:) * (uy(:) + vx(:)) ! tau
    txy    = flux4(tmp(:))
    tmp(:) = interpolation6(v6(:)) * tmp(:) ! vtau
    vtxy   = flux4(tmp(:))
  end subroutine tauxy4

  !dir$ inline
  attributes(device) function heat_conduction6(mu, T, dx) result(ans)
    real(8), intent(in), device :: mu(3), T(6)
    real(8), intent(in), value  :: dx
    real(8), dimension(3) ::  kTx
    real(8) :: ans
    kTx(:) = Cp_over_Pr * mu(:) * dx6(T(:), dx)
    ans = flux4(kTx(:))
  end function heat_conduction6
end module calc_visc_common

