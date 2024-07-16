module calc_common_dim
  implicit none
contains
  attributes(device) function vecsum(V1,V2) result(sum)
    real(8), intent(in), dimension(1), device :: V1, V2
    real(8) sum
    sum = V1(1) * V2(1)
  end function vecsum

  attributes(device) function q2(Vl,Vr) result(ans)
    real(8), intent(in), dimension(1), device :: Vl, Vr
    real(8) ans
    ans = Vl(1)**2 + Vl(2)**2
  end function q2
end module calc_common_dim

