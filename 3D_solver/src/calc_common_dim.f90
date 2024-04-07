module calc_common_dim
  implicit none
  interface cumatmul
    module procedure cumatmul55, cumatmul51
  end interface
contains
  attributes(device) function vecsum(V1,V2) result(sum)
    real(8), intent(in), dimension(3), device :: V1, V2
    real(8) sum
    sum = V1(1) * V2(1) + V1(2) * V2(2) + V1(3) * V2(3)
  end function vecsum

  attributes(device) function q2(Vl,Vr) result(ans)
    real(8), intent(in), dimension(3), device :: Vl, Vr
    real(8) ans
    ans = Vl(1)**2 + Vl(2)**2 + Vl(3)**2 + Vr(1)**2 + Vr(2)**2 + Vr(3)**2
  end function q2
  
  attributes(device) function Cij(A,B) result(ans)
    real(8), intent(in), dimension(5), device :: A, B
    real(8) ans
    ans = A(1) * B(1) + A(2) * B(2) + A(3) * B(3) + A(4) * B(4) + A(5) * B(5)
  end function

  attributes(device) function cumatmul55(A,B) result(C)
    real(8), intent(in), dimension(5,5), device :: A, B
    real(8), dimension(5,5) :: C
    real(8), dimension(5) :: A1, A2, A3, A4, A5, B1, B2, B3, B4, B5
    A1(:) = A(1,:)
    A2(:) = A(2,:)
    A3(:) = A(3,:)
    A4(:) = A(4,:)
    A5(:) = A(5,:)
    B1(:) = B(:,1)
    B2(:) = B(:,2)
    B3(:) = B(:,3)
    B4(:) = B(:,4)
    B5(:) = B(:,5)
    C(1,:) = (/Cij(A1,B1), Cij(A1,B2), Cij(A1,B3), Cij(A1,B4), Cij(A1,B5)/)
    C(2,:) = (/Cij(A2,B1), Cij(A2,B2), Cij(A2,B3), Cij(A2,B4), Cij(A2,B5)/)
    C(3,:) = (/Cij(A3,B1), Cij(A3,B2), Cij(A3,B3), Cij(A3,B4), Cij(A3,B5)/)
    C(4,:) = (/Cij(A4,B1), Cij(A4,B2), Cij(A4,B3), Cij(A4,B4), Cij(A4,B5)/)
    C(5,:) = (/Cij(A5,B1), Cij(A5,B2), Cij(A5,B3), Cij(A5,B4), Cij(A5,B5)/)
  end function cumatmul55
  
  attributes(device) function cumatmul51(A,B) result(C)
    real(8), intent(in), dimension(5,5), device :: A
    real(8), intent(in), dimension(5), device :: B
    real(8), dimension(5) :: C
    real(8), dimension(5) :: A1, A2, A3, A4, A5
    A1(:) = A(1,:)
    A2(:) = A(2,:)
    A3(:) = A(3,:)
    A4(:) = A(4,:)
    A5(:) = A(5,:)
    C(1) = Cij(A1,B)
    C(2) = Cij(A2,B)
    C(3) = Cij(A3,B)
    C(4) = Cij(A4,B)
    C(5) = Cij(A5,B)
  end function cumatmul51
end module calc_common_dim

