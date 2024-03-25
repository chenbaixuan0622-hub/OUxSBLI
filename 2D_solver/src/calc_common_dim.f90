module calc_common_dim
  implicit none
  interface cumatmul
    module procedure cumatmul44, cumatmul41
  end interface
contains
  attributes(device) function vecsum(V1,V2) result(sum)
    real(8), intent(in), dimension(2), device :: V1, V2
    real(8) sum
    sum = V1(1) * V2(1) + V1(2) * V2(2)
  end function vecsum
  
  attributes(device) function Cij(A,B) result(ans)
    real(8), intent(in), dimension(4), device :: A, B
    real(8) ans
    ans = A(1) * B(1) + A(2) * B(2) + A(3) * B(3) + A(4) * B(4)
  end function

  attributes(device) function cumatmul44(A,B) result(C)
    real(8), intent(in), dimension(4,4), device :: A, B
    real(8), dimension(4,4) :: C
    real(8), dimension(4) :: A1, A2, A3, A4, B1, B2, B3, B4
    A1(:) = A(1,:)
    A2(:) = A(2,:)
    A3(:) = A(3,:)
    A4(:) = A(4,:)
    B1(:) = B(:,1)
    B2(:) = B(:,2)
    B3(:) = B(:,3)
    B4(:) = B(:,4)
    C(1,:) = (/Cij(A1,B1), Cij(A1,B2), Cij(A1,B3), Cij(A1,B4)/)
    C(2,:) = (/Cij(A2,B1), Cij(A2,B2), Cij(A2,B3), Cij(A2,B4)/)
    C(3,:) = (/Cij(A3,B1), Cij(A3,B2), Cij(A3,B3), Cij(A3,B4)/)
    C(4,:) = (/Cij(A4,B1), Cij(A4,B2), Cij(A4,B3), Cij(A4,B4)/)
  end function cumatmul44
  
  attributes(device) function cumatmul41(A,B) result(C)
    real(8), intent(in), dimension(4,4), device :: A
    real(8), intent(in), dimension(4), device :: B
    real(8), dimension(4) :: C
    real(8), dimension(4) :: A1, A2, A3, A4
    A1(:) = A(1,:)
    A2(:) = A(2,:)
    A3(:) = A(3,:)
    A4(:) = A(4,:)
    C(1) = Cij(A1,B)
    C(2) = Cij(A2,B)
    C(3) = Cij(A3,B)
    C(4) = Cij(A4,B)
  end function cumatmul41
end module calc_common_dim

