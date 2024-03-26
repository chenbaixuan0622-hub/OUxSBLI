module calc_mat
  use mod_globals, only : gamma
  use calc_common_dim
  implicit none
contains
  attributes(device) function calc_AB(id,rho,H,c,V) result(mat)
    integer, intent(in), value :: id
    real(8), intent(in), value :: rho, H, c
    real(8), intent(in), device :: V(3)
    real(8), dimension(5,5) :: mat, lambda, R, Rinv, Rlambda
    real(8) q2, b1, b2
    real(8) :: vec(3) = 0.d0
    real(8) d(3,3)
    integer id1, id2, id3
    !!!!!!!!!!!!!!!!!
    ! id  ! 1  2  3 !  
    !!!!!!!!!!!!!!!!!
    ! id1 ! 3  1  2 !
    ! id2 ! 1  2  3 !
    ! id3 ! 2  3  1 !
    !!!!!!!!!!!!!!!!!
    id1 = mod(id+1,3) + 1
    id2 = mod(id+2,3) + 1
    id3 = mod(id+3,3) + 1
    vec(id) = 1.d0
    d(1,:) = (/1.d0, 0.d0, 0.d0/)
    d(2,:) = (/0.d0, 1.d0, 0.d0/)
    d(3,:) = (/0.d0, 0.d0, 1.d0/)
    q2 = V(1)**2 + V(2)**2 + V(3)**2
    b2 = (gamma - 1.d0) / c**2
    b1 = 0.5d0 * q2 * b2

    lambda(1,:) = (/abs(V(id)-c), 0.d0, 0.d0, 0.d0, 0.d0/)
    lambda(2,:) = (/0.d0, abs(V(id)), 0.d0, 0.d0, 0.d0/)
    lambda(3,:) = (/0.d0, 0.d0, abs(V(id)), 0.d0, 0.d0/)
    lambda(4,:) = (/0.d0, 0.d0, 0.d0, abs(V(id)), 0.d0/)
    lambda(5,:) = (/0.d0, 0.d0, 0.d0, 0.d0, abs(V(id)+c)/)

    R(:,1) = (/1.d0, V(1) - c * vec(1), V(2) - c * vec(2), V(3) - c * vec(3), H - c * V(id)/)
    R(:,5) = (/1.d0, V(1) + c * vec(1), V(2) + c * vec(2), V(3) + c * vec(3), H + c * V(id)/)
    R(:,1+id1) = (/1.d0, V(1), V(2), V(3), 0.5d0 * q2/)
    R(:,1+id2) = (/0.d0, - c * d(1,id1), - c * d(2,id1), - c * d(3,id1), - c * V(id1)/)
    R(:,1+id3) = (/0.d0,   c * d(1,id3),   c * d(2,id3),   c * d(3,id3),   c * V(id3)/)

    Rinv(1,:) = (/0.5d0*(b1+V(id)/c), -0.5d0*(d(1,id)/c+b2*V(1)), -0.5d0*(d(2,id)/c+b2*V(2)), -0.5d0*(d(3,id)/c+b2*V(3)), 0.5d0*b2/)
    Rinv(1+id2,:) = (/V(id1)/c,-d(1,id1)/c,-d(2,id1)/c,-d(3,id1)/c, 0.d0/)
    Rinv(1+id3,:) = (/-V(id3)/c,d(1,id3)/c,d(2,id3)/c,d(3,id3)/c, 0.d0/)
    Rinv(1+id1,:) = (/1.d0 - b1, b2 * V(1), b2 * V(2), b2 * V(3), -b2/)
    Rinv(5,:) = (/0.5d0*(b1-V(id)/c), -0.5d0*(-d(1,id)/c+b2*V(1)), -0.5d0*(-d(2,id)/c+b2*V(2)), -0.5d0*(-d(3,id)/c+b2*V(3)), 0.5d0*b2/)

    Rlambda(:,:) = cumatmul(R(:,:),lambda(:,:))
    mat(:,:) = cumatmul(Rlambda(:,:),Rinv(:,:))
  end function calc_AB
end module calc_mat

