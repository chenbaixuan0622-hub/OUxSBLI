module calc_mat
  use calc_common_dim
  use mod_globals, only : gamma
  implicit none
contains
  attributes(device) function calc_AB(id,rho,H,c,V) result(mat)
    integer, intent(in), value :: id
    real(8), intent(in), value :: rho, H, c
    real(8), intent(in), device :: V(2)
    real(8), dimension(4,4) :: mat, lambda, R, Rinv, Rlambda
    real(8) q2, b1, b2
    real(8) :: vec(2) = 0.d0
    real(8), dimension(4) :: R2, R3
    integer idi
    idi = mod(id,2) + 1 
    vec(id) = 1.d0
    q2 = V(1)**2 + V(2)**2
    b2 = (gamma - 1.d0) / c**2
    b1 = 0.5d0 * q2 * b2

    lambda(1,:) = (/abs(V(id)-c), 0.d0, 0.d0, 0.d0/)
    lambda(2,:) = (/0.d0, abs(V(id)), 0.d0, 0.d0/)
    lambda(3,:) = (/0.d0, 0.d0, abs(V(id)+c), 0.d0/)
    lambda(4,:) = (/0.d0, 0.d0, 0.d0, abs(V(id))/)
    
    R(1,:) = (/1.d0,                    1.d0,              1.d0,   0.d0/)
    R(2,:) = (/V(1) - vec(1) * c,       V(1), V(1) + vec(1) * c, vec(2)/)
    R(3,:) = (/V(2) - vec(2) * c,       V(2), V(2) + vec(2) * c, vec(1)/)
    R(4,:) = (/H - c * V(id),     0.5d0 * q2,     H + c * V(id), V(idi)/)
    
    Rinv(1,:) = (/0.5d0 * (b1 + V(id) / c), -0.5d0 * (1.d0 / c + b2 * V(id)), -0.5d0 * b2 * V(idi), 0.5d0 * b2/)
    Rinv(2,:) = (/1.d0 - b1,                                      b2 * V(id),          b2 * V(idi),        -b2/)
    Rinv(3,:) = (/0.5d0 * (b1 - V(id) / c),  0.5d0 * (1.d0 / c - b2 * V(id)), -0.5d0 * b2 * V(idi), 0.5d0 * b2/)
    Rinv(4,:) = (/-V(idi),                                              0.d0,                 1.d0,       0.d0/)
    
    ! swap columns of arrays
    R2 = Rinv(:,1+id)
    R3 = Rinv(:,1+idi)
    Rinv(:,2) = R2
    Rinv(:,3) = R3
    Rlambda(:,:) = cumatmul(R(:,:), lambda(:,:))
    mat(:,:) = cumatmul(Rlambda(:,:),Rinv(:,:))
  end function calc_AB
end module calc_mat

