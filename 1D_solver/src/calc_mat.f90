module calc_mat
  implicit none
contains
  subroutine calc_A(gamma,rho,u,H,c,A)
    real(8), intent(in) :: gamma, rho, u, H, c
    real(8), intent(out) :: A(3,3)
    real(8) b1, b2 
    real(8), dimension(3,3) :: R, Rinv, lambda
    b2 = (gamma - 1.d0) / c**2
    b1 = 0.5d0 * u**2 * b2
    lambda(1,:) = (/abs(u - c),   0.d0,       0.d0/)
    lambda(2,:) = (/0.d0,       abs(u),       0.d0/)
    lambda(3,:) = (/0.d0,         0.d0, abs(u + c)/)
    R(1,:) = (/1.d0,              1.d0,      1.d0/)
    R(2,:) = (/u - c,                u,     u + c/)
    R(3,:) = (/H - u * c, 0.5d0 * u**2, H + u * c/)
    Rinv(1,:) = (/0.5d0 * (b1 + u / c), -0.5d0 * (1.d0 / c + b2 * u), 0.5d0 * b2/)
    Rinv(2,:) = (/1.d0 - b1,                                  b2 * u,        -b2/)
    Rinv(3,:) = (/0.5d0 * (b1 - u / c),  0.5d0 * (1.d0 / c - b2 * u), 0.5d0 * b2/)
    A = matmul(matmul(R,lambda),Rinv)
  end subroutine calc_A
end module calc_mat

