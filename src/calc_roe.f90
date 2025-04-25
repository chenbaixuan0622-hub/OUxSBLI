module calc_roe
  use mod_globals, only : dimension, gamma
  use calc_common_dim
  implicit none
contains
  attributes(device) subroutine calc_A(u, v, w, c, q2, H, b1, b2, lambda, mat)
    real(8), intent(in)             :: u, v, w, c, q2, H, b1, b2
    real(8), intent(in), device     :: lambda(5,5)
    real(8), intent(out), device    :: mat(5,5)
    real(8), dimension(5,5), device :: R
    R(:,1) = (/     1.d0,    0.d0,  0.d0,       1.d0,      1.d0/)
    R(:,2) = (/    u - c,    0.d0,  0.d0,          u,     u + c/)
    R(:,3) = (/        v,    0.d0,     c,          v,         v/)
    R(:,4) = (/        w,      -c,  0.d0,          w,         w/)
    R(:,5) = (/H - c * u, - w * c, v * c, 0.5d0 * q2, H + c * u/)
    mat = cumatmul(R, lambda)
    R(:,1) = (/0.5d0 * (b1 + u / c), -0.5d0 * (1.d0 / c + b2 * u), -0.5d0 * b2 * v, -0.5d0 * b2 * w, 0.5d0 * b2/)
    R(:,2) = (/               w / c,                         0.d0,           0.d0 ,       -1.d0 / c,       0.d0/)
    R(:,3) = (/              -v / c,                         0.d0,        1.d0 / c,            0.d0,       0.d0/)
    R(:,4) = (/           1.d0 - b1,                       b2 * u,          b2 * v,          b2 * w,        -b2/)
    R(:,5) = (/0.5d0 * (b1 - u / c),  0.5d0 * (1.d0 / c - b2 * u), -0.5d0 * b2 * v, -0.5d0 * b2 * w, 0.5d0 * b2/)
    mat = cumatmul(mat, R)
  end subroutine calc_A


  attributes(device) subroutine calc_B(u, v, w, c, q2, H, b1, b2, lambda, mat)
    real(8), intent(in)             :: u, v, w, c, q2, H, b1, b2
    real(8), intent(in), device     :: lambda(5,5)
    real(8), intent(out), device    :: mat(5,5)
    real(8), dimension(5,5), device :: R
    R(:,1) = (/     1.d0,       1.d0,   0.d0,  0.d0,      1.d0/)
    R(:,2) = (/        u,          u,     -c,  0.d0,         u/)
    R(:,3) = (/    v - c,          v,   0.d0,  0.d0,     v + c/)
    R(:,4) = (/        w,          w,   0.d0,     c,         w/)
    R(:,5) = (/H - c * v, 0.5d0 * q2, -u * c, w * c, H + c * v/)
    mat = cumatmul(R, lambda)
    R(:,1) = (/0.5d0 * (b1 + v / c), -0.5d0 * b2 * u, -0.5d0 * (1.d0 / c + b2 * v), -0.5d0 * b2 * w, 0.5d0 * b2/)
    R(:,2) = (/           1.d0 - b1,          b2 * u,                       b2 * v,          b2 * w,        -b2/)
    R(:,3) = (/               u / c,       -1.d0 / c,                         0.d0,            0.d0,       0.d0/)
    R(:,4) = (/              -w / c,           0.d0 ,                         0.d0,        1.d0 / c,       0.d0/)
    R(:,5) = (/0.5d0 * (b1 - v / c), -0.5d0 * b2 * u,  0.5d0 * (1.d0 / c - b2 * v), -0.5d0 * b2 * w, 0.5d0 * b2/)
    mat = cumatmul(mat, R)
  end subroutine calc_B


  attributes(device) subroutine calc_C(u, v, w, c, q2, H, b1, b2, lambda, mat)
    real(8), intent(in)             :: u, v, w, c, q2, H, b1, b2
    real(8), intent(in), device     :: lambda(5,5)
    real(8), intent(out), device    :: mat(5,5)
    real(8), dimension(5,5), device :: R
    R(:,1) = (/     1.d0,  0.d0,       1.d0,   0.d0,      1.d0/)
    R(:,2) = (/        u,     c,          u,   0.d0,         u/)
    R(:,3) = (/        v,  0.d0,          v,     -c,         v/)
    R(:,4) = (/    w - c,  0.d0,          w,   0.d0,     w + c/)
    R(:,5) = (/H - c * w, u * c, 0.5d0 * q2, -v * c, H + c * w/)
    mat = cumatmul(R, lambda)
    R(:,1) = (/0.5d0 * (b1 + w / c), -0.5d0 * b2 * u, -0.5d0 * b2 * v, -0.5d0 * (1.d0 / c + b2 * w), 0.5d0 * b2/)
    R(:,2) = (/              -u / c,        1.d0 / c,            0.d0,                         0.d0,       0.d0/)
    R(:,3) = (/           1.d0 - b1,          b2 * u,          b2 * v,                       b2 * w,        -b2/)
    R(:,4) = (/               v / c,            0.d0,       -1.d0 / c,                         0.d0,       0.d0/)
    R(:,5) = (/0.5d0 * (b1 - w / c), -0.5d0 * b2 * u, -0.5d0 * b2 * v,  0.5d0 * (1.d0 / c - b2 * w), 0.5d0 * b2/)
    mat = cumatmul(mat, R)
  end subroutine calc_C


  attributes(device) subroutine calc_mat(u, v, w, c, q2, H, b1, b2, lambda, mat)
    real(8), intent(in)             :: u, v, w, c, q2, H, b1, b2
    real(8), intent(in), device     :: lambda(5,5)
    real(8), intent(out), device    :: mat(5,5)
    real(8), dimension(5,5), device :: R
    R(:,1) = (/     1.d0, 0.d0, 0.d0,       1.d0,      1.d0/)
    R(:,2) = (/    u - c, 0.d0, 0.d0,          u,     u + c/)
    R(:,3) = (/        v, 1.d0, 0.d0,          v,         v/)
    R(:,4) = (/        w, 0.d0, 1.d0,          w,         w/)
    R(:,5) = (/H - c * u,    v,    w, 0.5d0 * q2, H + c * u/)
    mat = cumatmul(R, lambda)
    R(:,1) = (/0.5d0 * (b1 + u / c), -0.5d0 * (b2 * u + 1.d0 / c), -0.5d0 * b2 * v, -0.5d0 * b2 * w, 0.5d0 * b2/)
    R(:,2) = (/                  -v,                         0.d0,            1.d0,            0.d0,       0.d0/)
    R(:,3) = (/                  -w,                         0.d0,            0.d0,            1.d0,       0.d0/)
    R(:,4) = (/          -b1 + 1.d0,                       b2 * u,          b2 * v,          b2 * w,        -b2/)
    R(:,5) = (/0.5d0 * (b1 - u / c), -0.5d0 * (b2 * u - 1.d0 / c), -0.5d0 * b2 * v, -0.5d0 * b2 * w, 0.5d0 * b2/)
    mat = cumatmul(mat, R)
  end subroutine calc_mat

  ! 3D only
  attributes(device) function Roe(id, rho, p, V, Normal) result(F)
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(2), device           :: rho, p
    real(8), intent(in), dimension(2,dimension), device :: V
    real(8), intent(in), dimension(dimension+2), device :: Normal
    real(8), dimension(dimension+2), device :: F, Fl, Fr, dQ
    real(8) el, er, Hl, Hr, q2, rho_ave, H_ave, c, H, b1, b2
    real(8), dimension(dimension), device               :: V_ave
    real(8), dimension(dimension+2,dimension+2), device :: lambda, mat
    el = p(1) / (gamma - 1.d0) + 0.5d0 * rho(1) * (V(1,1)**2 + V(1,2)**2 + V(1,3)**2)
    er = p(2) / (gamma - 1.d0) + 0.5d0 * rho(2) * (V(2,1)**2 + V(2,2)**2 + V(2,3)**2)
    Hl = (el + p(1)) / rho(1)
    Hr = (er + p(2)) / rho(2)

    rho_ave  = sqrt(rho(1) * rho(2))
    V_ave(:) = (sqrt(rho(1)) * V(1,:) + sqrt(rho(2)) * V(2,:)) / (sqrt(rho(1)) + sqrt(rho(2)))
    q2       = V_ave(1)**2 + V_ave(2)**2 + V_ave(3)**2
    H_ave    = (sqrt(rho(1)) * Hl + sqrt(rho(2)) * Hr) / (sqrt(rho(1)) + sqrt(rho(2)))
    c        = sqrt((gamma - 1.d0) * abs(H_ave - 0.5d0 * q2))
    H        = c**2 / (gamma - 1.d0) + 0.5d0 * q2
    b1       = 0.5d0 * q2 * (gamma - 1.d0) / c**2
    b2       = (gamma - 1.d0) / c**2

    lambda(:,1) = (/abs(V_ave(1) - c), 0.d0, 0.d0, 0.d0,  0.d0/)
    lambda(:,2) = (/ 0.d0,    abs(V_ave(1)), 0.d0, 0.d0,  0.d0/)
    lambda(:,3) = (/ 0.d0, 0.d0,    abs(V_ave(1)), 0.d0,  0.d0/)
    lambda(:,4) = (/ 0.d0, 0.d0, 0.d0,    abs(V_ave(1)),  0.d0/)
    lambda(:,5) = (/ 0.d0, 0.d0, 0.d0, 0.d0, abs(V_ave(1) + c)/)
    if (id == 1) then
      call calc_A(V_ave(1), V_ave(2), V_ave(3), c, q2, H, b1, b2, lambda, mat)
    elseif (id == 2) then
      call calc_B(V_ave(1), V_ave(2), V_ave(3), c, q2, H, b1, b2, lambda, mat)
    elseif (id == 3) then
      call calc_C(V_ave(1), V_ave(2), V_ave(3), c, q2, H, b1, b2, lambda, mat)
    endif
    !call calc_mat(V_ave(1), V_ave(2), V_ave(3), c, q2, H, b1, b2, lambda, mat)

    Fl(1) = rho(1) * V(1,id)
    Fr(1) = rho(2) * V(2,id)
    Fl(2:dimension+1) = Fl(1) * V(1,:)
    Fr(2:dimension+1) = Fr(1) * V(2,:)
    Fl(dimension+2) = (el + p(1)) * V(1,id)
    Fr(dimension+2) = (er + p(2)) * V(2,id)
    Fl(:) = Fl(:) + p(1) * Normal(:)
    Fr(:) = Fr(:) + p(2) * Normal(:)

    dQ(1)             = -rho(1)          + rho(2)
    dQ(2:dimension+1) = -rho(1) * V(1,:) + rho(2) * V(2,:)
    dQ(dimension+2)   = -el              + er
    F(:) = 0.5d0 * (Fl(:) + Fr(:) - cumatmul(mat(:,:), dQ(:)))
  end function Roe
end module calc_roe

