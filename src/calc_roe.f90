module calc_roe
  use cutensorex
  use mod_globals, only : dimension, gamma
  use calc_common
  use calc_common_dim
  use calc_physical_quantities
  use calc_mat
  implicit none
contains
  attributes(device) function Roe(id_dim,Qsl,Qsr,Normal) result(F)
    integer, intent(in), value :: id_dim
    real(8), intent(in), dimension(dimension+2), device :: Qsl, Qsr
    real(8), intent(in), dimension(dimension+2), device :: Normal
    real(8), dimension(dimension+2), device :: F, Fl, Fr, Ql, Qr, dQ
    real(8) rhol, rhor, pl, pr, el, er, Hl, Hr, rho_ave, H_ave, c_ave
    real(8), dimension(dimension) :: Vl, Vr, V_ave
    real(8) mat(dimension+2,dimension+2)
    call set_q(Qsl,Qsr,rhol,rhor,pl,pr,Vl,Vr)
    el = energy(pl,rhol,Vl(:))
    er = energy(pr,rhor,Vr(:))
    Hl = ENTHALPY(el,pl,rhol)
    Hr = ENTHALPY(er,pr,rhor)

    rho_ave = sqrt(rhol * rhor)
    V_ave(:) = (sqrt(rhol) * Vl(:) + sqrt(rhor) * Vr(:)) / (sqrt(rhol) + sqrt(rhor))
    H_ave = (sqrt(rhol) * Hl + sqrt(rhor) * Hr) / (sqrt(rhol) + sqrt(rhor))
    c_ave = sqrt((gamma - 1.d0) * (H_ave - 0.5d0 * sum(V_ave(:)**2)))

    Ql(1) = rhol
    Ql(2:dimension+1) = rhol * Vl(:)
    Ql(dimension+2) = el
    Qr(1) = rhor
    Qr(2:dimension+1) = rhor * Vr(:)
    Qr(dimension+2) = er

    Fl(1) = rhol * Vl(id_dim)
    Fr(1) = rhor * Vr(id_dim)
    Fl(2:dimension+1) = Fl(1) * Vl(:)
    Fr(2:dimension+1) = Fr(1) * Vr(:)
    Fl(dimension+2) = (el + pl) * Vl(id_dim)
    Fr(dimension+2) = (er + pr) * Vr(id_dim)
    Fl(:) = Fl(:) + pl * Normal(:)
    Fr(:) = Fr(:) + pr * Normal(:)

    mat(:,:) = calc_AB(id_dim, rho_ave, H_ave, c_ave, V_ave)
    dQ(:) = Qr(:) - Ql(:)
    F(:) = 0.5d0 * (Fl(:) + Fr(:) - cumatmul(mat(:,:), dQ(:)))
  end function Roe
end module calc_roe

