module calc_roe
  use mod_globals, only : gamma
  use calc_common
  use calc_physical_quantities
  use calc_mat
  implicit none
contains
  attributes(device) function Roe(id_dim,Qsl,Qsr,Normal) result(F)
    integer, intent(in), value :: id_dim
    real(8), intent(in), dimension(4), device :: Qsl, Qsr
    real(8), intent(in), dimension(4), device :: Normal
    real(8), dimension(4), device :: F, Fl, Fr, Ql, Qr
    real(8) rhol, rhor, pl, pr, el, er, Hl, Hr, rho_ave, H_ave, c_ave
    real(8), dimension(2) :: Vl, Vr, V_ave
    real(8) mat(4,4)
    call set_q(Qsl,Qsr,rhol,rhor,pl,pr,Vl,Vr)
    el = energy(pl,rhol,Vl(1),Vl(2))
    er = energy(pr,rhor,Vr(1),Vr(2))
    Hl = ENTHALPY(el,pl,rhol)
    Hr = ENTHALPY(er,pr,rhor)

    rho_ave = sqrt(rhol * rhor)
    V_ave(1) = (sqrt(rhol) * Vl(1) + sqrt(rhor) * Vr(1)) / (sqrt(rhol) + sqrt(rhor))
    V_ave(2) = (sqrt(rhol) * Vl(2) + sqrt(rhor) * Vr(2)) / (sqrt(rhol) + sqrt(rhor))
    H_ave = (sqrt(rhol) * Hl + sqrt(rhor) * Hr) / (sqrt(rhol) + sqrt(rhor))
    c_ave = sqrt((gamma - 1.d0) * (H_ave - 0.5d0 * (V_ave(1)**2 + V_ave(2)**2)))

    Ql(:) = (/rhol, rhol * Vl(1), rhol * Vl(2), el/)
    Qr(:) = (/rhor, rhor * Vr(1), rhor * Vr(2), er/)

    Fl(1) = rhol * Vl(id_dim)
    Fr(1) = rhor * Vr(id_dim)
    Fl(2:3) = (/Fl(1) * Vl(1), Fl(1) * Vl(2)/)
    Fr(2:3) = (/Fr(1) * Vr(1), Fr(1) * Vr(2)/)
    Fl(4) = (el + pl) * Vl(id_dim)
    Fr(4) = (er + pr) * Vr(id_dim)
    Fl(:) = Fl(:) + pl * Normal(:)
    Fr(:) = Fr(:) + pr * Normal(:)

    mat = calc_AB(id_dim, rho_ave, H_ave, c_ave, V_ave)
    F(:) = 0.5d0 * (Fl(:) + Fr(:) - matmul(mat(:,:), Qr(:) - Ql(:)))
  end function Roe
end module calc_roe

