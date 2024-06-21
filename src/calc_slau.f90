module calc_slau
  use mod_globals, only : id_slau, dimension, dp_max
  use calc_common
  use calc_common_dim
  use calc_physical_quantities
  implicit none
  interface f_slau
    module procedure sd_slau, basic_slau
  end interface
contains
  attributes(device) function basic_slau(id_slau,p,dp,dp_max,V,c,x) result(fslau)
    integer(kind=2), intent(in), value :: id_slau
    real(8), intent(in), value :: p, dp, dp_max, V, c, x
    real(8) fslau
    fslau = x
  end function basic_slau

  attributes(device) function sd_slau(id_slau,p,dp,dp_max,V,c,x) result(fslau)
    integer(kind=4), intent(in), value :: id_slau
    real(8), intent(in), value :: p, dp, dp_max, V, c, x
    real(8) :: Csd1 = 0.1d0
    real(8) :: Csd2 = 10.d0
    real(8) fslau, M, theta
    M = V / c
    theta = min(1.d0, ((Csd2 * abs(dp) / p + Csd1) / (abs(dp_max) / p + Csd1))**2)
    fslau = theta * 0.5d0 * (abs(M + 1.d0) + abs(M - 1.d0) - 2.d0 * abs(M))
  end function sd_slau

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) subroutine calc_quantities_AUSM(pl,pr,rhol,rhor,Vl,Vr,el,er,Hl,Hr,cl,cr)
    real(8), intent(in), value :: pl, pr, rhol, rhor
    real(8), intent(in), dimension(dimension), device :: Vl, Vr
    real(8), intent(out) :: el, er, Hl, Hr, cl, cr
    el = energy(pl,rhol,Vl(:))
    er = energy(pr,rhor,Vr(:))
    Hl = ENTHALPY(el,pl,rhol)
    Hr = ENTHALPY(er,pr,rhor)
    cl = speed_of_sound(pl,rhol)
    cr = speed_of_sound(pr,rhor)
  end subroutine calc_quantities_AUSM

  attributes(device) subroutine calc_beta(Mp,Mm,bp,bm)
    real(8), intent(in), value :: Mp, Mm
    real(8), intent(out) :: bp, bm
    if (abs(Mp) < 1.d0) then
      bp = 0.25d0 * (2.d0 - Mp) * (Mp + 1.d0) ** 2
    else
      bp = 0.5d0 * (1.d0 + sign(1.d0, Mp))
    endif
    if (abs(Mm) < 1.d0) then
      bm = 0.25d0 * (2.d0 + Mm) * (Mm - 1.d0) ** 2
    else
      bm = 0.5d0 * (1.d0 + sign(1.d0, -Mm))
    endif
  end subroutine calc_beta

  attributes(device) function flux_AUSM(mass,pressure,Hl,Hr,fd,Vl,Vr,Normal) result(F)
    real(8), intent(in), value                          :: mass, pressure, Hl, Hr, fd
    real(8), intent(in), dimension(dimension), device   :: Vl, Vr
    real(8), intent(in), dimension(dimension+2), device :: Normal
    real(8), dimension(dimension+2) :: F, phil, phir
    phil(1) = 1.d0
    phil(2:dimension+1) = Vl(:)
    phil(dimension+2) = Hl
    phir(1) = 1.d0
    phir(2:dimension+1) = Vr(:)
    phir(dimension+2) = Hr
    !F(:) = 0.5d0 * ((mass + abs(mass)) * phil(:) + (mass - abs(mass)) * phir(:)) + pressure * Normal(:)
    F(:) = 0.5d0 * mass * (phil(:) + phir(:)) - 0.5d0 * min(1.d0, 5.d0 * fd) * abs(mass) * (phir(:) - phil(:)) + pressure * Normal(:)
 end function flux_AUSM

  attributes(device) function SLAU(id_dim,Ql,Qr,Normal,fd) result(F)
    integer, intent(in), value                          :: id_dim
    real(8), intent(in), dimension(dimension+2), device :: Ql, Qr, Normal
    real(8), intent(in), value                          :: fd
    real(8) rhol, rhor, pl, pr, el, er, Hl, Hr, cl, cr, c
    real(8) Vp, Vm, Vt, Vtp, Vtm, Mp, Mm, M, x, fslau, g, p, dp, mass, bp, bm, Pressure
    real(8), dimension(dimension) :: Vl, Vr, xy
    real(8), dimension(dimension+2) :: F

    call set_q(Ql,Qr,rhol,rhor,pl,pr,Vl,Vr)
    call calc_quantities_AUSM(pl,pr,rhol,rhor,Vl,Vr,el,er,Hl,Hr,cl,cr)

    c = 0.5d0 * (cl + cr)
    ! contravariant velocity
    xy = Normal(2:dimension+1)
    Vp = vecsum(Vl, xy)
    Vm = vecsum(Vr, xy)
    !Vp = Vl(id_dim)
    !Vm = Vr(id_dim)
    Mp = Vp / c
    Mm = Vm / c
    M = min(1.d0, sqrt(0.5d0 * q2(Vl(:), Vr(:))) / c)
    x = (1.d0 - M) ** 2
    g = -max(min(Mp, 0.d0), -1.d0) * min(max(Mm, 0.d0), 1.d0)
    Vt = (rhol * abs(Vp) + rhor * abs(Vm)) / (rhol + rhor)
    Vtp = (1.d0 - g) * Vt + g * abs(Vp)
    Vtm = (1.d0 - g) * Vt + g * abs(Vm)
    dp = pr - pl
    call calc_beta(Mp,Mm,bp,bm)
    p = 0.5d0 * (pl + pr)
    fslau = f_slau(id_slau,p,dp,dp_max,Vt,c,x)
    ! mass flux
    mass = 0.5d0 * (rhol * (Vl(id_dim) + Vtp) + rhor * (Vr(id_dim) - Vtm) - fslau * dp / c)
    ! pressure flux
    Pressure = 0.5d0 * (pl + pr + (bp - bm) * (pl - pr) + (1.d0 - x) * (bp + bm - 1.d0) * (pl + pr))
    F(:) = flux_AUSM(mass,Pressure,Hl,Hr,1.d0,Vl,Vr,Normal)
  end function SLAU

  attributes(device) function simpleSLAU(id_dim,Ql,Qr,Normal,fd) result(F)
    integer, intent(in), value                          :: id_dim
    real(8), intent(in), dimension(dimension+2), device :: Ql, Qr, Normal
    real(8), intent(in), value                          :: fd
    real(8) rhol, rhor, pl, pr, el, er, Hl, Hr, cl, cr
    real(8) Vp, Vm, Vt, mass, bp, bm, Pressure
    real(8), dimension(dimension)   :: Vl, Vr, xy
    real(8), dimension(dimension+2) :: F

    call set_q(Ql,Qr,rhol,rhor,pl,pr,Vl,Vr)
    call calc_quantities_AUSM(pl,pr,rhol,rhor,Vl,Vr,el,er,Hl,Hr,cl,cr)

    Vp = Vl(id_dim)
    Vm = Vr(id_dim)
    Vt = (rhol * abs(Vp) + rhor * abs(Vm)) / (rhol + rhor)
    ! mass flux
    mass = 0.5d0 * (rhol * Vl(id_dim) + rhor * Vr(id_dim) - Vt * (rhor - rhol))
    ! pressure flux
    bp = 0.5d0 * (1.d0 + sign(1.d0,  Vp)) 
    bm = 0.5d0 * (1.d0 + sign(1.d0, -Vm))
    Pressure = bp * pl + bm * pr
    F(:) = flux_AUSM(mass,Pressure,Hl,Hr,fd,Vl,Vr,Normal)
  end function simpleSLAU
end module calc_slau

