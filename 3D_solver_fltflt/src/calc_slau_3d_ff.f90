  ! ================================================================
  ! fltflt (FP32-pair) counterparts of SLAU_common/phi/HRSLAU2, bound
  ! into the same `SLAU` generic interface (declared in the enclosing
  ! module) alongside the untouched real(8) originals above -- Fortran
  ! resolves by actual argument type, so calc_slau_kernel.f90.fypp and
  ! the Hybrid kernel files (still real(8)) keep resolving to those.
  ! ================================================================

  attributes(device) subroutine SLAU_common_ff(rho1, rho2, over_rho1, over_rho2, u1, u2, v1, v2, w1, w2, &
                             un1, un2, p1, p2, c, over_c, Mp, Mm, bp, bm, dp, Vtp, Vtm)
    use fltflt
    type(fltflt), intent(in)  :: rho1, rho2, over_rho1, over_rho2, u1, u2, v1, v2, w1, w2, un1, un2, p1, p2
    type(fltflt), intent(out) :: c, over_c, Mp, Mm, bp, bm, dp, Vtp, Vtm
    block
      type(fltflt) cl, cr
      cl = fltflt_sqrt(gamma * p1 * over_rho1)
      cr = fltflt_sqrt(gamma * p2 * over_rho2)
      c  = 0.5d0 * (cl + cr)
    end block
    over_c = 1.d0 / c
    Mp  = un1 * over_c
    Mm  = un2 * over_c
    block
      type(fltflt) g, one_g_Vt
      g   = -fltflt_max(fltflt_min(Mp, 0.0), -1.0) * fltflt_min(fltflt_max(Mm, 0.0), 1.0)
      one_g_Vt = (1.d0 - g) * (rho1 * fltflt_abs(un1) + rho2 * fltflt_abs(un2)) / (rho1 + rho2)
      Vtp = one_g_Vt + g * fltflt_abs(un1)
      Vtm = one_g_Vt + g * fltflt_abs(un2)
    end block
    if (fltflt_abs(Mp) < 1.0) then
      bp = 0.25d0 * (2.d0 - Mp) * fltflt_square(Mp + 1.d0)
    else
      bp = 0.5d0 * (1.d0 + fltflt_sign(fltflt_init(1.0_4), Mp))
    endif
    if (fltflt_abs(Mm) < 1.0) then
      bm = 0.25d0 * (2.d0 + Mm) * fltflt_square(Mm - 1.d0)
    else
      bm = 0.5d0 * (1.d0 + fltflt_sign(fltflt_init(1.0_4), -Mm))
    endif
    dp = -p1 + p2
  end subroutine SLAU_common_ff


  attributes(device) function phi_ff(rho, k, p, over_rho) result(ans)
    use mod_globals, only : gamma
    use mod_constant, only : over_gamma_1_ff
    use fltflt
    type(fltflt), intent(in) :: rho, k, p, over_rho
    type(fltflt) :: ans
    type(fltflt) :: gamma_over_gamma_1
    gamma_over_gamma_1 = gamma * over_gamma_1_ff
    ans = (p * gamma_over_gamma_1 + rho * k) * over_rho
  end function phi_ff


  attributes(device) subroutine HRSLAU2_ff(id_slau, rho1, rho2, u1, u2, v1, v2, w1, w2, &
                         un1, un2, p1, p2, Norm, HR, F1, F2, F3, F4, F5)
    use fltflt
    integer(4), intent(in), value :: id_slau
    type(fltflt), intent(in)      :: rho1, rho2, u1, u2, v1, v2, w1, w2, un1, un2, p1, p2
    type(fltflt), intent(in)      :: Norm(5)
    real(sp), intent(in), value   :: HR
    real(8), intent(out)          :: F1, F2, F3, F4, F5
    type(fltflt) :: c, over_c, Mp, Mm, Vtp, Vtm, dp, bp, bm
    type(fltflt) :: mass, mass1, mass2, Vec2, over_rho1, over_rho2, k1, k2
    type(fltflt) :: pres, phi1, phi2, ff_F1, ff_F2, ff_F3, ff_F4, ff_F5
    over_rho1 = 1.d0 / rho1
    over_rho2 = 1.d0 / rho2
    call SLAU_common_ff(rho1, rho2, over_rho1, over_rho2, u1, u2, v1, v2, w1, w2, un1, un2, &
                        p1, p2, c, over_c, Mp, Mm, bp, bm, dp, Vtp, Vtm)
    k1   = 0.5d0 * (u1*u1 + v1*v1 + w1*w1)
    k2   = 0.5d0 * (u2*u2 + v2*v2 + w2*w2)
    Vec2 = fltflt_sqrt(k1 + k2)
    block
      type(fltflt) M, x
      M  = fltflt_min(Vec2 * over_c, 1.0)
      x  = fltflt_square(1.d0 - M)
      mass  = 0.25d0 * (rho1 * (un1 + Vtp) + rho2 * (un2 - Vtm) - x * dp * over_c)
      mass1 = mass + fltflt_abs(mass)
      mass2 = mass - fltflt_abs(mass)
    end block
    pres = 0.5d0 * (p1 + p2 + (bp - bm) * (-dp) + HR * Vec2 * (bp + bm - 1.d0) * 0.5d0 * (rho1 + rho2) * c)
    phi1 = phi_ff(rho1, k1, p1, over_rho1)
    phi2 = phi_ff(rho2, k2, p2, over_rho2)
    ff_F1 = mass1 + mass2
    ff_F2 = fltflt_fma(mass1, u1,   mass2 * u2)   + pres * Norm(2)
    ff_F3 = fltflt_fma(mass1, v1,   mass2 * v2)   + pres * Norm(3)
    ff_F4 = fltflt_fma(mass1, w1,   mass2 * w2)   + pres * Norm(4)
    ff_F5 = fltflt_fma(mass1, phi1, mass2 * phi2)
    F1 = real(ff_F1%hi, 8) + real(ff_F1%lo, 8)
    F2 = real(ff_F2%hi, 8) + real(ff_F2%lo, 8)
    F3 = real(ff_F3%hi, 8) + real(ff_F3%lo, 8)
    F4 = real(ff_F4%hi, 8) + real(ff_F4%lo, 8)
    F5 = real(ff_F5%hi, 8) + real(ff_F5%lo, 8)
  end subroutine HRSLAU2_ff


  attributes(device) subroutine SLAU1_ff(id_slau, rho1, rho2, u1, u2, v1, v2, w1, w2, &
                                         un1, un2, p1, p2, Norm, HR, F1, F2, F3, F4, F5)
    use fltflt
    integer(2), intent(in), value :: id_slau
    type(fltflt), intent(in)      :: rho1, rho2, u1, u2, v1, v2, w1, w2, un1, un2, p1, p2
    type(fltflt), intent(in)      :: Norm(5)
    real(sp), intent(in), value   :: HR
    real(8), intent(out)          :: F1, F2, F3, F4, F5
    type(fltflt) :: c, over_c, Mp, Mm, Vtp, Vtm, dp, bp, bm
    type(fltflt) :: mass, mass1, mass2, over_rho1, over_rho2, k1, k2
    type(fltflt) :: M, x, pres, phi1, phi2, ff_F1, ff_F2, ff_F3, ff_F4, ff_F5
    over_rho1 = 1.d0 / rho1
    over_rho2 = 1.d0 / rho2
    call SLAU_common_ff(rho1, rho2, over_rho1, over_rho2, u1, u2, v1, v2, w1, w2, un1, un2, &
                        p1, p2, c, over_c, Mp, Mm, bp, bm, dp, Vtp, Vtm)
    k1 = 0.5d0 * (u1*u1 + v1*v1 + w1*w1)
    k2 = 0.5d0 * (u2*u2 + v2*v2 + w2*w2)
    M  = fltflt_min(fltflt_sqrt(k1 + k2) * over_c, 1.0)
    x  = fltflt_square(1.d0 - M)
    mass  = 0.25d0 * (rho1 * (un1 + Vtp) + rho2 * (un2 - Vtm) - x * dp * over_c)
    mass1 = mass + fltflt_abs(mass)
    mass2 = mass - fltflt_abs(mass)
    pres = 0.5d0 * (p1 + p2 + (bp - bm) * (-dp) + (1.d0 - x) * (bp + bm - 1.d0) * (p1 + p2))
    phi1 = phi_ff(rho1, k1, p1, over_rho1)
    phi2 = phi_ff(rho2, k2, p2, over_rho2)
    ff_F1 = mass1 + mass2
    ff_F2 = fltflt_fma(mass1, u1,   mass2 * u2)   + pres * Norm(2)
    ff_F3 = fltflt_fma(mass1, v1,   mass2 * v2)   + pres * Norm(3)
    ff_F4 = fltflt_fma(mass1, w1,   mass2 * w2)   + pres * Norm(4)
    ff_F5 = fltflt_fma(mass1, phi1, mass2 * phi2)
    F1 = real(ff_F1%hi, 8) + real(ff_F1%lo, 8)
    F2 = real(ff_F2%hi, 8) + real(ff_F2%lo, 8)
    F3 = real(ff_F3%hi, 8) + real(ff_F3%lo, 8)
    F4 = real(ff_F4%hi, 8) + real(ff_F4%lo, 8)
    F5 = real(ff_F5%hi, 8) + real(ff_F5%lo, 8)
  end subroutine SLAU1_ff

