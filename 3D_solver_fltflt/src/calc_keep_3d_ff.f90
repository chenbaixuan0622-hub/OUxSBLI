  ! ================================================================
  ! fltflt (FP32-pair) counterpart of KEEP6, bound into the same `KEEP`
  ! generic interface (declared in the enclosing module) alongside the
  ! untouched real(8) original above -- resolved by argument type, so
  ! calc_keep_kernel.f90.fypp and the Hybrid kernel files (still
  ! real(8)) keep resolving to the original KEEP6.
  ! ================================================================

  attributes(device) function KEEP6_ff(id_accuracy, rho, u, v, w, uu, p, T, Normal) result(F)
    integer(8), intent(in), value     :: id_accuracy
    type(fltflt), intent(in), dimension(6), contiguous :: rho, u, v, w, uu, p, T
    type(fltflt), intent(in), dimension(5) :: Normal
    real(8) F(5)
    block
      type(fltflt) :: rho1,rho2,rho3,rho4,rho5,rho6
      type(fltflt) :: u1,u2,u3,u4,u5,u6, v1,v2,v3,v4,v5,v6, w1,w2,w3,w4,w5,w6
      type(fltflt) :: uu1,uu2,uu3,uu4,uu5,uu6, p1,p2,p3,p4,p5,p6, T1,T2,T3,T4,T5,T6
      type(fltflt) :: ff_RV1, ff_RV2, ff_RV3, ff_RV4, ff_RV5, ff_RV6
      type(fltflt) :: ff_RV23, ff_RV456
      type(fltflt) :: ff_RV3_RV5, ff_RV1_RV2_RV4, ff_RV1_RV3_RV6, ff_RV2_RV5
      type(fltflt) :: ff_F1, ff_ruu, ff_ruv, ff_ruw, ff_ene, ff_F5
      type(fltflt) :: pres, p34
      type(fltflt) :: d1, d2, d3, d4, d5, d6
      rho1=rho(1); rho2=rho(2); rho3=rho(3); rho4=rho(4); rho5=rho(5); rho6=rho(6)
      u1=u(1); u2=u(2); u3=u(3); u4=u(4); u5=u(5); u6=u(6)
      v1=v(1); v2=v(2); v3=v(3); v4=v(4); v5=v(5); v6=v(6)
      w1=w(1); w2=w(2); w3=w(3); w4=w(4); w5=w(5); w6=w(6)
      uu1=uu(1); uu2=uu(2); uu3=uu(3); uu4=uu(4); uu5=uu(5); uu6=uu(6)
      p1=p(1); p2=p(2); p3=p(3); p4=p(4); p5=p(5); p6=p(6)
      T1=T(1); T2=T(2); T3=T(3); T4=T(4); T5=T(5); T6=T(6)
      ff_RV1 = (rho3 + rho4) * (uu3 + uu4)
      ff_RV2 = (rho3 + rho5) * (uu3 + uu5)
      ff_RV3 = (rho2 + rho4) * (uu2 + uu4)
      ff_RV4 = (rho3 + rho6) * (uu3 + uu6)
      ff_RV5 = (rho2 + rho5) * (uu2 + uu5)
      ff_RV6 = (rho1 + rho4) * (uu1 + uu4)
      ! one_120 is a real(8) parameter (not fltflt) in the Hybrid-kernel copies
      ! of this file, so it can't be a fltflt_dot* operand there -- fold it in
      ! via a plain multiply instead of fusing it into the dot2 below.
      ff_RV23  = ff_RV2 + ff_RV3
      ff_RV456 = fltflt_add3(ff_RV4, ff_RV5, ff_RV6)
      ff_F1 = fltflt_dot2(c0_375_ff, ff_RV1, -c0_075_ff, ff_RV23) + ff_RV456 * one_120
      ff_RV1 = c0_1875_ff * ff_RV1
      ff_RV2 = c0_0375_ff * ff_RV2
      ff_RV3 = c0_0375_ff * ff_RV3
      ff_RV4 = ff_RV4 * one_240
      ff_RV5 = ff_RV5 * one_240
      ff_RV6 = ff_RV6 * one_240
      ff_RV3_RV5     = -ff_RV3 + ff_RV5
      ff_RV1_RV2_RV4 = fltflt_add3(ff_RV1, -ff_RV2, ff_RV4)
      ff_RV1_RV3_RV6 = fltflt_add3(ff_RV1, -ff_RV3, ff_RV6)
      ff_RV2_RV5     = -ff_RV2 + ff_RV5
      ! pres: same one_60 caveat as one_120 above
      p34  = p3 + p4
      pres = fltflt_add3(p1, fltflt_dot3(-c8_ff, p2, c37_ff, p34, -c8_ff, p5), p6) * one_60
      ff_ruu = fltflt_fma(pres, Normal(2), &
                 fltflt_dot6(ff_RV6, u1, ff_RV3_RV5, u2, ff_RV1_RV2_RV4, u3, &
                             ff_RV1_RV3_RV6, u4, ff_RV2_RV5, u5, ff_RV4, u6))
      ff_ruv = fltflt_fma(pres, Normal(3), &
                 fltflt_dot6(ff_RV6, v1, ff_RV3_RV5, v2, ff_RV1_RV2_RV4, v3, &
                             ff_RV1_RV3_RV6, v4, ff_RV2_RV5, v5, ff_RV4, v6))
      ff_ruw = fltflt_fma(pres, Normal(4), &
                 fltflt_dot6(ff_RV6, w1, ff_RV3_RV5, w2, ff_RV1_RV2_RV4, w3, &
                             ff_RV1_RV3_RV6, w4, ff_RV2_RV5, w5, ff_RV4, w6))
      ! internal + kinetic energy
      ff_ene = fltflt_dot6(ff_RV6, T1, ff_RV3_RV5, T2, ff_RV1_RV2_RV4, T3, &
                           ff_RV1_RV3_RV6, T4, ff_RV2_RV5, T5, ff_RV4, T6) * R_over_gamma_1_ff
      d1 = fltflt_dot3(u3,u4, v3,v4, w3,w4)
      d2 = fltflt_dot3(u3,u5, v3,v5, w3,w5)
      d3 = fltflt_dot3(u2,u4, v2,v4, w2,w4)
      d4 = fltflt_dot3(u3,u6, v3,v6, w3,w6)
      d5 = fltflt_dot3(u2,u5, v2,v5, w2,w5)
      d6 = fltflt_dot3(u1,u4, v1,v4, w1,w4)
      ff_ene = ff_ene + fltflt_dot6(ff_RV1,d1, -ff_RV2,d2, -ff_RV3,d3, ff_RV4,d4, ff_RV5,d5, ff_RV6,d6)
      ! pressure diffusion
      ff_F5 = fltflt_add3(ff_ene, &
                fltflt_dot2(c0_75_ff, fltflt_dot2(uu3,p4, uu4,p3), &
                            -c0_15_ff, fltflt_dot4(uu3,p5, uu5,p3, uu2,p4, uu4,p2)), &
                fltflt_dot6(uu3,p6, uu6,p3, uu2,p5, uu5,p2, uu1,p4, uu4,p1) * one_60)
      F(1) = real(ff_F1%hi,8) + real(ff_F1%lo,8)
      F(2) = real(ff_ruu%hi,8) + real(ff_ruu%lo,8)
      F(3) = real(ff_ruv%hi,8) + real(ff_ruv%lo,8)
      F(4) = real(ff_ruw%hi,8) + real(ff_ruw%lo,8)
      F(5) = real(ff_F5%hi,8) + real(ff_F5%lo,8)
    end block
  end function KEEP6_ff


  attributes(device) function KEEP2_ff(id_accuracy, rho, u, v, w, uu, p, T, Normal) result(F)
    integer(2), intent(in), value     :: id_accuracy
    type(fltflt), intent(in), dimension(2), contiguous :: rho, u, v, w, uu, p, T
    type(fltflt), intent(in), dimension(5) :: Normal
    real(8) F(5)
    block
      type(fltflt) :: rho1, rho2, u1, u2, v1, v2, w1, w2, uu1, uu2, p1, p2, T1, T2
      type(fltflt) :: ff_F1, ff_F2, ff_F3, ff_F4, ff_F5
      type(fltflt) :: u12, v12, w12, p12, coef
      rho1=rho(1); rho2=rho(2)
      u1=u(1); u2=u(2)
      v1=v(1); v2=v(2)
      w1=w(1); w2=w(2)
      uu1=uu(1); uu2=uu(2)
      p1=p(1); p2=p(2)
      T1=T(1); T2=T(2)
      ff_F1 = c0_25_ff * (rho1 + rho2) * (uu1 + uu2)
      u12 = u1 + u2;  v12 = v1 + v2;  w12 = w1 + w2;  p12 = p1 + p2
      ff_F2 = c0_5_ff * fltflt_dot2(ff_F1, u12, p12, Normal(2))
      ff_F3 = c0_5_ff * fltflt_dot2(ff_F1, v12, p12, Normal(3))
      ff_F4 = c0_5_ff * fltflt_dot2(ff_F1, w12, p12, Normal(4))
      ff_F5 = ff_F1 * c0_5_ff * (T1 + T2) * R_over_gamma_1_ff ! internal energy
      ff_F5 = fltflt_fma(c0_5_ff, fltflt_dot2(uu1, p2, uu2, p1), ff_F5) ! pressure diffusion
      coef = c0_5_ff * ff_F1
      ff_F5 = fltflt_fma(coef, fltflt_dot3(u1, u2, v1, v2, w1, w2), ff_F5) ! kinetic energy
      F(1) = real(ff_F1%hi,8) + real(ff_F1%lo,8)
      F(2) = real(ff_F2%hi,8) + real(ff_F2%lo,8)
      F(3) = real(ff_F3%hi,8) + real(ff_F3%lo,8)
      F(4) = real(ff_F4%hi,8) + real(ff_F4%lo,8)
      F(5) = real(ff_F5%hi,8) + real(ff_F5%lo,8)
    end block
  end function KEEP2_ff


  attributes(device) function KEEP4_ff(id_accuracy, rho, u, v, w, uu, p, T, Normal) result(F)
    integer(4), intent(in), value     :: id_accuracy
    type(fltflt), intent(in), dimension(4), contiguous :: rho, u, v, w, uu, p, T
    type(fltflt), intent(in), dimension(5) :: Normal
    real(8) F(5)
    block
      type(fltflt) :: rho1,rho2,rho3,rho4
      type(fltflt) :: u1,u2,u3,u4, v1,v2,v3,v4, w1,w2,w3,w4
      type(fltflt) :: uu1,uu2,uu3,uu4, p1,p2,p3,p4, T1,T2,T3,T4
      type(fltflt) :: ff_RV1, ff_RV2, ff_RV3, ff_RV1_RV2, ff_RV1_RV3, ff_RV23
      type(fltflt) :: ff_F1, ff_F2, ff_F3, ff_F4, ff_ene, ff_F5
      type(fltflt) :: pres, p23, e1, e2, e3
      rho1=rho(1); rho2=rho(2); rho3=rho(3); rho4=rho(4)
      u1=u(1); u2=u(2); u3=u(3); u4=u(4)
      v1=v(1); v2=v(2); v3=v(3); v4=v(4)
      w1=w(1); w2=w(2); w3=w(3); w4=w(4)
      uu1=uu(1); uu2=uu(2); uu3=uu(3); uu4=uu(4)
      p1=p(1); p2=p(2); p3=p(3); p4=p(4)
      T1=T(1); T2=T(2); T3=T(3); T4=T(4)
      ff_RV1 = (rho2 + rho3) * (uu2 + uu3)
      ff_RV2 = (rho2 + rho4) * (uu2 + uu4)
      ff_RV3 = (rho1 + rho3) * (uu1 + uu3)
      ! one_24 is a real(8) parameter (not fltflt) in the Hybrid-kernel copies
      ! of this file, so it can't be a fltflt_dot* operand there.
      ff_RV23 = ff_RV2 + ff_RV3
      ff_F1 = one_third_ff * ff_RV1 - ff_RV23 * one_24
      ff_RV1 = one_sixth_ff * ff_RV1
      ff_RV2 = one_48 * ff_RV2
      ff_RV3 = one_48 * ff_RV3
      ff_RV1_RV2 = ff_RV1 - ff_RV2
      ff_RV1_RV3 = ff_RV1 - ff_RV3
      ! pres: same one_24-style caveat applies here to seven_twelfth
      p23  = p2 + p3
      pres = fltflt_dot2(-one_twelfth_ff, p1, -one_twelfth_ff, p4) + seven_twelfth * p23
      ff_F2 = fltflt_fma(pres, Normal(2), &
                fltflt_dot4(-ff_RV3, u1, ff_RV1_RV2, u2, ff_RV1_RV3, u3, -ff_RV2, u4))
      ff_F3 = fltflt_fma(pres, Normal(3), &
                fltflt_dot4(-ff_RV3, v1, ff_RV1_RV2, v2, ff_RV1_RV3, v3, -ff_RV2, v4))
      ff_F4 = fltflt_fma(pres, Normal(4), &
                fltflt_dot4(-ff_RV3, w1, ff_RV1_RV2, w2, ff_RV1_RV3, w3, -ff_RV2, w4))
      ! internal energy
      ff_ene = fltflt_dot4(-ff_RV3, T1, ff_RV1_RV2, T2, ff_RV1_RV3, T3, -ff_RV2, T4) * R_over_gamma_1_ff
      ! kinetic energy
      e1 = fltflt_dot3(u2,u3, v2,v3, w2,w3)
      e2 = fltflt_dot3(u2,u4, v2,v4, w2,w4)
      e3 = fltflt_dot3(u1,u3, v1,v3, w1,w3)
      ff_ene = ff_ene + fltflt_dot3(ff_RV1, e1, -ff_RV2, e2, -ff_RV3, e3)
      ! pressure diffusion
      ff_F5 = ff_ene + fltflt_dot2(two_third_ff, fltflt_dot2(uu2,p3, uu3,p2), &
                                   -one_twelfth_ff, fltflt_dot4(uu2,p4, uu4,p2, uu1,p3, uu3,p1))
      F(1) = real(ff_F1%hi,8) + real(ff_F1%lo,8)
      F(2) = real(ff_F2%hi,8) + real(ff_F2%lo,8)
      F(3) = real(ff_F3%hi,8) + real(ff_F3%lo,8)
      F(4) = real(ff_F4%hi,8) + real(ff_F4%lo,8)
      F(5) = real(ff_F5%hi,8) + real(ff_F5%lo,8)
    end block
  end function KEEP4_ff

