  pure attributes(device) function KEEP2(id_accuracy, rho, u, v,  uu, p, T, Normal) result(F)
    integer(2), intent(in), value     :: id_accuracy
    real(8), intent(in), dimension(2) :: rho, u, v, uu, p, T
    real(8), intent(in), dimension(5) :: Normal
    real(8) F(4)
    F(1) = 0.25d0 * (rho(1) + rho(2)) * (uu(1) + uu(2))
    F(2) = 0.5d0 * (F(1) * (u(1) + u(2)) + (p(1) + p(2)) * Normal(2))
    F(3) = 0.5d0 * (F(1) * (v(1) + v(2)) + (p(1) + p(2)) * Normal(3))
    F(4) = F(1) * 0.5d0 * (T(1) + T(2)) * R_over_gamma_1 ! internal energy
    F(4) = F(4) + 0.5d0 * (uu(1) * p(2) + uu(2) * p(1)) ! pressure diffusion
    F(4) = F(4) + 0.5d0 * F(1) * (u(1) * u(2) + v(1) * v(2) ) ! kinetic energy
  end function KEEP2


  pure attributes(device) function KEEP4(id_accuracy, rho, u, v, uu, p, T, Normal) result(F)
    integer(4), intent(in), value     :: id_accuracy
    real(8), intent(in), dimension(4) :: rho, u, v, uu, p, T
    real(8), intent(in), dimension(5) :: Normal
    real(8) F(4), RV1, RV2, RV3, RV1_RV2, RV1_RV3 
    real(8) u1, u2, u3, u4, v1, v2, v3, v4
    RV1 = (rho(2) + rho(3)) * (uu(2) + uu(3))
    RV2 = (rho(2) + rho(4)) * (uu(2) + uu(4))
    RV3 = (rho(1) + rho(3)) * (uu(1) + uu(3))
    F(1) = one_third * RV1 - (RV2 + RV3) * one_24
    RV1 = one_sixth * RV1
    RV2 = one_48 * RV2
    RV3 = one_48 * RV3
    RV1_RV2 = RV1 - RV2
    RV1_RV3 = RV1 - RV3
    u1 = u(1); u2 = u(2); u3 = u(3); u4 = u(4)
    v1 = v(1); v2 = v(2); v3 = v(3); v4 = v(4)
    
    block
      real(8) pres
      pres = -one_twelfth * p(1) + seven_twelfth * (p(2) + p(3)) - one_twelfth * p(4)
      F(2) = -RV3 * u1 + RV1_RV2 * u2 + RV1_RV3 * u3 - RV2 * u4 + pres * Normal(2)
      F(3) = -RV3 * v1 + RV1_RV2 * v2 + RV1_RV3 * v3 - RV2 * v4 + pres * Normal(3)
      
    end block
    block
      real(8) ene
      ! internal energy
      ene  = (-RV3 * T(1) + RV1_RV2 * T(2) + RV1_RV3 * T(3) - RV2 * T(4)) * R_over_gamma_1
      ! kinetic energy
      ene  = ene + ((RV1 * (u2*u3 + v2*v3 )) &
                  - (RV2 * (u2*u4 + v2*v4 ) + &
                    RV3 * (u1*u3 + v1*v3 )))
      ! pressure diffusion
      F(4) = ene + (two_third * (uu(2)*p(3)+uu(3)*p(2)) &
                - one_twelfth * (uu(2)*p(4)+uu(4)*p(2) &
                                +uu(1)*p(3)+uu(3)*p(1)))
    end block
  end function KEEP4

 
 pure attributes(device) function KEEP6(id_accuracy, rho, u, v, uu, p, T, Normal) result(F)
    integer(8), intent(in), value     :: id_accuracy
    real(8), intent(in), dimension(6) :: rho, u, v, uu, p, T
    real(8), intent(in), dimension(5) :: Normal
    real(8) F(4), RV1, RV2, RV3, RV4, RV5, RV6, RV3_RV5, RV2_RV5, RV1_RV2_RV4, RV1_RV3_RV6
    RV1 = (rho(3) + rho(4)) * (uu(3) + uu(4))
    RV2 = (rho(3) + rho(5)) * (uu(3) + uu(5))
    RV3 = (rho(2) + rho(4)) * (uu(2) + uu(4))
    RV4 = (rho(3) + rho(6)) * (uu(3) + uu(6))
    RV5 = (rho(2) + rho(5)) * (uu(2) + uu(5))
    RV6 = (rho(1) + rho(4)) * (uu(1) + uu(4))
    F(1) = 0.375d0 * RV1 - 0.075d0 * (RV2 + RV3) + (RV4 + RV5 + RV6) * one_120
    RV1 = 0.1875d0 * RV1
    RV2 = 0.0375d0 * RV2
    RV3 = 0.0375d0 * RV3
    RV4 = RV4 * one_240
    RV5 = RV5 * one_240
    RV6 = RV6 * one_240
    RV3_RV5 = -RV3 + RV5
    RV1_RV2_RV4 = RV1 - RV2 + RV4
    RV1_RV3_RV6 = RV1 - RV3 + RV6
    RV2_RV5 = -RV2 + RV5
    block
      real(8) pres, ruu, ruv, ruw
      pres = (p(1) - 8.d0 * p(2) + 37.d0 * (p(3) + p(4)) - 8.d0 * p(5) + p(6)) * one_60
      ruu = RV6*u(1) + RV3_RV5*u(2) + RV1_RV2_RV4*u(3) + RV1_RV3_RV6*u(4) + RV2_RV5*u(5) + RV4*u(6)
      ruu = ruu + pres * Normal(2)
      ruv = RV6*v(1) + RV3_RV5*v(2) + RV1_RV2_RV4*v(3) + RV1_RV3_RV6*v(4) + RV2_RV5*v(5) + RV4*v(6)
      ruv = ruv + pres * Normal(3)
      F(2) = ruu
      F(3) = ruv
    
    end block
    block
      real(8) ene
      ! internal energy
      ene = (RV6 * T(1) + RV3_RV5 * T(2) + RV1_RV2_RV4 * T(3) + RV1_RV3_RV6 * T(4) + RV2_RV5 * T(5) + RV4 * T(6)) * R_over_gamma_1
      ! kinetic energy
      ene = ene + RV1 * (u(3)*u(4) + v(3)*v(4)) - &
                 (RV2 * (u(3)*u(5) + v(3)*v(5)) + &
                  RV3 * (u(2)*u(4) + v(2)*v(4))) + &
                 (RV4 * (u(3)*u(6) + v(3)*v(6)) + &
                  RV5 * (u(2)*u(5) + v(2)*v(5)) + &
                  RV6 * (u(1)*u(4) + v(1)*v(4)))
      ! pressure diffusion
      F(4) = ene + (0.75d0 * (uu(3)*p(4) + uu(4)*p(3)) &
                  - 0.15d0 * (uu(3)*p(5) + uu(5)*p(3) + &
                              uu(2)*p(4) + uu(4)*p(2)) + &
                             (uu(3)*p(6) + uu(6)*p(3) + &
                              uu(2)*p(5) + uu(5)*p(2) + &
                              uu(1)*p(4) + uu(4)*p(1)) * one_60)
   end block
  end function KEEP6

