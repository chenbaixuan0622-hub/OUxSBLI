module calc_keep
  use mod_globals, only : dimension, gamma
  use calc_common_dim
  use calc_term
  use calc_mat
  implicit none
contains
  ! KEEP PE scheme
  attributes(device) function KEEP2(id,rho,p,V,Normal) result(F)
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(2), device           :: rho, p
    real(8), intent(in), dimension(2,dimension), device :: V
    real(8), intent(in), dimension(dimension), device   :: Normal
    real(8), dimension(dimension+2) :: F
    real(8)                         :: Rho_m, P_m, PRho
    real(8), dimension(dimension)   :: V_m, V1, V2
    V1(:) = V(1,:)
    V2(:) = V(2,:)
    Rho_m = 0.5d0 * (rho(1) + rho(2))
    ! contravariant velocity
    V_m(:) = 0.5d0 * (V1(:) + V2(:)) * Normal(id)
    P_m = 0.5d0 * (p(1) + p(2))
    PRho = 0.5d0 * (p(1) / rho(1) + p(2) / rho(2))
    F(1) = Rho_m * V_m(id)
    F(2:dimension+1) = F(1) * V_m(:) + P_m * Normal(:)
    F(dimension+2) = V_m(id) * P_m / (gamma - 1.d0) &
    & + 0.5d0 * F(1) * vecsum(V1, V2) &
    & + 0.5d0 * (V(1,id) * p(2) + V(2,id) * p(1)) * Normal(id)
  end function KEEP2

  ! KEEP PE scheme
  attributes(device) function KEEP4(id,rho,p,V,Normal) result(F)
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(4), device           :: rho, p
    real(8), intent(in), dimension(4,dimension), device :: V
    real(8), intent(in), dimension(dimension), device   :: Normal
    real(8), dimension(4)           :: P_over_Rho
    real(8), dimension(dimension+2) :: F
    real(8), dimension(3)           :: V_m, RhoV, RhoVIE, RhoVKE, VP, Energy
    real(8), dimension(3,dimension) :: RhoVV_P
    integer i
    RhoV(:) = RhoPhi(rho(:), V(:,id))
    ! energy equation
    !P_over_Rho(:) = p(:) / rho(:)
    V_m(:) = Phi(V(:,id))
    RhoVIE(:) = RhoPhiU(V_m(:), p(:)) / (gamma - 1.d0)
    RhoVKE(:) = RhoUPhiPhi(RhoV(:), V(:,:))
    VP(:) = PhiPsi(V(:,id), p(:))

    Energy(:) = RhoVIE(:) + RhoVKE(:) + VP(:)
    F(1) = Flux(RhoV(:))
    do i = 1, dimension
      RhoVV_P(:,i) = RhoPhiU(RhoV(:), V(:,i)) + Phi(p(:)) * Normal(i)
      F(i+1) = Flux(RhoVV_P(:,i))
    enddo
    F(dimension+2) = Flux(Energy(:))
  end function KEEP4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function KEEPFVS2(id,rho,p,V,Normal,fd) result(F)
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(2), device           :: rho, p
    real(8), intent(in), dimension(2,dimension), device :: V
    real(8), intent(in), dimension(dimension), device   :: Normal
    real(8), intent(in), value                          :: fd
    real(8), dimension(dimension+2) :: F
    real(8)                         :: C, IE, KE, PG, PD
    real(8), dimension(dimension)   :: V1, V2, M
    real(8) Rho_m, P_m, PRho
    real(8), device :: V_m(dimension)
    ! FVS
    real(8) e1, e2, H1, H2, c1, c2, cave, uave, Mach
    real(8), dimension(dimension+2), device             :: Q1, Q2
    real(8), dimension(dimension+2,dimension+2), device :: A1, A2
    V1(:) = V(1,:)
    V2(:) = V(2,:)
    ! mass convection
    !C = 0.25d0 * (rho(1) + rho(2)) * (V1(id) + V2(id))
    ! momentum convection
    !M(:) = C * 0.5d0 * (V1(:) + V2(:))
    ! internal energy
    !IE = C * (p(1) / rho(1) + p(2) / rho(2)) / (2.d0 * (gamma - 1.d0))
    ! kinetic energy
    !KE = C * 0.5d0 * vecsum(V1, V2)
    ! pressure gradient
    !PG = 0.5d0 * (p(1) + p(2))
    ! pressure diffusion
    !PD = 0.5d0 * (V1(id) * p(2) + V2(id) * p(1))
    !F(1) = C
    !F(2:dimension+1) = M(:) + PG * Normal(:)
    !F(dimension+2) = IE + KE + PD
    Rho_m = 0.5d0 * (rho(1) + rho(2))
    ! contravariant velocity
    V_m(:) = 0.5d0 * (V1(:) + V2(:)) * Normal(id)
    P_m = 0.5d0 * (p(1) + p(2))
    PRho = 0.5d0 * (p(1) / rho(1) + p(2) / rho(2))
    F(1) = Rho_m * V_m(id)
    F(2:dimension+1) = F(1) * V_m(:) + P_m * Normal(:)
    F(dimension+2) = V_m(id) * P_m / (gamma - 1.d0) &
    & + 0.5d0 * F(1) * vecsum(V1, V2) &
    & + 0.5d0 * (V(1,id) * p(2) + V(2,id) * p(1)) * Normal(id)
    ! Flux Vector Spliting
    e1 = p(1) / (gamma - 1.d0) + 0.5d0 * rho(1) * vecsum(V1(:), V1(:))
    e2 = p(2) / (gamma - 1.d0) + 0.5d0 * rho(2) * vecsum(V2(:), V2(:))
    H1 = e1 + p(1) 
    H2 = e2 + p(2)
    c1 = sqrt(gamma * p(1) / rho(1))
    c2 = sqrt(gamma * p(2) / rho(2))
    cave = (sqrt(rho(1)) * c1     + sqrt(rho(2)) * c2)     / (sqrt(rho(1)) + sqrt(rho(2)))
    uave = (sqrt(rho(1)) * V1(id) + sqrt(rho(2)) * V2(id)) / (sqrt(rho(1)) + sqrt(rho(2)))
    Mach = abs(uave) / cave
    A1(:,:) = calc_AB(id,rho(1),H1,c1,V1(:))
    A2(:,:) = calc_AB(id,rho(2),H2,c2,V2(:))
    Q1(1) = rho(1)
    Q2(1) = rho(2)
    Q1(2:dimension+1) = rho(1) * V1(:)
    Q2(2:dimension+1) = rho(2) * V2(:)
    Q1(dimension+2) = e1
    Q2(dimension+2) = e2
    F(:) = F(:) - 0.5d0 * min(1.d0, Mach**2) * fd * (cumatmul(A2(:,:), Q2(:)) - cumatmul(A1(:,:), Q1(:)))
  end function KEEPFVS2

  attributes(device) function KEEPFVS4(id,rho,p,V,Normal,fd) result(F)
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(4), device           :: rho, p
    real(8), intent(in), dimension(4,dimension), device :: V
    real(8), intent(in), dimension(dimension), device   :: Normal
    real(8), intent(in), value                          :: fd
    real(8), dimension(4)           :: P_over_Rho
    real(8), dimension(dimension+2) :: F
    real(8), dimension(3)           :: RhoV, RhoVIE, RhoVKE, VP, Energy
    real(8), dimension(3,dimension) :: RhoVV_P
    integer i
    ! FVS
    real(8) e1, e2, H1, H2, c1, c2, cave, uave, Mach
    real(8), dimension(dimension), device               :: V1, V2
    real(8), dimension(dimension+2), device             :: Q1, Q2
    real(8), dimension(dimension+2,dimension+2), device :: A1, A2
    RhoV(:) = RhoPhi(rho(:), V(:,id))
    ! energy equation
    P_over_Rho(:) = p(:) / rho(:)
    RhoVIE(:) = RhoPhiU(RhoV(:), P_over_Rho(:)) / (gamma - 1.d0)
    RhoVKE(:) = RhoUPhiPhi(RhoV(:), V(:,:))
    VP(:) = PhiPsi(V(:,id), p(:))

    Energy(:) = RhoVIE(:) + RhoVKE(:) + VP(:)
    F(1) = Flux(RhoV(:))
    do i = 1, dimension
      RhoVV_P(:,i) = RhoPhiU(RhoV(:), V(:,i)) + Phi(p(:)) * Normal(i)
      F(i+1) = Flux(RhoVV_P(:,i))
    enddo
    F(dimension+2) = Flux(Energy(:))
    ! Flux Vector Spliting
    V1(:) = V(2,:)
    V2(:) = V(3,:)
    e1 = p(2) / (gamma - 1.d0) + 0.5d0 * rho(2) * vecsum(V1(:), V1(:))
    e2 = p(3) / (gamma - 1.d0) + 0.5d0 * rho(3) * vecsum(V2(:), V2(:))
    H1 = e1 + p(2)
    H2 = e2 + p(3)
    c1 = sqrt(gamma * p(2) / rho(2))
    c2 = sqrt(gamma * p(3) / rho(3))
    cave = (sqrt(rho(2)) * c1     + sqrt(rho(3)) * c2)     / (sqrt(rho(2)) + sqrt(rho(3)))
    uave = (sqrt(rho(2)) * V1(id) + sqrt(rho(3)) * V2(id)) / (sqrt(rho(2)) + sqrt(rho(3)))
    Mach = abs(uave) / cave
    A1(:,:) = calc_AB(id,rho(2),H1,c1,V1(:))
    A2(:,:) = calc_AB(id,rho(3),H2,c2,V2(:))
    Q1(1) = rho(2)
    Q2(1) = rho(3)
    Q1(2:dimension+1) = rho(2) * V(2,:)
    Q2(2:dimension+1) = rho(3) * V(3,:)
    Q1(dimension+2) = e1
    Q2(dimension+2) = e2
    F(:) = F(:) - 0.5d0 * min(1.d0, Mach**2) * fd * (cumatmul(A2(:,:), Q2(:)) - cumatmul(A1(:,:), Q1(:)))
  end function KEEPFVS4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function KEEPUP(id,rho,p,V,Normal,Ma) result(F)
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(2), device           :: rho, p
    real(8), intent(in), dimension(2,dimension), device :: V
    real(8), intent(in), dimension(dimension), device   :: Normal
    real(8), intent(in), value                          :: Ma
    real(8), dimension(dimension+2) :: F
    real(8)                         :: C, IE, KE, PG, PD, uave, phi, psi
    real(8), dimension(dimension)   :: V1, V2, M
    V1(:) = V(1,:)
    V2(:) = V(2,:)
    uave = (sqrt(rho(1)) * V1(id) + sqrt(rho(2)) * V2(id)) / (sqrt(rho(1)) + sqrt(rho(2)))
    ! blend central and upwind
    phi = 0.5d0
    psi = 0.d0!min(1.d0, Ma)
    ! mass convection
    C = 0.5d0 * (phi * (rho(1) * V1(id) + rho(2) * V2(id)) + (1.d0 - phi) * (rho(1) * V2(id) + rho(2) * V1(id)) &
    & - psi * abs(uave) * (-rho(1) + rho(2)))
    ! momentum convection
    M(:) = 0.5d0 * (C * (V1(:) + V2(:)) - psi * abs(C) * (-V1(:) + V2(:)))
    ! internal energy
    IE = (C * (p(1) / rho(1) + p(2) / rho(2)) - psi * abs(C) * (-p(1) / rho(1) + p(2) / rho(2))) / (2.d0 * (gamma - 1.d0))
    ! kinetic energy
    KE = (1.d0 - psi) * C * 0.5d0 * vecsum(V1, V2) &
    & + 0.25d0 * psi * (C * (vecsum(V1, V1) + vecsum(V2,V2)) - abs(C) * (-vecsum(V1, V1) + vecsum(V2, V2)))
    ! pressure gradient
    PG = 0.5d0 * (p(1) + p(2))
    !PG = 0.5d0 * ((p(1) + p(2)) - psi * (-p(1) + p(2)) * C / abs(C))
    ! pressure diffusion
    PD = 0.5d0 * (V1(id) * p(2) + V2(id) * p(1))
    !PD = 0.5d0 * ((1.d0 - psi) * (V1(id) * p(2) + V2(id) * p(1)) &
    !& + psi * (V1(id) * p(1) + V2(id) * p(2) - (-V1(id) * p(1) + V2(id) * p(2)) * C / abs(C)))

    F(1) = C
    F(2:dimension+1) = M(:) + PG * Normal(:)
    F(dimension+2) = IE + KE + PD
  end function KEEPUP
end module calc_keep

