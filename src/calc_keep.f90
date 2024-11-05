module calc_keep
  use mod_globals, only : dimension, gamma
  use calc_common_dim
  use calc_term
  use calc_mat
  implicit none
  interface Et2
    module procedure EtKEEP2, EtKEEPPE2, EtKEP2
  end interface

  interface Et4
    module procedure EtKEEP4, EtKEEPPE4, EtKEP4
  end interface

  interface Et6
    module procedure EtKEEP6, EtKEEPPE6, EtKEP6
  end interface
contains
  !KEEP energy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 2nd-order accuracy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function EtKEEP2(id_keep,id,C,rho,p,V1,V2) result(Et)
    integer(kind=2), intent(in), value                :: id_keep
    integer, intent(in), value                        :: id
    real(8), intent(in), value                        :: C
    real(8), intent(in), dimension(2), device         :: rho, p
    real(8), intent(in), dimension(dimension), device :: V1, V2
    real(8) :: IE, PV, KE, Et
    IE = C * 0.5d0 * (p(1) / rho(1) + p(2) / rho(2)) / (gamma - 1.d0)
    PV = 0.5d0 * (V1(id) * p(2) + V2(id) * p(1))
    KE = 0.5d0 * C * vecsum(V1, V2)
    Et = IE + PV + KE
  end function EtKEEP2

  attributes(device) function EtKEEPPE2(id_keep,id,C,rho,p,V1,V2) result(Et)
    integer(kind=4), intent(in), value                :: id_keep
    integer, intent(in), value                        :: id
    real(8), intent(in), value                        :: C
    real(8), intent(in), dimension(2), device         :: rho, p
    real(8), intent(in), dimension(dimension), device :: V1, V2
    real(8) :: IE, PV, KE, Et
    IE = 0.25d0 * (V1(id) + V2(id)) * (p(1) + p(2)) / (gamma - 1.d0)
    PV = 0.5d0 * (V1(id) * p(2) + V2(id) * p(1))
    KE = 0.5d0 * C * vecsum(V1, V2)
    Et = IE + PV + KE
  end function EtKEEPPE2

  attributes(device) function EtKEP2(id_keep,id,C,rho,p,V1,V2) result(Et)
    integer(kind=8), intent(in), value                :: id_keep
    integer, intent(in), value                        :: id
    real(8), intent(in), value                        :: C
    real(8), intent(in), dimension(2), device         :: rho, p
    real(8), intent(in), dimension(dimension), device :: V1, V2
    real(8) :: H, KE, Et
    H  = C * 0.5d0 * gamma / (gamma - 1.d0) * (p(1) / rho(1) + p(2) / rho(2))
    KE = 0.5d0 * C * vecsum(V1, V2)
    Et = H + KE
  end function EtKEP2

  ! 4th-order accuracy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function  EtKEEP4(id_keep,id,rho,p,V,RhoV) result(Et)
    integer(kind=2), intent(in), value                  :: id_keep
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(4), device           :: rho, p
    real(8), intent(in), dimension(4,dimension), device :: V
    real(8), intent(in), dimension(3), device           :: RhoV
    real(8), dimension(4) :: P_over_Rho
    real(8), dimension(3) :: RhoVIE, RhoVKE, VP, Et
    P_over_Rho(:) = p(:) / rho(:)
    RhoVIE(:)     = RhoPhiU4(RhoV(:), P_over_Rho(:)) / (gamma - 1.d0)
    RhoVKE(:)     = RhoUPhiPhi4(RhoV(:), V(:,:))
    VP(:)         = PhiPsi4(V(:,id), p(:))
    Et(:)         = RhoVIE(:) + RhoVKE(:) + VP(:)
  end function EtKEEP4

  attributes(device) function EtKEEPPE4(id_keep,id,rho,p,V,RhoV) result(Et)
    integer(kind=4), intent(in), value                  :: id_keep
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(4), device           :: rho, p
    real(8), intent(in), dimension(4,dimension), device :: V
    real(8), intent(in), dimension(3), device           :: RhoV
    real(8), dimension(3) :: Vm, RhoVIE, RhoVKE, VP, Et
    Vm(:)     = Phi4(V(:,id))
    RhoVIE(:) = RhoPhiU4(Vm(:), p(:)) / (gamma - 1.d0)    
    RhoVKE(:) = RhoUPhiPhi4(RhoV(:), V(:,:))
    VP(:)     = PhiPsi4(V(:,id), p(:))
    Et(:)     = RhoVIE(:) + RhoVKE(:) + VP(:)
  end function EtKEEPPE4

  attributes(device) function EtKEP4(id_keep,id,rho,p,V,RhoV) result(Et)
    integer(kind=8), intent(in), value                  :: id_keep
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(4), device           :: rho, p
    real(8), intent(in), dimension(4,dimension), device :: V
    real(8), intent(in), dimension(3), device           :: RhoV
    real(8), dimension(3) :: RhoVH, RhoVKE, Et
    real(8), dimension(4) :: P_over_Rho
    integer i
    P_over_Rho(:) = p(:) / rho(:)
    RhoVH(:)      = RhoPhiU4(RhoV(:), P_over_Rho(:)) * gamma / (gamma - 1.d0)
    RhoVKE(:)     = RhoUPhiPhi4(RhoV(:), V(:,:))
    Et(:)         = RhoVH(:) + RhoVKE(:)
  end function EtKEP4

  ! 6th-order accuracy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function  EtKEEP6(id_keep,id,rho,p,V,RhoV) result(Et)
    integer(kind=2), intent(in), value                  :: id_keep
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(6), device           :: rho, p
    real(8), intent(in), dimension(6,dimension), device :: V
    real(8), intent(in), dimension(6), device           :: RhoV
    real(8), dimension(6) :: P_over_Rho
    real(8), dimension(6) :: RhoVIE, RhoVKE, VP, Et
    P_over_Rho(:) = p(:) / rho(:)
    RhoVIE(:)     = RhoPhiU6(RhoV(:), P_over_Rho(:)) / (gamma - 1.d0)
    RhoVKE(:)     = RhoUPhiPhi6(RhoV(:), V(:,:))
    VP(:)         = PhiPsi6(V(:,id), p(:))
    Et(:)         = RhoVIE(:) + RhoVKE(:) + VP(:)
  end function EtKEEP6

  attributes(device) function EtKEEPPE6(id_keep,id,rho,p,V,RhoV) result(Et)
    integer(kind=4), intent(in), value                  :: id_keep
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(6), device           :: rho, p
    real(8), intent(in), dimension(6,dimension), device :: V
    real(8), intent(in), dimension(6), device           :: RhoV
    real(8), dimension(6) :: Vm, RhoVIE, RhoVKE, VP, Et
    Vm(:)     = Phi6(V(:,id))
    RhoVIE(:) = RhoPhiU6(Vm(:), p(:)) / (gamma - 1.d0)    
    RhoVKE(:) = RhoUPhiPhi6(RhoV(:), V(:,:))
    VP(:)     = PhiPsi6(V(:,id), p(:))
    Et(:)     = RhoVIE(:) + RhoVKE(:) + VP(:)
  end function EtKEEPPE6

  attributes(device) function EtKEP6(id_keep,id,rho,p,V,RhoV) result(Et)
    integer(kind=8), intent(in), value                  :: id_keep
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(6), device           :: rho, p
    real(8), intent(in), dimension(6,dimension), device :: V
    real(8), intent(in), dimension(6), device           :: RhoV
    real(8), dimension(6) :: RhoVH, RhoVKE, Et
    real(8), dimension(6) :: P_over_Rho
    integer i
    P_over_Rho(:) = p(:) / rho(:)
    RhoVH(:)      = RhoPhiU6(RhoV(:), P_over_Rho(:)) * gamma / (gamma - 1.d0)
    RhoVKE(:)     = RhoUPhiPhi6(RhoV(:), V(:,:))
    Et(:)         = RhoVH(:) + RhoVKE(:)
  end function EtKEP6

  !KEEP main!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function KEEP2(id,rho,p,V,Normal) result(F)
    use mod_globals, only : id_keep
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(2), device           :: rho, p
    real(8), intent(in), dimension(2,dimension), device :: V
    real(8), intent(in), dimension(dimension+2), device :: Normal
    real(8), dimension(dimension+2) :: F
    real(8), dimension(dimension)   :: Vm, V1, V2
    V1(:) = V(1,:)
    V2(:) = V(2,:)
    Vm(:) = 0.5d0 * (V1(:) + V2(:))
    F(1)  = 0.5d0 * (rho(1) + rho(2)) * Vm(id)
    F(2:dimension+1) = F(1) * Vm(:) + 0.5d0 * (p(1) + p(2)) * Normal(2:dimension+1)
    F(dimension+2) = Et2(id_keep,id,F(1),rho,p,V1,V2)
  end function KEEP2

  attributes(device) function KEEP4(id,rho,p,V,Normal) result(F)
    use mod_globals, only : id_keep
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(4), device           :: rho, p
    real(8), intent(in), dimension(4,dimension), device :: V
    real(8), intent(in), dimension(dimension+2), device :: Normal
    real(8), dimension(dimension+2) :: F
    real(8), dimension(3)           :: RhoV, Energy
    real(8), dimension(3,dimension) :: RhoVV_P
    integer i
    RhoV(:) = RhoPhi4(rho(:), V(:,id))
    F(1)    = Flux4(RhoV(:))
    do i = 1, dimension
      RhoVV_P(:,i) = RhoPhiU4(RhoV(:), V(:,i)) + Phi4(p(:)) * Normal(i+1)
      F(i+1)       = Flux4(RhoVV_P(:,i))
    enddo
    Energy(:) =  Et4(id_keep,id,rho,p,V,RhoV)
    F(dimension+2) = Flux4(Energy(:))
  end function KEEP4

  attributes(device) function KEEP6(id,rho,p,V,Normal) result(F)
    use mod_globals, only : id_keep
    integer, intent(in), value                          :: id
    real(8), intent(in), dimension(6), device           :: rho, p
    real(8), intent(in), dimension(6,dimension), device :: V
    real(8), intent(in), dimension(dimension+2), device :: Normal
    real(8), dimension(dimension+2) :: F
    real(8), dimension(6)           :: RhoV, Energy
    real(8), dimension(6,dimension) :: RhoVV_P
    integer i
    RhoV(:) = RhoPhi6(rho(:), V(:,id))
    F(1)    = Flux6(RhoV(:))
    do i = 1, dimension
      RhoVV_P(:,i) = RhoPhiU6(RHoV(:), V(:,i)) + Phi6(p(:)) * Normal(i+1)
      F(i+1)       = Flux6(RhoVV_P(:,i))
    enddo
    Energy(:) = Et6(id_keep,id,rho,p,V,RhoV)
    F(dimension+2) = Flux6(Energy(:))
  end function KEEP6
end module calc_keep

