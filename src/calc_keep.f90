module calc_keep
  use mod_globals, only : dim => dimension, gamma
  use mod_constant, only : over_gamma_1
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

  attributes(device) function EtKEEP2(id_keep, id, C, rho, p, V1, V2) result(Et)
    integer(kind=2), intent(in), value  :: id_keep
    integer, intent(in), value          :: id
    real(8), intent(in), value          :: C
    real(8), intent(in), dimension(2)   :: rho, p
    real(8), intent(in), dimension(dim) :: V1, V2
    real(8) Et
    Et = C * 0.5d0 * (p(1) / rho(1) + p(2) / rho(2)) * over_gamma_1 ! internal energy
    Et = Et + 0.5d0 * (V1(id) * p(2) + V2(id) * p(1)) ! pressure diffusion
    Et = Et + 0.5d0 * C * vecsum(V1, V2) ! kinetic energy
  end function EtKEEP2

  attributes(device) function EtKEEPPE2(id_keep, id, C, rho, p, V1, V2) result(Et)
    integer(kind=4), intent(in), value  :: id_keep
    integer, intent(in), value          :: id
    real(8), intent(in), value          :: C
    real(8), intent(in), dimension(2)   :: rho, p
    real(8), intent(in), dimension(dim) :: V1, V2
    real(8) :: Et
    Et = 0.25d0 * (V1(id) + V2(id)) * (p(1) + p(2)) * over_gamma_1 ! internal energy
    Et = Et + 0.5d0 * (V1(id) * p(2) + V2(id) * p(1)) ! pressure diffusion
    Et = Et + 0.5d0 * C * vecsum(V1, V2) ! kinetic energy
  end function EtKEEPPE2

  attributes(device) function EtKEP2(id_keep, id, C, rho, p, V1, V2) result(Et)
    integer(kind=8), intent(in), value  :: id_keep
    integer, intent(in), value          :: id
    real(8), intent(in), value          :: C
    real(8), intent(in), dimension(2)   :: rho, p
    real(8), intent(in), dimension(dim) :: V1, V2
    real(8) :: Et
    Et = C * 0.5d0 * gamma * over_gamma_1 * (p(1) / rho(1) + p(2) / rho(2)) ! enthalpy
    Et = Et + 0.5d0 * C * vecsum(V1, V2) ! kinetic energy
  end function EtKEP2

  ! 4th-order accuracy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function  EtKEEP4(id_keep, id, rho, p, V, RhoV) result(Et)
    integer(kind=2), intent(in), value    :: id_keep
    integer, intent(in), value            :: id
    real(8), intent(in), dimension(4)     :: rho, p
    real(8), intent(in), dimension(4,dim) :: V
    real(8), intent(in), dimension(3)     :: RhoV
    real(8), dimension(4) :: P_over_Rho
    real(8), dimension(3) :: Et
    P_over_Rho(:) = p(:) / rho(:)
    Et(:) = RhoPhiU4(RhoV(:), P_over_Rho(:)) * over_gamma_1 ! internal energy
    Et(:) = Et(:) + RhoUPhiPhi4(RhoV(:), V(:,:)) ! kinetic energy
    Et(:) = Et(:) + PhiPsi4(V(:,id), p(:)) ! pressure diffusion
  end function EtKEEP4

  attributes(device) function EtKEEPPE4(id_keep, id, rho, p, V, RhoV) result(Et)
    integer(kind=4), intent(in), value    :: id_keep
    integer, intent(in), value            :: id
    real(8), intent(in), dimension(4)     :: rho, p
    real(8), intent(in), dimension(4,dim) :: V
    real(8), intent(in), dimension(3)     :: RhoV
    real(8), dimension(3) :: Vm, Et
    Vm(:) = Phi4(V(:,id))
    Et(:) = RhoPhiU4(Vm(:), p(:)) * over_gamma_1 ! internal energy    
    Et(:) = Et(:) + RhoUPhiPhi4(RhoV(:), V(:,:)) ! kinetic energy
    Et(:) = Et(:) + PhiPsi4(V(:,id), p(:)) ! pressure diffusion
  end function EtKEEPPE4

  attributes(device) function EtKEP4(id_keep, id, rho, p, V, RhoV) result(Et)
    integer(kind=8), intent(in), value    :: id_keep
    integer, intent(in), value            :: id
    real(8), intent(in), dimension(4)     :: rho, p
    real(8), intent(in), dimension(4,dim) :: V
    real(8), intent(in), dimension(3)     :: RhoV
    real(8), dimension(3) :: Et
    real(8), dimension(4) :: P_over_Rho
    P_over_Rho(:) = p(:) / rho(:)
    Et(:) = RhoPhiU4(RhoV(:), P_over_Rho(:)) * gamma * over_gamma_1 ! enthalpy
    Et(:) = Et(:) + RhoUPhiPhi4(RhoV(:), V(:,:)) ! kinetic energy
  end function EtKEP4

  ! 6th-order accuracy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function  EtKEEP6(id_keep, id, rho, p, V, RhoV) result(Et)
    integer(kind=2), intent(in), value    :: id_keep
    integer, intent(in), value            :: id
    real(8), intent(in), dimension(6)     :: rho, p
    real(8), intent(in), dimension(6,dim) :: V
    real(8), intent(in), dimension(6)     :: RhoV
    real(8), dimension(6) :: P_over_Rho, Et
    P_over_Rho(:) = p(:) / rho(:)
    Et(:) = RhoPhiU6(RhoV(:), P_over_Rho(:)) * over_gamma_1 ! internal energy
    Et(:) = Et(:) + RhoUPhiPhi6(RhoV(:), V(:,:)) ! kinetic energy
    Et(:) = Et(:) + PhiPsi6(V(:,id), p(:)) ! pressure diffusion
  end function EtKEEP6

  attributes(device) function EtKEEPPE6(id_keep, id, rho, p, V, RhoV) result(Et)
    integer(kind=4), intent(in), value    :: id_keep
    integer, intent(in), value            :: id
    real(8), intent(in), dimension(6)     :: rho, p
    real(8), intent(in), dimension(6,dim) :: V
    real(8), intent(in), dimension(6)     :: RhoV
    real(8), dimension(6) :: Vm, Et
    Vm(:) = Phi6(V(:,id))
    Et(:) = RhoPhiU6(Vm(:), p(:)) * over_gamma_1 ! internal energy    
    Et(:) = Et(:) + RhoUPhiPhi6(RhoV(:), V(:,:)) ! kinetic energy
    Et(:) = Et(:) + PhiPsi6(V(:,id), p(:)) ! pressure diffusion
  end function EtKEEPPE6

  attributes(device) function EtKEP6(id_keep, id, rho, p, V, RhoV) result(Et)
    integer(kind=8), intent(in), value    :: id_keep
    integer, intent(in), value            :: id
    real(8), intent(in), dimension(6)     :: rho, p
    real(8), intent(in), dimension(6,dim) :: V
    real(8), intent(in), dimension(6)     :: RhoV
    real(8), dimension(6) :: P_over_Rho, Et
    P_over_Rho(:) = p(:) / rho(:)
    Et(:) = RhoPhiU6(RhoV(:), P_over_Rho(:)) * gamma * over_gamma_1 ! enthalpy
    Et(:) = Et(:) + RhoUPhiPhi6(RhoV(:), V(:,:)) ! kinetic energy
  end function EtKEP6

  !KEEP main!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function KEEP2(id, rho, u, v, w, p, Normal) result(F)
    use mod_globals, only : id_keep
    integer, intent(in), value            :: id
    real(8), intent(in), dimension(2)     :: rho, u, v, w, p
    real(8), intent(in), dimension(dim+2) :: Normal
    real(8), dimension(dim+2) :: F
    real(8), dimension(dim)   :: Vm, V1, V2
    V1(:) = (/u(1), v(1), w(1)/)
    V2(:) = (/u(2), v(2), w(2)/)
    Vm(:) = 0.5d0 * (V1(:) + V2(:))
    F(1)  = 0.5d0 * (rho(1) + rho(2)) * Vm(id)
    F(2:dim+1) = F(1) * Vm(:) + 0.5d0 * (p(1) + p(2)) * Normal(2:dim+1)
    F(dim+2) = Et2(id_keep,id,F(1),rho,p,V1,V2)
  end function KEEP2

  attributes(device) function KEEP4(id, rho, u, v, w, p, Normal) result(F)
    use mod_globals, only : id_keep
    integer, intent(in), value            :: id
    real(8), intent(in), dimension(4)     :: rho, u, v, w, p
    real(8), intent(in), dimension(dim+2) :: Normal
    real(8), dimension(dim+2) :: F
    real(8), dimension(3)     :: RhoV
    real(8), dimension(4,dim) :: V4
    integer i
    V4(:,1) = u
    V4(:,2) = v
    V4(:,3) = w
    RhoV(:) = RhoPhi4(rho(:), V4(:,id))
    F(1)    = Flux4(RhoV(:))
    block
      real(8) RhoVV_P(3)
      do i = 1, dim
        RhoVV_P(:) = RhoPhiU4(RhoV(:), V4(:,i)) + Phi4(p(:)) * Normal(i+1)
        F(i+1)     = Flux4(RhoVV_P(:))
      enddo
    end block
    block
      real(8) Energy(3)
      Energy(:) =  Et4(id_keep, id, rho, p, V4, RhoV)
      F(dim+2) = Flux4(Energy(:))
    end block
  end function KEEP4

  attributes(device) function KEEP6(id, rho, u, v, w, p, Normal) result(F)
    use mod_globals, only : id_keep
    integer, intent(in), value            :: id
    real(8), intent(in), dimension(6)     :: rho, u, v, w, p
    real(8), intent(in), dimension(dim+2) :: Normal
    real(8), dimension(dim+2) :: F
    real(8), dimension(6)     :: RhoV
    real(8), dimension(6,dim) :: V6
    integer i
    V6(:,1) = u
    V6(:,2) = v
    V6(:,3) = w
    RhoV(:) = RhoPhi6(rho(:), V6(:,id))
    F(1)    = Flux6(RhoV(:))
    block
      real(8) RhoVV_P(6)
      do i = 1, dim
        RhoVV_P(:) = RhoPhiU6(RhoV(:), V6(:,i)) + Phi6(p(:)) * Normal(i+1)
        F(i+1)     = Flux6(RhoVV_P(:))
      enddo
    end block
    block
      real(8) Energy(6)
      Energy(:) = Et6(id_keep, id, rho, p, V6, RhoV)
      F(dim+2)  = Flux6(Energy(:))
    end block
  end function KEEP6
end module calc_keep

