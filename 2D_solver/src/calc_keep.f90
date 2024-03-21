module calc_keep
  use mod_globals, only : accuracy, offset, nx, ny, gamma
  use calc_term
  implicit none
  interface KEEP
    module procedure KEEP_2nd, KEEP_4th
  end interface
contains
  attributes(device) function KEEP_2nd(id_accuracy,id,rho,p,V,Normal) result(F)
    integer(kind=2), intent(in), value :: id_accuracy
    integer, intent(in), value :: id
    real(8), intent(in), dimension(2), device :: rho, p, Normal
    real(8), intent(in), dimension(2,2), device :: V
    real(8), dimension(4) :: F
    real(8) Rho_m, P_m, PRho
    real(8), dimension(2) :: V_m
    Rho_m = 0.5d0 * sum(rho(:))
    V_m(1) = 0.5d0 * sum(V(:,1))
    V_m(2) = 0.5d0 * sum(V(:,2))
    P_m = 0.5d0 * sum(p(:))
    PRho = 0.5d0 * sum(p(:) / rho(:))
    F(1) = Rho_m * V_m(id)
    F(2) = F(1) * V_m(1) + P_m * Normal(1)
    F(3) = F(1) * V_m(2) + P_m * Normal(2)
    F(4) = F(1) * PRho / (gamma - 1.d0) &
    & + 0.5d0 * F(1) * (V(1,1) * V(2,1) + V(1,2) * V(2,2)) &
    & + 0.5d0 * (V(1,id) * p(2) + V(2,id) * p(1))
  end function KEEP_2nd

  attributes(device) function KEEP_4th(id_accuracy,id,rho,p,V,Normal) result(F)
    integer(kind=4), intent(in), value :: id_accuracy
    integer, intent(in), value :: id
    real(8), intent(in), dimension(4), device :: rho, p
    real(8), intent(in), dimension(4,2), device :: V
    real(8), intent(in), dimension(2), device :: Normal
    real(8), dimension(4) :: F, PRho
    real(8), dimension(3) :: RhoV, RhoVU_P, RhoVV_P, IE, Energy
    RhoV(:) = RhoPhi(rho(:), V(:,id))
    RhoVU_P(:) = RhoPhiU(RhoV(:), V(:,1)) + Phi(p(:)) * Normal(1)
    RhoVV_P(:) = RHoPhiU(RhoV(:), V(:,2)) + Phi(p(:)) * Normal(2)
    PRho(:) = p(:) / rho(:)
    IE(:) = Phi(PRho(:)) / (gamma - 1.d0)
    Energy(:) = IE(:) + RhoUPhiPhi(RhoV(:), V(:,1), V(:,2)) + PhiPsi(V(:,id), p(:))
    F(1) = Flux(RhoV(:))
    F(2) = Flux(RhoVU_P(:))
    F(3) = Flux(RhoVV_P(:))
    F(4) = Flux(Energy(:))
  end function KEEP_4th

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(global) subroutine calc_E_NoMUSCL(id_accuracy, rho, u, v, p, E)
    integer(kind=accuracy), intent(in), value :: id_accuracy
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(out), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    integer i, j
    real(8), dimension(accuracy) :: rhos, ps
    real(8), dimension(2) :: Normal = (/1.d0, 0.d0/)
    real(8), dimension(accuracy,2) :: Vs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    
    rhos = rho(i:i+accuracy-1,j)
    ps = p(i:i+accuracy-1,j)
    Vs(:,1) = u(i:i+accuracy-1,j)
    Vs(:,2) = v(i:i+accuracy-1,j)
    E(i,j-offset,:) = KEEP(id_accuracy,1,rhos,ps,Vs,Normal)
  end subroutine calc_E_NoMUSCL

  attributes(global) subroutine calc_F_NoMUSCL(id_accuracy, rho, u, v, p, F)
    integer(kind=accuracy), intent(in), value :: id_accuracy
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    integer i, j
    real(8), dimension(accuracy) :: rhos, ps
    real(8), dimension(2) :: Normal = (/0.d0, 1.d0/)
    real(8), dimension(accuracy,2) :: Vs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    
    rhos = rho(i,j:j+accuracy-1)
    ps = p(i,j:j+accuracy-1)
    Vs(:,1) = u(i,j:j+accuracy-1)
    Vs(:,2) = v(i,j:j+accuracy-1)
    F(i-offset,j,:) = KEEP(id_accuracy,2,rhos,ps,Vs,Normal)
  end subroutine calc_F_NoMUSCL
end module calc_keep

