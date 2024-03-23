module calc_keep
  use mod_globals, only : accuracy, offset, nx, ny, nz, gamma
  use calc_term
  implicit none
  interface KEEP
    module procedure KEEP_2nd, KEEP_4th
  end interface
contains
  attributes(device) function KEEP_2nd(id_accuracy,id,rho,p,V,Normal) result(F)
    integer(kind=2), intent(in), value :: id_accuracy
    integer, intent(in), value :: id
    real(8), intent(in), dimension(2), device :: rho, p
    real(8), intent(in), dimension(2,3), device :: V
    real(8), intent(in), dimension(3), device :: Normal
    real(8), dimension(5) :: F
    real(8) Rho_m, P_m, PRho
    real(8), dimension(3) :: V_m
    Rho_m = 0.5d0 * sum(rho(:))
    V_m(1) = 0.5d0 * sum(V(:,1))
    V_m(2) = 0.5d0 * sum(V(:,2))
    V_m(3) = 0.5d0 * sum(V(:,3))
    p_m = 0.5d0 * sum(p(:))
    PRho = 0.5d0 * sum(p(:) / rho(:))
    F(1) = Rho_m * V_m(id)
    F(2) = F(1) * V_m(1) + P_m * Normal(1)
    F(3) = F(1) * V_m(2) + P_m * Normal(2)
    F(4) = F(1) * V_m(3) + P_m * Normal(3)
    F(5) = F(1) * PRho / (gamma - 1.d0) &
    & + 0.5d0 * F(1) * (V(1,1) * V(2,1) + V(1,2) * V(2,2) + V(1,3) * V(2,3)) &
    & + 0.5d0 * (V(1,id) * p(2) + V(2,id) * p(1))
  end function KEEP_2nd

  attributes(device) function KEEP_4th(id_accuracy,id,rho,p,V,Normal) result(F)
    integer(kind=4), intent(in), value :: id_accuracy
    integer, intent(in), value :: id
    real(8), intent(in), dimension(4), device :: rho, p
    real(8), intent(in), dimension(4,3), device :: V
    real(8), intent(in), dimension(3), device :: Normal
    real(8), dimension(5) :: F
    real(8), dimension(4) :: Vr, Vx, Vy, Vz, PRho
    real(8), dimension(3) :: RhoV, RhoVU_P, RhoVV_P, RhoVW_P, IE, Energy
    Vr = V(:,id)
    Vx = V(:,1)
    Vy = V(:,2)
    Vz = V(:,3)
    RhoV(:) = RhoPhi(rho(:), Vr(:))
    RhoVU_P(:) = RhoPhiU(RhoV(:), Vx(:)) + Phi(p(:)) * Normal(1)
    RhoVV_P(:) = RhoPhiU(RhoV(:), Vy(:)) + Phi(p(:)) * Normal(2)
    RhoVW_P(:) = RhoPhiU(RhoV(:), Vz(:)) + Phi(p(:)) * Normal(3)
    PRho(:) = p(:) / rho(:)
    IE(:) = Phi(PRho(:)) / (gamma - 1.d0)
    Energy(:) = IE(:) + RhoUPhiPhi(RhoV(:), Vx(:), Vy(:), Vz(:)) + PhiPsi(Vr(:), p(:))
    F(1) = Flux(RhoV(:))
    F(2) = Flux(RhoVU_P(:))
    F(3) = Flux(RhoVV_P(:))
    F(4) = Flux(RhoVW_P(:))
    F(5) = Flux(Energy(:))
  end function KEEP_4th

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(global) subroutine calc_E(id_accuracy, rho, u, v, w, p, E)
    integer(kind=2), intent(in), value :: id_accuracy
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w, p
    real(8), intent(out), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    integer i, j, k
    real(8), dimension(accuracy) :: rho_keep, p_keep
    real(8), dimension(3) :: Normal = (/1.d0, 0.d0, 0.d0/)
    real(8), dimension(accuracy,3) :: V_keep
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    
    rho_keep = rho(i:i+accuracy-1,j,k)
    p_keep = p(i:i+accuracy-1,j,k)
    V_keep(:,1) = u(i:i+accuracy-1,j,k)
    V_keep(:,2) = v(i:i+accuracy-1,j,k)
    V_keep(:,3) = w(i:i+accuracy-1,j,k)
    E(i,j-offset,k-offset,:) = KEEP(id_accuracy,1,rho_keep,p_keep,V_keep,Normal)
  end subroutine calc_E

  attributes(global) subroutine calc_F(id_accuracy, rho, u, v, w, p, F)
    integer(kind=2), intent(in), value :: id_accuracy
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w, p
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    integer i, j, k
    real(8), dimension(accuracy) :: rho_keep, p_keep
    real(8), dimension(3) :: Normal = (/0.d0, 1.d0, 0.d0/)
    real(8), dimension(accuracy,3) :: V_keep
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset

    rho_keep = rho(i,j:j+accuracy-1,k)
    p_keep = p(i,j:j+accuracy-1,k)
    V_keep(:,1) = u(i,j:j+accuracy-1,k)
    V_keep(:,2) = v(i,j:j+accuracy-1,k)
    V_keep(:,3) = w(i,j:j+accuracy-1,k)
    F(i-offset,j,k-offset,:) = KEEP(id_accuracy,2,rho_keep,p_keep,V_keep,Normal)
  end subroutine calc_F

  attributes(global) subroutine calc_G(id_accuracy, rho, u, v, w, p, G)
    integer(kind=2), intent(in), value :: id_accuracy
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w, p
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    integer i, j, k
    real(8), dimension(accuracy) :: rho_keep, p_keep
    real(8), dimension(3) :: Normal = (/0.d0, 0.d0, 1.d0/)
    real(8), dimension(accuracy,3) :: V_keep
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z

    rho_keep = rho(i,j,k:k+accuracy-1)
    p_keep = p(i,j,k:k+accuracy-1)
    V_keep(:,1) = u(i,j,k:k+accuracy-1)
    V_keep(:,2) = v(i,j,k:k+accuracy-1)
    V_keep(:,3) = w(i,j,k:k+accuracy-1)
    G(i-offset,j-offset,k,:) = KEEP(id_accuracy,3,rho_keep,p_keep,V_keep,Normal)
  end subroutine calc_G
end module calc_keep

