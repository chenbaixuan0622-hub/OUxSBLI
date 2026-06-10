!> Curvilinear SLAU flux kernels (2nd-order, no MUSCL reconstruction).
!> E and F are area-scaled; G uses Normal_z (not area-scaled).
module calc_slau_kernel_curv
  use mod_globals, only : threadsE, threadsF, threadsG
  use mod_constant, only : over_gamma_1, id_slau
  use calc_hybrid
  implicit none
  private
  public calc_slau_xi_curv, calc_slau_eta_curv, calc_slau_z_curv
  interface SLAU
    module procedure SLAU1, HRSLAU2
  end interface SLAU
contains
  include '../../3D_solver/src/calc_slau_3d.f90'

  !> SLAU 2nd-order flux at xi-faces (i+1/2, j, k). Area-scaled by |S_xi|.
  attributes(global) subroutine calc_slau_xi_curv(id_accuracy, nx, ny, nz, n_xi_x, n_xi_y, Q, sensor, E)
    integer(2), intent(in), value               :: id_accuracy
    integer,    intent(in), value               :: nx, ny, nz
    real(8),    intent(in), device, contiguous  :: n_xi_x(nx-1,ny-2), n_xi_y(nx-1,ny-2)
    real(8),    intent(in), device, contiguous  :: Q(nx,5,ny,nz)
    real(sp),   intent(in), device, contiguous  :: sensor(nx,ny,nz)
    real(8),    intent(out), device, contiguous :: E(5,nx-1,ny-2,nz-2)
    integer :: i, j, k
    real(8) :: nxx, nxy, S, Normal(5)
    real(8) :: rhol, rhor, ul, ur, vl, vr, wl, wr, pl, pr, unl, unr
    real(sp) :: fdx
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    nxx = n_xi_x(i,j-1); nxy = n_xi_y(i,j-1)
    S   = sqrt(nxx*nxx + nxy*nxy)
    Normal = (/ 0.d0, nxx/S, nxy/S, 0.d0, 0.d0 /)
    rhol = Q(i,1,j,k); rhor = Q(i+1,1,j,k)
    ul   = Q(i,2,j,k); ur   = Q(i+1,2,j,k)
    vl   = Q(i,3,j,k); vr   = Q(i+1,3,j,k)
    wl   = Q(i,4,j,k); wr   = Q(i+1,4,j,k)
    pl   = Q(i,5,j,k); pr   = Q(i+1,5,j,k)
    unl  = ul*Normal(2) + vl*Normal(3)
    unr  = ur*Normal(2) + vr*Normal(3)
    fdx  = 0.5_sp * (sensor(i,j,k) + sensor(i+1,j,k))
    call SLAU(id_slau, rhol, rhor, ul, ur, vl, vr, wl, wr, unl, unr, pl, pr, Normal, fdx, &
              E(1,i,j-1,k-1), E(2,i,j-1,k-1), E(3,i,j-1,k-1), E(4,i,j-1,k-1), E(5,i,j-1,k-1))
    E(:,i,j-1,k-1) = E(:,i,j-1,k-1) * S
  end subroutine calc_slau_xi_curv


  !> SLAU 2nd-order flux at eta-faces (i, j+1/2, k). Area-scaled by |S_eta|.
  attributes(global) subroutine calc_slau_eta_curv(id_accuracy, nx, ny, nz, n_eta_x, n_eta_y, Q, sensor, F)
    integer(2), intent(in), value               :: id_accuracy
    integer,    intent(in), value               :: nx, ny, nz
    real(8),    intent(in), device, contiguous  :: n_eta_x(nx-2,ny-1), n_eta_y(nx-2,ny-1)
    real(8),    intent(in), device, contiguous  :: Q(nx,5,ny,nz)
    real(sp),   intent(in), device, contiguous  :: sensor(nx,ny,nz)
    real(8),    intent(out), device, contiguous :: F(5,nx-2,ny-1,nz-2)
    integer :: i, j, k
    real(8) :: nex, ney, S, Normal(5)
    real(8) :: rhol, rhor, ul, ur, vl, vr, wl, wr, pl, pr, unl, unr
    real(sp) :: fdy
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    nex = n_eta_x(i-1,j); ney = n_eta_y(i-1,j)
    S   = sqrt(nex*nex + ney*ney)
    Normal = (/ 0.d0, nex/S, ney/S, 0.d0, 0.d0 /)
    rhol = Q(i,1,j,k); rhor = Q(i,1,j+1,k)
    ul   = Q(i,2,j,k); ur   = Q(i,2,j+1,k)
    vl   = Q(i,3,j,k); vr   = Q(i,3,j+1,k)
    wl   = Q(i,4,j,k); wr   = Q(i,4,j+1,k)
    pl   = Q(i,5,j,k); pr   = Q(i,5,j+1,k)
    unl  = ul*Normal(2) + vl*Normal(3)
    unr  = ur*Normal(2) + vr*Normal(3)
    fdy  = 0.5_sp * (sensor(i,j,k) + sensor(i,j+1,k))
    call SLAU(id_slau, rhol, rhor, ul, ur, vl, vr, wl, wr, unl, unr, pl, pr, Normal, fdy, &
              F(1,i-1,j,k-1), F(2,i-1,j,k-1), F(3,i-1,j,k-1), F(4,i-1,j,k-1), F(5,i-1,j,k-1))
    F(:,i-1,j,k-1) = F(:,i-1,j,k-1) * S
  end subroutine calc_slau_eta_curv


  !> SLAU 2nd-order flux at zeta-faces. Normal_z; G not area-scaled.
  attributes(global) subroutine calc_slau_z_curv(id_accuracy, nx, ny, nz, Q, sensor, G)
    use mod_constant, only : Normal_z
    integer(2), intent(in), value               :: id_accuracy
    integer,    intent(in), value               :: nx, ny, nz
    real(8),    intent(in), device, contiguous  :: Q(nx,5,ny,nz)
    real(sp),   intent(in), device, contiguous  :: sensor(nx,ny,nz)
    real(8),    intent(out), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    integer :: i, j, k
    real(8) :: rhol, rhor, ul, ur, vl, vr, wl, wr, pl, pr
    real(sp) :: fdz
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rhol = Q(i,1,j,k); rhor = Q(i,1,j,k+1)
    ul   = Q(i,2,j,k); ur   = Q(i,2,j,k+1)
    vl   = Q(i,3,j,k); vr   = Q(i,3,j,k+1)
    wl   = Q(i,4,j,k); wr   = Q(i,4,j,k+1)
    pl   = Q(i,5,j,k); pr   = Q(i,5,j,k+1)
    fdz  = 0.5_sp * (sensor(i,j,k) + sensor(i,j,k+1))
    call SLAU(id_slau, rhol, rhor, ul, ur, vl, vr, wl, wr, wl, wr, pl, pr, Normal_z, fdz, &
              G(1,i-1,j-1,k), G(2,i-1,j-1,k), G(3,i-1,j-1,k), G(4,i-1,j-1,k), G(5,i-1,j-1,k))
  end subroutine calc_slau_z_curv
end module calc_slau_kernel_curv
