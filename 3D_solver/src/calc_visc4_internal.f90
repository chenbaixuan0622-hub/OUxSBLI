!> Module containing 4th-order viscous flux computation kernels for interior-only regions
!> Optimized variants without boundary condition handling for periodic domains
!> Uses centered difference stencils; assumes all threads execute interior 4th-order path
module calc_visc4_internal
  use mod_globals, only : id_visc, gamma, R, Pr, Prt, dt, threadsEv, threadsFv, threadsGv
  use mod_constant, only : Cp, gamma_1, Cp_over_Pr, one_third, two_third, one_twelfth
  use load_smem_visc4
  implicit none
  private
  public calc_Ev4_in, calc_Fv4_in, calc_Gv4_in
  real(8), parameter :: one_24 = 1.d0 / 24.d0 !< coefficient for 4th-order flux (1/24)
contains
  include 'calc_visc_me4_base.f90'

  !> CUDA Fortran kernel for 4th-order viscous flux in x direction (interior only)
  !> High-order accurate computation of viscous stress and heat flux
  !> No boundary condition handling — all threads execute interior 4th-order stencil
  attributes(global) subroutine calc_Ev4_in(nx, ny, nz, dx, dy, dz, Q, T, mu, E)
    integer, intent(in), value                 :: nx                  !< number of grid points in x direction
    integer, intent(in), value                 :: ny                  !< number of grid points in y direction
    integer, intent(in), value                 :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous    :: dx(nx-1)            !< inverse grid spacing in x (1/dx)
    real(8), intent(in), device, contiguous    :: dy(ny-1)            !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous    :: dz(nz-1)            !< inverse grid spacing in z (1/dz)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)         !< temperature at grid points
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)        !< molecular viscosity coefficient
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2) !< viscous flux components in x direction
    integer, parameter :: io_v = 2
    integer, parameter :: sx = threadsEv%x + 2*io_v + 1
    integer, parameter :: sy = threadsEv%y
    integer, parameter :: sz = threadsEv%z
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared ::  u
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared ::  v
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared ::  w
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: uy
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: vy
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: uz
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: wz
    integer i, j, k, it, jt, kt, ii, idx, offset_yz
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    offset_yz = (jt-1)*sx + (kt-1)*sx*sy
    call load_smem_visc4_x(it, jt, kt, j, k, nx, ny, nz, dy, dz, Q, u, v, w, uy, vy, uz, wz)
    i  = (blockIdx%x-1)*blockDim%x + it
    idx = it + offset_yz
    if (3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8) :: mu3(3)
        mu3(:) = 0.0625d0 * (9.d0 * (mu(i-1:i+1,j,k) + mu(i:i+2,j,k)) - (mu(i-2:i,j,k) + mu(i+1:i+3,j,k)))
        block ! dQdx
          real(8) :: kTx3(3)
          kTx3(:) = Cp_over_Pr * mu3(:) * &
                    (1.125d0 * (-T(i-1:i+1,j,k) + T(i:i+2,j,k)) - (-T(i-2:i,j,k) + T(i+1:i+3,j,k)) * one_24) * dx(i)
          kTx     = flux4(kTx3(:))
        end block
        call calc_tau_straight(mu3, u(idx-2:idx+3), vy(idx-2:idx+3), wz(idx-2:idx+3), dx(i), txx, utxx)
        call calc_tau_cross(mu3, v(idx-2:idx+3), uy(idx-2:idx+3), dx(i), txy, vtxy)
        call calc_tau_cross(mu3, w(idx-2:idx+3), uz(idx-2:idx+3), dx(i), txz, wtxz)
      end block
      E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
      E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
      E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
      E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx)
    endif
  end subroutine calc_Ev4_in

  !> CUDA Fortran kernel for 4th-order viscous flux in y direction (interior only)
  !> High-order accurate computation of viscous stress and heat flux
  !> No boundary condition handling — all threads execute interior 4th-order stencil
  attributes(global) subroutine calc_Fv4_in(nx, ny, nz, dy, dx, dz, Q, T, mu, F)
    integer, intent(in), value                 :: nx                  !< number of grid points in x direction
    integer, intent(in), value                 :: ny                  !< number of grid points in y direction
    integer, intent(in), value                 :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous    :: dy(ny-1)            !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous    :: dx(nx-1)            !< inverse grid spacing in x (1/dx)
    real(8), intent(in), device, contiguous    :: dz(nz-1)            !< inverse grid spacing in z (1/dz)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)         !< temperature at grid points
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)        !< molecular viscosity coefficient
    real(8), intent(inout), device, contiguous :: F(5,nx-2,ny-1,nz-2) !< viscous flux components in y direction
    integer, parameter :: io_v = 2
    integer, parameter :: sx = threadsFv%x
    integer, parameter :: sy = threadsFv%y + 2*io_v + 1
    integer, parameter :: sz = threadsFv%z
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared ::  u
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared ::  v
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared ::  w
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: ux
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: vx
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: vz
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: wz
    integer i, j, k, it, jt, kt, jj, idx, offset_xz
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    offset_xz = (it-1)*sy + (kt-1)*sy*sx
    call load_smem_visc4_y(it, jt, kt, i, k, nx, ny, nz, dx, dz, Q, u, v, w, ux, vx, vz, wz)
    j  = (blockIdx%y-1)*blockDim%y + jt
    idx = jt + offset_xz
    if (3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8) :: mu3(3)
        mu3(:) = 0.0625d0 * (9.d0 * (mu(i,j-1:j+1,k) + mu(i,j:j+2,k)) - (mu(i,j-2:j,k) + mu(i,j+1:j+3,k)))
        block ! dQdy
          real(8) :: kTy3(3)
          kTy3(:) = Cp_over_Pr * mu3(:) * &
                    (1.125d0 * (-T(i,j-1:j+1,k) + T(i,j:j+2,k)) - (-T(i,j-2:j,k) + T(i,j+1:j+3,k)) * one_24) * dy(j)
          kTy     = flux4(kTy3(:))
        end block
        call calc_tau_straight(mu3, v(idx-2:idx+3), wz(idx-2:idx+3), ux(idx-2:idx+3), dy(j), tyy, vtyy)
        call calc_tau_cross(mu3, u(idx-2:idx+3), vx(idx-2:idx+3), dy(j), tyx, utyx)
        call calc_tau_cross(mu3, w(idx-2:idx+3), vz(idx-2:idx+3), dy(j), tyz, wtyz)
      end block
      F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
      F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
      F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
      F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy)
    endif
  end subroutine calc_Fv4_in

  !> CUDA Fortran kernel for 4th-order viscous flux in z direction (interior only)
  !> High-order accurate computation of viscous stress and heat flux
  !> No boundary condition handling — all threads execute interior 4th-order stencil
  attributes(global) subroutine calc_Gv4_in(nx, ny, nz, dx, dy, dz, Q, T, mu, G)
    integer, intent(in), value                 :: nx                  !< number of grid points in x direction
    integer, intent(in), value                 :: ny                  !< number of grid points in y direction
    integer, intent(in), value                 :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous    :: dx(nx-1)            !< inverse grid spacing in x (1/dx)
    real(8), intent(in), device, contiguous    :: dy(ny-1)            !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous    :: dz(nz-1)            !< inverse grid spacing in z (1/dz)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)         !< temperature at grid points
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)        !< molecular viscosity coefficient
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1) !< viscous flux components in z direction
    integer, parameter :: io_v = 2
    integer, parameter :: sx = threadsGv%x
    integer, parameter :: sy = threadsGv%y
    integer, parameter :: sz = threadsGv%z + 2*io_v + 1
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared ::  u
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared ::  v
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared ::  w
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: wx
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: wy
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: ux
    real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: vy
    integer i, j, k, it, jt, kt, kk, idx, offset_xy
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    offset_xy = (jt-1)*sz + (it-1)*sz*sy
    call load_smem_visc4_z(it, jt, kt, i, j, nx, ny, nz, dx, dy, Q, u, v, w, ux, wx, vy, wy)
    k  = (blockIdx%z-1)*blockDim%z + kt
    idx = kt + offset_xy
    if (3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      block
        real(8) :: mu3(3)
        mu3(:) = 0.0625d0 * (9.d0 * (mu(i,j,k-1:k+1) + mu(i,j,k:k+2)) - (mu(i,j,k-2:k) + mu(i,j,k+1:k+3)))
        block ! dQdz
          real(8) :: kTz3(3)
          kTz3(:) = Cp_over_Pr * mu3(:) * &
                    (1.125d0 * (-T(i,j,k-1:k+1) + T(i,j,k:k+2)) - (-T(i,j,k-2:k) + T(i,j,k+1:k+3)) * one_24) * dz(k)
          kTz     = flux4(kTz3(:))
        end block
        call calc_tau_straight(mu3, w(idx-2:idx+3), ux(idx-2:idx+3), vy(idx-2:idx+3), dz(k), tzz, wtzz)
        call calc_tau_cross(mu3, u(idx-2:idx+3), wx(idx-2:idx+3), dz(k), tzx, utzx)
        call calc_tau_cross(mu3, v(idx-2:idx+3), wy(idx-2:idx+3), dz(k), tzy, vtzy)
      end block
      G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
      G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
      G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
      G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz)
    endif
  end subroutine calc_Gv4_in
end module calc_visc4_internal
