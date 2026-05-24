!> Module containing 4th-order viscous flux computation kernels
!> Uses centered difference stencils to compute viscous stresses and heat flux
!> Generally more accurate but requires larger stencils than 2nd-order
module calc_visc4
  use mod_globals, only : id_visc, id_bc_x, id_bc_y, id_bc_z, gamma, R, Pr, Prt, dt, threadsEv, threadsFv, threadsGv
  use mod_constant, only : Cp, gamma_1, Cp_over_Pr, one_third, two_third, one_twelfth
  use load_smem_visc4
  implicit none
  private
  public calc_Ev4, calc_Ev_LES4, calc_Fv4, calc_Fv_LES4, calc_Gv4, calc_Gv_LES4, &
         calc_Ev4_koff, calc_Fv4_koff, calc_Gv4_koff, calc_Gv_LES4_koff, &
         calc_Ev_LES4_koff, calc_Fv_LES4_koff
  real(8), parameter :: one_24 = 1.d0 / 24.d0 !< coefficient for 4th-order flux (1/24)
contains
  include 'calc_visc_me4_base.f90'

  !> CUDA Fortran kernel for 4th-order viscous flux in x direction
  !> High-order accurate computation of viscous stress and heat flux
  attributes(global) subroutine calc_Ev4(nx, ny, nz, dx, dy, dz, Q, T, mu, E)
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
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
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
    else
      block
        real(8) mx, mux, mvx, mwx, muy, mvy, muz, mwz
        mx  = 0.5d0 * (mu(i,j,k) + mu(i+1,j,k))
        kTx = Cp_over_Pr * mx * (-T(i,j,k) + T(i+1,j,k)) * dx(i)
        block
          real(8) :: my(2)
          my(1) = 0.25d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i+1,j-1,k) + mu(i+1,j,  k))
          my(2) = 0.25d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i+1,j,  k) + mu(i+1,j+1,k))
          muy = 0.25d0 * (my(1) * (-Q(i,2,j-1,k) + Q(i,2,j,k) - Q(i+1,2,j-1,k) + Q(i+1,2,j,k)) &
                        + my(2) * (-Q(i,2,j,k) + Q(i,2,j+1,k) - Q(i+1,2,j,k) + Q(i+1,2,j+1,k))) * dy(j)
          mvy = 0.25d0 * (my(1) * (-Q(i,3,j-1,k) + Q(i,3,j,k) - Q(i+1,3,j-1,k) + Q(i+1,3,j,k)) &
                        + my(2) * (-Q(i,3,j,k) + Q(i,3,j+1,k) - Q(i+1,3,j,k) + Q(i+1,3,j+1,k))) * dy(j)
        end block
        block
          real(8) :: mz(2)
          mz(1) = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i+1,j,k-1) + mu(i+1,j,k  ))
          mz(2) = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i+1,j,k  ) + mu(i+1,j,k+1))
          muz = 0.25d0 * (mz(1) * (-Q(i,2,j,k-1) + Q(i,2,j,k) - Q(i+1,2,j,k-1) + Q(i+1,2,j,k)) &
                        + mz(2) * (-Q(i,2,j,k) + Q(i,2,j,k+1) - Q(i+1,2,j,k) + Q(i+1,2,j,k+1))) * dz(k)
          mwz = 0.25d0 * (mz(1) * (-Q(i,4,j,k-1) + Q(i,4,j,k) - Q(i+1,4,j,k-1) + Q(i+1,4,j,k)) &
                        + mz(2) * (-Q(i,4,j,k) + Q(i,4,j,k+1) - Q(i+1,4,j,k) + Q(i+1,4,j,k+1))) * dz(k)
        end block
        mux  = mx * (-u(idx) + u(idx+1)) * dx(i)
        mvx  = mx * (-v(idx) + v(idx+1)) * dx(i)
        mwx  = mx * (-w(idx) + w(idx+1)) * dx(i)
        txx  = two_third * (2.d0 * mux - mvy - mwz)
        txy  = muy + mvx
        txz  = mwx + muz
        utxx = 0.5d0 * (u(idx) + u(idx+1)) * txx
        vtxy = 0.5d0 * (v(idx) + v(idx+1)) * txy
        wtxz = 0.5d0 * (w(idx) + w(idx+1)) * txz
      end block
    endif
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx)
  end subroutine calc_Ev4


  !> CUDA Fortran kernel for 4th-order viscous flux with LES SGS model in x direction
  attributes(global) subroutine calc_Ev_LES4(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, E)
    integer, intent(in), value                 :: nx                  !< number of grid points in x direction
    integer, intent(in), value                 :: ny                  !< number of grid points in y direction
    integer, intent(in), value                 :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous    :: dx(nx-1)            !< inverse grid spacing in x (1/dx)
    real(8), intent(in), device, contiguous    :: dy(ny-1)            !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous    :: dz(nz-1)            !< inverse grid spacing in z (1/dz)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)         !< temperature at grid points
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)        !< molecular viscosity coefficient
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)       !< turbulent eddy viscosity (LES model)
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)       !< quadratic constitutive relation correction
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2) !< viscous + SGS flux in x direction
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
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx, Hsgs
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    offset_yz = (jt-1)*sx + (kt-1)*sx*sy
    call load_smem_visc4_x(it, jt, kt, j, k, nx, ny, nz, dy, dz, Q, u, v, w, uy, vy, uz, wz)
    i  = (blockIdx%x-1)*blockDim%x + it
    idx = it + offset_yz
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8) :: mu3(3), mut3(3)
        mu3(:)  = 0.0625d0 * (9.d0 * ( mu(i-1:i+1,j,k)  + mu(i:i+2,j,k)) - ( mu(i-2:i,j,k) +  mu(i+1:i+3,j,k)))
        mut3(:) = 0.0625d0 * (9.d0 * (mut(i-1:i+1,j,k) + mut(i:i+2,j,k)) - (mut(i-2:i,j,k) + mut(i+1:i+3,j,k)))
        block ! dQdx
          real(8) :: kTx3(3)
          kTx3(:) = Cp_over_Pr * mu3(:) * &
                    (1.125d0 * (-T(i-1:i+1,j,k) + T(i:i+2,j,k)) - (-T(i-2:i,j,k) + T(i+1:i+3,j,k)) * one_24) * dx(i)
          kTx     = flux4(kTx3(:))
        end block
        call calc_tau_straight_LES(mu3, mut3, u(idx-2:idx+3), vy(idx-2:idx+3), wz(idx-2:idx+3), dx(i), txx, utxx)
        call calc_tau_cross_LES(mu3, mut3, v(idx-2:idx+3), uy(idx-2:idx+3), dx(i), txy, vtxy)
        call calc_tau_cross_LES(mu3, mut3, w(idx-2:idx+3), uz(idx-2:idx+3), dx(i), txz, wtxz)
        block
          real(8) :: H(4)
          H(:) = Cp * T(i-1:i+2,j,k) + 0.5d0 * (u(idx-1:idx+2)**2 + v(idx-1:idx+2)**2 + w(idx-1:idx+2)**2) + qc2(i-1:i+2,j,k)
          Hsgs = -flux4(mut3) * (1.125d0 * (-H(2) + H(3)) - (-H(1) + H(4)) * one_24) * dx(i) / Prt
        end block
      end block
    else
      block
        real(8), dimension(2) :: my, mysgs, mz, mzsgs
        real(8) mx, mxsgs, mux, muxsgs, mvx, mvxsgs, mwx, mwxsgs, muy, muysgs, mvy, mvysgs, muz, muzsgs, mwz, mwzsgs
        ! SGS
        mxsgs    = 0.5d0 * (mut(i,j,k) + mut(i+1,j,k))
        mysgs(:) = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i+1,j-1,k) + mut(i+1,j,k)), &
                     0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i+1,j,k) + mut(i+1,j+1,k))/)
        mzsgs(:) = (/0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i+1,j,k-1) + mut(i+1,j,k)), &
                     0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i+1,j,k) + mut(i+1,j,k+1))/)
        mx  = 0.5d0 * (mu(i,j,k) + mu(i+1,j,k))
        kTx = Cp_over_Pr * mx * (-T(i,j,k) + T(i+1,j,k)) * dx(i)
        block
          real(8) :: my(2)
          my(1)  = 0.25d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i+1,j-1,k) + mu(i+1,j,  k))
          my(2)  = 0.25d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i+1,j,  k) + mu(i+1,j+1,k))
          muy    = 0.25d0 * (my(1)    * (-Q(i,2,j-1,k) + Q(i,2,j,k) - Q(i+1,2,j-1,k) + Q(i+1,2,j,k)) &
                           + my(2)    * (-Q(i,2,j,k) + Q(i,2,j+1,k) - Q(i+1,2,j,k) + Q(i+1,2,j+1,k))) * dy(j)
          muysgs = 0.25d0 * (mysgs(1) * (-Q(i,2,j-1,k) + Q(i,2,j,k) - Q(i+1,2,j-1,k) + Q(i+1,2,j,k)) &
                           + mysgs(2) * (-Q(i,2,j,k) + Q(i,2,j+1,k) - Q(i+1,2,j,k) + Q(i+1,2,j+1,k))) * dy(j)
          mvy    = 0.25d0 * (my(1)    * (-Q(i,3,j-1,k) + Q(i,3,j,k) - Q(i+1,3,j-1,k) + Q(i+1,3,j,k)) &
                           + my(2)    * (-Q(i,3,j,k) + Q(i,3,j+1,k) - Q(i+1,3,j,k) + Q(i+1,3,j+1,k))) * dy(j)
          mvysgs = 0.25d0 * (mysgs(1) * (-Q(i,3,j-1,k) + Q(i,3,j,k) - Q(i+1,3,j-1,k) + Q(i+1,3,j,k)) &
                           + mysgs(2) * (-Q(i,3,j,k) + Q(i,3,j+1,k) - Q(i+1,3,j,k) + Q(i+1,3,j+1,k))) * dy(j)
        end block
        block
          real(8) :: mz(2)
          mz(1)  = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i+1,j,k-1) + mu(i+1,j,k  ))
          mz(2)  = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i+1,j,k  ) + mu(i+1,j,k+1))
          muz    = 0.25d0 * (mz(1)    * (-Q(i,2,j,k-1) + Q(i,2,j,k) - Q(i+1,2,j,k-1) + Q(i+1,2,j,k)) &
                           + mz(2)    * (-Q(i,2,j,k) + Q(i,2,j,k+1) - Q(i+1,2,j,k) + Q(i+1,2,j,k+1))) * dz(k)
          muzsgs = 0.25d0 * (mzsgs(1) * (-Q(i,2,j,k-1) + Q(i,2,j,k) - Q(i+1,2,j,k-1) + Q(i+1,2,j,k)) &
                           + mzsgs(2) * (-Q(i,2,j,k) + Q(i,2,j,k+1) - Q(i+1,2,j,k) + Q(i+1,2,j,k+1))) * dz(k)
          mwz    = 0.25d0 * (mz(1)    * (-Q(i,4,j,k-1) + Q(i,4,j,k) - Q(i+1,4,j,k-1) + Q(i+1,4,j,k)) &
                           + mz(2)    * (-Q(i,4,j,k) + Q(i,4,j,k+1) - Q(i+1,4,j,k) + Q(i+1,4,j,k+1))) * dz(k)
          mwzsgs = 0.25d0 * (mzsgs(1) * (-Q(i,4,j,k-1) + Q(i,4,j,k) - Q(i+1,4,j,k-1) + Q(i+1,4,j,k)) &
                           + mzsgs(2) * (-Q(i,4,j,k) + Q(i,4,j,k+1) - Q(i+1,4,j,k) + Q(i+1,4,j,k+1))) * dz(k)
        end block
        mux    = mx    * (-u(idx) + u(idx+1)) * dx(i)
        muxsgs = mxsgs * (-u(idx) + u(idx+1)) * dx(i)
        mvx    = mx    * (-v(idx) + v(idx+1)) * dx(i)
        mvxsgs = mxsgs * (-v(idx) + v(idx+1)) * dx(i)
        mwx    = mx    * (-w(idx) + w(idx+1)) * dx(i)
        mwxsgs = mxsgs * (-w(idx) + w(idx+1)) * dx(i)
        txx    = two_third * (2.d0 * mux - mvy - mwz)
        txy    = muy + mvx
        txz    = mwx + muz
        utxx   = 0.5d0 * (u(idx) + u(idx+1)) * txx
        vtxy   = 0.5d0 * (v(idx) + v(idx+1)) * txy
        wtxz   = 0.5d0 * (w(idx) + w(idx+1)) * txz
        txx    = txx + two_third * (2.d0 * muxsgs - mvysgs - mwzsgs)
        txy    = txy + muysgs + mvxsgs
        txz    = txz + mwxsgs + muzsgs
        block
          real(8) :: H(2)
          H(:) = Cp * T(i:i+1,j,k) + 0.5d0 * (u(idx:idx+1)**2 + v(idx:idx+1)**2 + w(idx:idx+1)**2) + qc2(i:i+1,j,k)
          Hsgs = -mxsgs * (-H(1) + H(2)) * dx(i) / Prt
        end block
      end block
    endif
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx + Hsgs)
  end subroutine calc_Ev_LES4


  !> CUDA Fortran kernel for 4th-order viscous flux in y direction
  attributes(global) subroutine calc_Fv4(nx, ny, nz, dy, dx, dz, Q, T, mu, F)
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
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
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
    else
      block
        real(8) my, muy, mvy, mwy, mvz, mwz, mux, mvx
        my  = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k))
        kTy = Cp_over_Pr * my * (-T(i,j,k) + T(i,j+1,k)) * dy(j)
        block
          real(8) :: mx(2)
          mx(1) = 0.25d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j+1,k) + mu(i,  j+1,k))
          mx(2) = 0.25d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j+1,k) + mu(i+1,j+1,k))
          mux = 0.25d0 * (mx(1) * (-Q(i-1,2,j,k) + Q(i,2,j,k) - Q(i-1,2,j+1,k) + Q(i,2,j+1,k)) &
                        + mx(2) * (-Q(i,2,j,k) + Q(i+1,2,j,k) - Q(i,2,j+1,k) + Q(i+1,2,j+1,k))) * dx(i)
          mvx = 0.25d0 * (mx(1) * (-Q(i-1,3,j,k) + Q(i,3,j,k) - Q(i-1,3,j+1,k) + Q(i,3,j+1,k)) &
                        + mx(2) * (-Q(i,3,j,k) + Q(i+1,3,j,k) - Q(i,3,j+1,k) + Q(i+1,3,j+1,k))) * dx(i)
        end block
        my  = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k))
        kTy = Cp_over_Pr * my * (-T(i,j,k) + T(i,j+1,k)) * dy(j)
        block
          real(8) :: mz(2)
          mz(1) = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i,j+1,k-1) + mu(i,j+1,k  ))
          mz(2) = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i,j+1,k  ) + mu(i,j+1,k+1))
          mvz = 0.25d0 * (mz(1) * (-Q(i,3,j,k-1) + Q(i,3,j,k) - Q(i,3,j+1,k-1) + Q(i,3,j+1,k)) &
                        + mz(2) * (-Q(i,3,j,k) + Q(i,3,j,k+1) - Q(i,3,j+1,k) + Q(i,3,j+1,k+1))) * dz(k)
          mwz = 0.25d0 * (mz(1) * (-Q(i,4,j,k-1) + Q(i,4,j,k) - Q(i,4,j+1,k-1) + Q(i,4,j+1,k)) &
                        + mz(2) * (-Q(i,4,j,k) + Q(i,4,j,k+1) - Q(i,4,j+1,k) + Q(i,4,j+1,k+1))) * dz(k)
        end block
        muy  = my * (-u(idx) + u(idx+1)) * dy(j)
        mvy  = my * (-v(idx) + v(idx+1)) * dy(j)
        mwy  = my * (-w(idx) + w(idx+1)) * dy(j)
        tyx  = muy + mvx
        tyy  = two_third * (2.d0 * mvy - mwz - mux)
        tyz  = mvz + mwy
        utyx = 0.5d0 * (u(idx) + u(idx+1)) * tyx
        vtyy = 0.5d0 * (v(idx) + v(idx+1)) * tyy
        wtyz = 0.5d0 * (w(idx) + w(idx+1)) * tyz
      end block
    endif
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy)
  end subroutine calc_Fv4


  !> CUDA Fortran kernel for 4th-order viscous flux with LES SGS model in y direction
  attributes(global) subroutine calc_Fv_LES4(nx, ny, nz, dy, dx, dz, Q, T, mu, mut, qc2, F)
    integer, intent(in), value                 :: nx                  !< number of grid points in x direction
    integer, intent(in), value                 :: ny                  !< number of grid points in y direction
    integer, intent(in), value                 :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous    :: dy(ny-1)            !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous    :: dx(nx-1)            !< inverse grid spacing in x (1/dx)
    real(8), intent(in), device, contiguous    :: dz(nz-1)            !< inverse grid spacing in z (1/dz)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)         !< temperature at grid points
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)        !< molecular viscosity coefficient
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)       !< turbulent eddy viscosity (LES model)
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)       !< quadratic constitutive relation correction
    real(8), intent(inout), device, contiguous :: F(5,nx-2,ny-1,nz-2) !< viscous + SGS flux in y direction
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
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, Hsgs
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    offset_xz = (it-1)*sy + (kt-1)*sy*sx
    call load_smem_visc4_y(it, jt, kt, i, k, nx, ny, nz, dx, dz, Q, u, v, w, ux, vx, vz, wz)
    j  = (blockIdx%y-1)*blockDim%y + jt
    idx = jt + offset_xz
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8) :: mu3(3), mut3(3)
        mu3(:)  = 0.0625d0 * (9.d0 * ( mu(i,j-1:j+1,k) +  mu(i,j:j+2,k)) - ( mu(i,j-2:j,k) +  mu(i,j+1:j+3,k)))
        mut3(:) = 0.0625d0 * (9.d0 * (mut(i,j-1:j+1,k) + mut(i,j:j+2,k)) - (mut(i,j-2:j,k) + mut(i,j+1:j+3,k)))
        block ! dQdy
          real(8) :: kTy3(3)
          kTy3(:) = Cp_over_Pr * mu3(:) * &
                    (1.125d0 * (-T(i,j-1:j+1,k) + T(i,j:j+2,k)) - (-T(i,j-2:j,k) + T(i,j+1:j+3,k)) * one_24) * dy(j)
          kTy     = flux4(kTy3(:))
        end block
        call calc_tau_straight_LES(mu3, mut3, v(idx-2:idx+3), wz(idx-2:idx+3), ux(idx-2:idx+3), dy(j), tyy, vtyy)
        call calc_tau_cross_LES(mu3, mut3, u(idx-2:idx+3), vx(idx-2:idx+3), dy(j), tyx, utyx)
        call calc_tau_cross_LES(mu3, mut3, w(idx-2:idx+3), vz(idx-2:idx+3), dy(j), tyz, wtyz)
        block
          real(8) :: H(4)
          H(:) = Cp * T(i,j-1:j+2,k) + 0.5d0 * (u(idx-1:idx+2)**2 + v(idx-1:idx+2)**2 + w(idx-1:idx+2)**2) + qc2(i,j-1:j+2,k)
          Hsgs = -flux4(mut3) * (1.125d0 * (-H(2) + H(3)) - (-H(1) + H(4)) * one_24) * dy(j) / Prt
        end block
      end block
    else
      block
        real(8), dimension(2) :: u2, v2, w2, mz, mzsgs, mx, mxsgs
        real(8) my, mysgs, muy, muysgs, mvy, mvysgs, mwy, mwysgs, mvz, mvzsgs, mwz, mwzsgs, mux, muxsgs, mvx, mvxsgs
        ! SGS
        mysgs    = 0.5d0 * (mut(i,j,k) + mut(i,j+1,k))
        mzsgs(:) = (/0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i,j+1,k-1) + mut(i,j+1,k)), &
                     0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i,j+1,k) + mut(i,j+1,k+1))/)
        mxsgs(:) = (/0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j+1,k) + mut(i,j+1,k)), &
                     0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j+1,k) + mut(i+1,j+1,k))/)
        block
          real(8) :: mx(2)
          mx(1)  = 0.25d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j+1,k) + mu(i,  j+1,k))
          mx(2)  = 0.25d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j+1,k) + mu(i+1,j+1,k))
          mux    = 0.25d0 * (mx(1)    * (-Q(i-1,2,j,k) + Q(i,2,j,k) - Q(i-1,2,j+1,k) + Q(i,2,j+1,k)) &
                           + mx(2)    * (-Q(i,2,j,k) + Q(i+1,2,j,k) - Q(i,2,j+1,k) + Q(i+1,2,j+1,k))) * dx(i)
          muxsgs = 0.25d0 * (mxsgs(1) * (-Q(i-1,2,j,k) + Q(i,2,j,k) - Q(i-1,2,j+1,k) + Q(i,2,j+1,k)) &
                           + mxsgs(2) * (-Q(i,2,j,k) + Q(i+1,2,j,k) - Q(i,2,j+1,k) + Q(i+1,2,j+1,k))) * dx(i)
          mvx    = 0.25d0 * (mx(1)    * (-Q(i-1,3,j,k) + Q(i,3,j,k) - Q(i-1,3,j+1,k) + Q(i,3,j+1,k)) &
                           + mx(2)    * (-Q(i,3,j,k) + Q(i+1,3,j,k) - Q(i,3,j+1,k) + Q(i+1,3,j+1,k))) * dx(i)
          mvxsgs = 0.25d0 * (mxsgs(1) * (-Q(i-1,3,j,k) + Q(i,3,j,k) - Q(i-1,3,j+1,k) + Q(i,3,j+1,k)) &
                           + mxsgs(2) * (-Q(i,3,j,k) + Q(i+1,3,j,k) - Q(i,3,j+1,k) + Q(i+1,3,j+1,k))) * dx(i)
        end block
        my  = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k))
        kTy = Cp_over_Pr * my * (-T(i,j,k) + T(i,j+1,k)) * dy(j)
        block
          real(8) :: mz(2)
          mz(1)  = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i,j+1,k-1) + mu(i,j+1,k  ))
          mz(2)  = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i,j+1,k  ) + mu(i,j+1,k+1))
          mvz    = 0.25d0 * (mz(1)    * (-Q(i,3,j,k-1) + Q(i,3,j,k) - Q(i,3,j+1,k-1) + Q(i,3,j+1,k)) &
                           + mz(2)    * (-Q(i,3,j,k) + Q(i,3,j,k+1) - Q(i,3,j+1,k) + Q(i,3,j+1,k+1))) * dz(k)
          mvzsgs = 0.25d0 * (mzsgs(1) * (-Q(i,3,j,k-1) + Q(i,3,j,k) - Q(i,3,j+1,k-1) + Q(i,3,j+1,k)) &
                           + mzsgs(2) * (-Q(i,3,j,k) + Q(i,3,j,k+1) - Q(i,3,j+1,k) + Q(i,3,j+1,k+1))) * dz(k)
          mwz    = 0.25d0 * (mz(1)    * (-Q(i,4,j,k-1) + Q(i,4,j,k) - Q(i,4,j+1,k-1) + Q(i,4,j+1,k)) &
                           + mz(2)    * (-Q(i,4,j,k) + Q(i,4,j,k+1) - Q(i,4,j+1,k) + Q(i,4,j+1,k+1))) * dz(k)
          mwzsgs = 0.25d0 * (mzsgs(1) * (-Q(i,4,j,k-1) + Q(i,4,j,k) - Q(i,4,j+1,k-1) + Q(i,4,j+1,k)) &
                           + mzsgs(2) * (-Q(i,4,j,k) + Q(i,4,j,k+1) - Q(i,4,j+1,k) + Q(i,4,j+1,k+1))) * dz(k)
        end block
        muy    = my    * (-u(idx) + u(idx+1)) * dy(j)
        muysgs = mysgs * (-u(idx) + u(idx+1)) * dy(j)
        mvy    = my    * (-v(idx) + v(idx+1)) * dy(j)
        mvysgs = mysgs * (-v(idx) + v(idx+1)) * dy(j)
        mwy    = my    * (-w(idx) + w(idx+1)) * dy(j)
        mwysgs = mysgs * (-w(idx) + w(idx+1)) * dy(j)
        tyx    = muy + mvx
        tyy    = two_third * (2.d0 * mvy - mwz - mux)
        tyz    = mvz + mwy
        utyx   = 0.5d0 * (u(idx) + u(idx+1)) * tyx
        vtyy   = 0.5d0 * (v(idx) + v(idx+1)) * tyy
        wtyz   = 0.5d0 * (w(idx) + w(idx+1)) * tyz
        tyx    = tyx + muysgs + mvxsgs
        tyy    = tyy + two_third * (2.d0 * mvysgs - mwzsgs - muxsgs)
        tyz    = tyz + mvzsgs + mwysgs
        block
          real(8) :: H(2)
          H(:) = Cp * T(i,j:j+1,k) + 0.5d0 * (u(idx:idx+1)**2 + v(idx:idx+1)**2 + w(idx:idx+1)**2) + qc2(i,j:j+1,k)
          Hsgs = -mysgs * (-H(1) + H(2)) * dy(j) / Prt
        end block
      end block
    endif
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy + Hsgs)
  end subroutine calc_Fv_LES4


  !> CUDA Fortran kernel for 4th-order viscous flux in z direction
  attributes(global) subroutine calc_Gv4(nx, ny, nz, dx, dy, dz, Q, T, mu, G)
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
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
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
    else
      block
        real(8) mz, muz, mvz, mwz, mwx, mux, mvy, mwy
        mz  = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1))
        kTz = Cp_over_Pr * mz * (-T(i,j,k) + T(i,j,k+1)) * dz(k)
        block
          real(8) :: mx(2)
          mx(1) = 0.25d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j,k+1) + mu(i,  j,k+1))
          mx(2) = 0.25d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j,k+1) + mu(i+1,j,k+1))
          mux = 0.25d0 * (mx(1) * (-Q(i-1,2,j,k) + Q(i,2,j,k) - Q(i-1,2,j,k+1) + Q(i,2,j,k+1)) &
                        + mx(2) * (-Q(i,2,j,k) + Q(i+1,2,j,k) - Q(i,2,j,k+1) + Q(i+1,2,j,k+1))) * dx(i)
          mwx = 0.25d0 * (mx(1) * (-Q(i-1,4,j,k) + Q(i,4,j,k) - Q(i-1,4,j,k+1) + Q(i,4,j,k+1)) &
                        + mx(2) * (-Q(i,4,j,k) + Q(i+1,4,j,k) - Q(i,4,j,k+1) + Q(i+1,4,j,k+1))) * dx(i)
        end block
        block
          real(8) :: my(2)
          my(1) = 0.25d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i,j-1,k+1) + mu(i,j,  k+1))
          my(2) = 0.25d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i,j,  k+1) + mu(i,j+1,k+1))
          mvy = 0.25d0 * (my(1) * (-Q(i,3,j-1,k) + Q(i,3,j,k) - Q(i,3,j-1,k+1) + Q(i,3,j,k+1)) &
                        + my(2) * (-Q(i,3,j,k) + Q(i,3,j+1,k) - Q(i,3,j,k+1) + Q(i,3,j+1,k+1))) * dy(j)
          mwy = 0.25d0 * (my(1) * (-Q(i,4,j-1,k) + Q(i,4,j,k) - Q(i,4,j-1,k+1) + Q(i,4,j,k+1)) &
                        + my(2) * (-Q(i,4,j,k) + Q(i,4,j+1,k) - Q(i,4,j,k+1) + Q(i,4,j+1,k+1))) * dy(j)
        end block
        mz   = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1))
        kTz  = Cp_over_Pr * mz * (-T(i,j,k) + T(i,j,k+1)) * dz(k)
        muz  = mz * (-u(idx) + u(idx+1)) * dz(k)
        mvz  = mz * (-v(idx) + v(idx+1)) * dz(k)
        mwz  = mz * (-w(idx) + w(idx+1)) * dz(k)
        tzx  = mwx + muz
        tzy  = mvz + mwy
        tzz  = two_third * (2.d0 * mwz - mux - mvy)
        utzx = 0.5d0 * (u(idx) + u(idx+1)) * tzx
        vtzy = 0.5d0 * (v(idx) + v(idx+1)) * tzy
        wtzz = 0.5d0 * (w(idx) + w(idx+1)) * tzz
      end block
    endif
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz)
  end subroutine calc_Gv4


  attributes(global) subroutine calc_Gv_LES4(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, G)
    integer, intent(in), value                 :: nx                  !< number of grid points in x direction
    integer, intent(in), value                 :: ny                  !< number of grid points in y direction
    integer, intent(in), value                 :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous    :: dx(nx-1)            !< inverse grid spacing in x (1/dx)
    real(8), intent(in), device, contiguous    :: dy(ny-1)            !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous    :: dz(nz-1)            !< inverse grid spacing in z (1/dz)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)         !< temperature at grid points
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)        !< molecular viscosity coefficient
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)       !< turbulent viscosity coefficient
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)       !< kinetic energy correction term for total enthalpy (0.5 * (u^2 + v^2 + w^2))
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
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, Hsgs
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    offset_xy = (jt-1)*sz + (it-1)*sz*sy
    call load_smem_visc4_z(it, jt, kt, i, j, nx, ny, nz, dx, dy, Q, u, v, w, ux, wx, vy, wy)
    k  = (blockIdx%z-1)*blockDim%z + kt
    idx = kt + offset_xy
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      block
        real(8) :: mu3(3), mut3(3)
        mu3(:)  = 0.0625d0 * (9.d0 * ( mu(i,j,k-1:k+1) +  mu(i,j,k:k+2)) - ( mu(i,j,k-2:k) +  mu(i,j,k+1:k+3)))
        mut3(:) = 0.0625d0 * (9.d0 * (mut(i,j,k-1:k+1) + mut(i,j,k:k+2)) - (mut(i,j,k-2:k) + mut(i,j,k+1:k+3)))
        block ! dQdz
          real(8) :: kTz3(3)
          kTz3(:) = Cp_over_Pr * mu3(:) * &
                    (1.125d0 * (-T(i,j,k-1:k+1) + T(i,j,k:k+2)) - (-T(i,j,k-2:k) + T(i,j,k+1:k+3)) * one_24) * dz(k)
          kTz     = flux4(kTz3(:))
        end block
        call calc_tau_straight_LES(mu3, mut3, w(idx-2:idx+3), ux(idx-2:idx+3), vy(idx-2:idx+3), dz(k), tzz, wtzz)
        call calc_tau_cross_LES(mu3, mut3, u(idx-2:idx+3), wx(idx-2:idx+3), dz(k), tzx, utzx)
        call calc_tau_cross_LES(mu3, mut3, v(idx-2:idx+3), wy(idx-2:idx+3), dz(k), tzy, vtzy)
        block
          real(8) :: H(4)
          H(:) = Cp * T(i,j,k-1:k+2) + 0.5d0 * (u(idx-1:idx+2)**2 + v(idx-1:idx+2)**2 + w(idx-1:idx+2)**2) + qc2(i,j,k-1:k+2)
          Hsgs = -flux4(mut3) * (1.125d0 * (-H(2) + H(3)) - (-H(1) + H(4)) * one_24) * dz(k) / Prt
        end block
      end block
    else
      block
        real(8), dimension(2) :: mx, mxsgs, my, mysgs
        real(8) mz, mzsgs, muz, muzsgs, mvz, mvzsgs, mwz, mwzsgs, mwx, mwxsgs, mux, muxsgs, mvy, mvysgs, mwy, mwysgs
        ! SGS
        mzsgs    = 0.5d0 * (mut(i,j,k) + mut(i,j,k+1))
        mxsgs(:) = (/0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j,k+1) + mut(i,j,k+1)), &
                     0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j,k+1) + mut(i+1,j,k+1))/)
        mysgs(:) = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i,j-1,k+1) + mut(i,j,k+1)), &
                     0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i,j,k+1) + mut(i,j+1,k+1))/)
        block
          real(8) :: mx(2)
          mx(1)  = 0.25d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j,k+1) + mu(i,  j,k+1))
          mx(2)  = 0.25d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j,k+1) + mu(i+1,j,k+1))
          mux    = 0.25d0 * (mx(1)    * (-Q(i-1,2,j,k) + Q(i,2,j,k) - Q(i-1,2,j,k+1) + Q(i,2,j,k+1)) &
                           + mx(2)    * (-Q(i,2,j,k) + Q(i+1,2,j,k) - Q(i,2,j,k+1) + Q(i+1,2,j,k+1))) * dx(i)
          muxsgs = 0.25d0 * (mxsgs(1) * (-Q(i-1,2,j,k) + Q(i,2,j,k) - Q(i-1,2,j,k+1) + Q(i,2,j,k+1)) &
                           + mxsgs(2) * (-Q(i,2,j,k) + Q(i+1,2,j,k) - Q(i,2,j,k+1) + Q(i+1,2,j,k+1))) * dx(i)
          mwx    = 0.25d0 * (mx(1)    * (-Q(i-1,4,j,k) + Q(i,4,j,k) - Q(i-1,4,j,k+1) + Q(i,4,j,k+1)) &
                           + mx(2)    * (-Q(i,4,j,k) + Q(i+1,4,j,k) - Q(i,4,j,k+1) + Q(i+1,4,j,k+1))) * dx(i)
          mwxsgs = 0.25d0 * (mxsgs(1) * (-Q(i-1,4,j,k) + Q(i,4,j,k) - Q(i-1,4,j,k+1) + Q(i,4,j,k+1)) &
                           + mxsgs(2) * (-Q(i,4,j,k) + Q(i+1,4,j,k) - Q(i,4,j,k+1) + Q(i+1,4,j,k+1))) * dx(i)
        end block
        block
          real(8) :: my(2)
          my(1)  = 0.25d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i,j-1,k+1) + mu(i,j,  k+1))
          my(2)  = 0.25d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i,j,  k+1) + mu(i,j+1,k+1))
          mvy    = 0.25d0 * (my(1)    * (-Q(i,3,j-1,k) + Q(i,3,j,k) - Q(i,3,j-1,k+1) + Q(i,3,j,k+1)) &
                           + my(2)    * (-Q(i,3,j,k) + Q(i,3,j+1,k) - Q(i,3,j,k+1) + Q(i,3,j+1,k+1))) * dy(j)
          mvysgs = 0.25d0 * (mysgs(1) * (-Q(i,3,j-1,k) + Q(i,3,j,k) - Q(i,3,j-1,k+1) + Q(i,3,j,k+1)) &
                           + mysgs(2) * (-Q(i,3,j,k) + Q(i,3,j+1,k) - Q(i,3,j,k+1) + Q(i,3,j+1,k+1))) * dy(j)
          mwy    = 0.25d0 * (my(1)    * (-Q(i,4,j-1,k) + Q(i,4,j,k) - Q(i,4,j-1,k+1) + Q(i,4,j,k+1)) &
                           + my(2)    * (-Q(i,4,j,k) + Q(i,4,j+1,k) - Q(i,4,j,k+1) + Q(i,4,j+1,k+1))) * dy(j)
          mwysgs = 0.25d0 * (mysgs(1) * (-Q(i,4,j-1,k) + Q(i,4,j,k) - Q(i,4,j-1,k+1) + Q(i,4,j,k+1)) &
                           + mysgs(2) * (-Q(i,4,j,k) + Q(i,4,j+1,k) - Q(i,4,j,k+1) + Q(i,4,j+1,k+1))) * dy(j)
        end block
        mz     = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1))
        kTz    = Cp_over_Pr * mz * (-T(i,j,k) + T(i,j,k+1)) * dz(k)
        muz    = mz    * (-u(idx) + u(idx+1)) * dz(k)
        muzsgs = mzsgs * (-u(idx) + u(idx+1)) * dz(k)
        mvz    = mz    * (-v(idx) + v(idx+1)) * dz(k)
        mvzsgs = mzsgs * (-v(idx) + v(idx+1)) * dz(k)
        mwz    = mz    * (-w(idx) + w(idx+1)) * dz(k)
        mwzsgs = mzsgs * (-w(idx) + w(idx+1)) * dz(k)
        tzx    = mwx + muz
        tzy    = mvz + mwy
        tzz    = two_third * (2.d0 * mwz - mux - mvy)
        utzx   = 0.5d0 * (u(idx) + u(idx+1)) * tzx
        vtzy   = 0.5d0 * (v(idx) + v(idx+1)) * tzy
        wtzz   = 0.5d0 * (w(idx) + w(idx+1)) * tzz
        tzx    = tzx + mwxsgs + muzsgs
        tzy    = tzy + mvzsgs + mwysgs
        tzz    = tzz + two_third * (2.d0 * mwzsgs - muxsgs - mvysgs)
        block
          real(8) :: H(2)
          H(:) = Cp * T(i,j,k:k+1) + 0.5d0 * (u(idx:idx+1)**2 + v(idx:idx+1)**2 + w(idx:idx+1)**2) + qc2(i,j,k:k+1)
          Hsgs = -mzsgs * (-H(1) + H(2)) * dz(k) / Prt
        end block
      end block
    endif
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz + Hsgs)
  end subroutine calc_Gv_LES4


  !> Like calc_Ev4 but restricted to z-range [k_lo, k_hi] for communication-overlap split
  attributes(global) subroutine calc_Ev4_koff(nx, ny, nz, dx, dy, dz, Q, T, mu, E, k_lo, k_hi)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), device, contiguous    :: dx(nx-1)
    real(8), intent(in), device, contiguous    :: dy(ny-1)
    real(8), intent(in), device, contiguous    :: dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2)
    integer, intent(in), value                 :: k_lo        !< z lower bound (global, 1-based)
    integer, intent(in), value                 :: k_hi        !< z upper bound (global, 1-based)
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
    k  = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt
    offset_yz = (jt-1)*sx + (kt-1)*sx*sy
    call load_smem_visc4_x(it, jt, kt, j, k, nx, ny, nz, dy, dz, Q, u, v, w, uy, vy, uz, wz)
    i  = (blockIdx%x-1)*blockDim%x + it
    idx = it + offset_yz
    if (nx-1 < i .or. ny-1 < j .or. k_hi < k) return
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
    else
      block
        real(8) mx, mux, mvx, mwx, muy, mvy, muz, mwz
        mx  = 0.5d0 * (mu(i,j,k) + mu(i+1,j,k))
        kTx = Cp_over_Pr * mx * (-T(i,j,k) + T(i+1,j,k)) * dx(i)
        block
          real(8) :: my(2)
          my(1) = 0.25d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i+1,j-1,k) + mu(i+1,j,  k))
          my(2) = 0.25d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i+1,j,  k) + mu(i+1,j+1,k))
          muy = 0.25d0 * (my(1) * (-Q(i,2,j-1,k) + Q(i,2,j,k) - Q(i+1,2,j-1,k) + Q(i+1,2,j,k)) &
                        + my(2) * (-Q(i,2,j,k) + Q(i,2,j+1,k) - Q(i+1,2,j,k) + Q(i+1,2,j+1,k))) * dy(j)
          mvy = 0.25d0 * (my(1) * (-Q(i,3,j-1,k) + Q(i,3,j,k) - Q(i+1,3,j-1,k) + Q(i+1,3,j,k)) &
                        + my(2) * (-Q(i,3,j,k) + Q(i,3,j+1,k) - Q(i+1,3,j,k) + Q(i+1,3,j+1,k))) * dy(j)
        end block
        block
          real(8) :: mz(2)
          mz(1) = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i+1,j,k-1) + mu(i+1,j,k  ))
          mz(2) = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i+1,j,k  ) + mu(i+1,j,k+1))
          muz = 0.25d0 * (mz(1) * (-Q(i,2,j,k-1) + Q(i,2,j,k) - Q(i+1,2,j,k-1) + Q(i+1,2,j,k)) &
                        + mz(2) * (-Q(i,2,j,k) + Q(i,2,j,k+1) - Q(i+1,2,j,k) + Q(i+1,2,j,k+1))) * dz(k)
          mwz = 0.25d0 * (mz(1) * (-Q(i,4,j,k-1) + Q(i,4,j,k) - Q(i+1,4,j,k-1) + Q(i+1,4,j,k)) &
                        + mz(2) * (-Q(i,4,j,k) + Q(i,4,j,k+1) - Q(i+1,4,j,k) + Q(i+1,4,j,k+1))) * dz(k)
        end block
        mux  = mx * (-u(idx) + u(idx+1)) * dx(i)
        mvx  = mx * (-v(idx) + v(idx+1)) * dx(i)
        mwx  = mx * (-w(idx) + w(idx+1)) * dx(i)
        txx  = two_third * (2.d0 * mux - mvy - mwz)
        txy  = muy + mvx
        txz  = mwx + muz
        utxx = 0.5d0 * (u(idx) + u(idx+1)) * txx
        vtxy = 0.5d0 * (v(idx) + v(idx+1)) * txy
        wtxz = 0.5d0 * (w(idx) + w(idx+1)) * txz
      end block
    endif
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx)
  end subroutine calc_Ev4_koff


  !> Like calc_Fv4 but restricted to z-range [k_lo, k_hi] for communication-overlap split
  attributes(global) subroutine calc_Fv4_koff(nx, ny, nz, dy, dx, dz, Q, T, mu, F, k_lo, k_hi)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), device, contiguous    :: dy(ny-1)
    real(8), intent(in), device, contiguous    :: dx(nx-1)
    real(8), intent(in), device, contiguous    :: dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: F(5,nx-2,ny-1,nz-2)
    integer, intent(in), value                 :: k_lo        !< z lower bound (global, 1-based)
    integer, intent(in), value                 :: k_hi        !< z upper bound (global, 1-based)
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
    k  = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt
    offset_xz = (it-1)*sy + (kt-1)*sy*sx
    call load_smem_visc4_y(it, jt, kt, i, k, nx, ny, nz, dx, dz, Q, u, v, w, ux, vx, vz, wz)
    j  = (blockIdx%y-1)*blockDim%y + jt
    idx = jt + offset_xz
    if (nx-1 < i .or. ny-1 < j .or. k_hi < k) return
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
    else
      block
        real(8) my, muy, mvy, mwy, mvz, mwz, mux, mvx
        my  = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k))
        kTy = Cp_over_Pr * my * (-T(i,j,k) + T(i,j+1,k)) * dy(j)
        block
          real(8) :: mx(2)
          mx(1) = 0.25d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j+1,k) + mu(i,  j+1,k))
          mx(2) = 0.25d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j+1,k) + mu(i+1,j+1,k))
          mux = 0.25d0 * (mx(1) * (-Q(i-1,2,j,k) + Q(i,2,j,k) - Q(i-1,2,j+1,k) + Q(i,2,j+1,k)) &
                        + mx(2) * (-Q(i,2,j,k) + Q(i+1,2,j,k) - Q(i,2,j+1,k) + Q(i+1,2,j+1,k))) * dx(i)
          mvx = 0.25d0 * (mx(1) * (-Q(i-1,3,j,k) + Q(i,3,j,k) - Q(i-1,3,j+1,k) + Q(i,3,j+1,k)) &
                        + mx(2) * (-Q(i,3,j,k) + Q(i+1,3,j,k) - Q(i,3,j+1,k) + Q(i+1,3,j+1,k))) * dx(i)
        end block
        my  = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k))
        kTy = Cp_over_Pr * my * (-T(i,j,k) + T(i,j+1,k)) * dy(j)
        block
          real(8) :: mz(2)
          mz(1) = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i,j+1,k-1) + mu(i,j+1,k  ))
          mz(2) = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i,j+1,k  ) + mu(i,j+1,k+1))
          mvz = 0.25d0 * (mz(1) * (-Q(i,3,j,k-1) + Q(i,3,j,k) - Q(i,3,j+1,k-1) + Q(i,3,j+1,k)) &
                        + mz(2) * (-Q(i,3,j,k) + Q(i,3,j,k+1) - Q(i,3,j+1,k) + Q(i,3,j+1,k+1))) * dz(k)
          mwz = 0.25d0 * (mz(1) * (-Q(i,4,j,k-1) + Q(i,4,j,k) - Q(i,4,j+1,k-1) + Q(i,4,j+1,k)) &
                        + mz(2) * (-Q(i,4,j,k) + Q(i,4,j,k+1) - Q(i,4,j+1,k) + Q(i,4,j+1,k+1))) * dz(k)
        end block
        muy  = my * (-u(idx) + u(idx+1)) * dy(j)
        mvy  = my * (-v(idx) + v(idx+1)) * dy(j)
        mwy  = my * (-w(idx) + w(idx+1)) * dy(j)
        tyx  = muy + mvx
        tyy  = two_third * (2.d0 * mvy - mwz - mux)
        tyz  = mvz + mwy
        utyx = 0.5d0 * (u(idx) + u(idx+1)) * tyx
        vtyy = 0.5d0 * (v(idx) + v(idx+1)) * tyy
        wtyz = 0.5d0 * (w(idx) + w(idx+1)) * tyz
      end block
    endif
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy)
  end subroutine calc_Fv4_koff


  !> Like calc_Gv4 but restricted to z-range [k_lo, k_hi] for communication-overlap split
  !> Uses load_smem_visc4_z_koff so the shared-memory loader applies the same z-offset internally.
  attributes(global) subroutine calc_Gv4_koff(nx, ny, nz, dx, dy, dz, Q, T, mu, G, k_lo, k_hi)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), device, contiguous    :: dx(nx-1)
    real(8), intent(in), device, contiguous    :: dy(ny-1)
    real(8), intent(in), device, contiguous    :: dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    integer, intent(in), value                 :: k_lo        !< z lower bound (global, 1-based)
    integer, intent(in), value                 :: k_hi        !< z upper bound (global, 1-based)
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
    call load_smem_visc4_z_koff(it, jt, kt, i, j, nx, ny, nz, dx, dy, Q, u, v, w, ux, wx, vy, wy, k_lo)
    k  = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt
    idx = kt + offset_xy
    if (nx-1 < i .or. ny-1 < j .or. k_hi < k) return
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
    else
      block
        real(8) mz, muz, mvz, mwz, mwx, mux, mvy, mwy
        mz  = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1))
        kTz = Cp_over_Pr * mz * (-T(i,j,k) + T(i,j,k+1)) * dz(k)
        block
          real(8) :: mx(2)
          mx(1) = 0.25d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j,k+1) + mu(i,  j,k+1))
          mx(2) = 0.25d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j,k+1) + mu(i+1,j,k+1))
          mux = 0.25d0 * (mx(1) * (-Q(i-1,2,j,k) + Q(i,2,j,k) - Q(i-1,2,j,k+1) + Q(i,2,j,k+1)) &
                        + mx(2) * (-Q(i,2,j,k) + Q(i+1,2,j,k) - Q(i,2,j,k+1) + Q(i+1,2,j,k+1))) * dx(i)
          mwx = 0.25d0 * (mx(1) * (-Q(i-1,4,j,k) + Q(i,4,j,k) - Q(i-1,4,j,k+1) + Q(i,4,j,k+1)) &
                        + mx(2) * (-Q(i,4,j,k) + Q(i+1,4,j,k) - Q(i,4,j,k+1) + Q(i+1,4,j,k+1))) * dx(i)
        end block
        block
          real(8) :: my(2)
          my(1) = 0.25d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i,j-1,k+1) + mu(i,j,  k+1))
          my(2) = 0.25d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i,j,  k+1) + mu(i,j+1,k+1))
          mvy = 0.25d0 * (my(1) * (-Q(i,3,j-1,k) + Q(i,3,j,k) - Q(i,3,j-1,k+1) + Q(i,3,j,k+1)) &
                        + my(2) * (-Q(i,3,j,k) + Q(i,3,j+1,k) - Q(i,3,j,k+1) + Q(i,3,j+1,k+1))) * dy(j)
          mwy = 0.25d0 * (my(1) * (-Q(i,4,j-1,k) + Q(i,4,j,k) - Q(i,4,j-1,k+1) + Q(i,4,j,k+1)) &
                        + my(2) * (-Q(i,4,j,k) + Q(i,4,j+1,k) - Q(i,4,j,k+1) + Q(i,4,j+1,k+1))) * dy(j)
        end block
        mz   = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1))
        kTz  = Cp_over_Pr * mz * (-T(i,j,k) + T(i,j,k+1)) * dz(k)
        muz  = mz * (-u(idx) + u(idx+1)) * dz(k)
        mvz  = mz * (-v(idx) + v(idx+1)) * dz(k)
        mwz  = mz * (-w(idx) + w(idx+1)) * dz(k)
        tzx  = mwx + muz
        tzy  = mvz + mwy
        tzz  = two_third * (2.d0 * mwz - mux - mvy)
        utzx = 0.5d0 * (u(idx) + u(idx+1)) * tzx
        vtzy = 0.5d0 * (v(idx) + v(idx+1)) * tzy
        wtzz = 0.5d0 * (w(idx) + w(idx+1)) * tzz
      end block
    endif
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz)
  end subroutine calc_Gv4_koff


  attributes(global) subroutine calc_Ev_LES4_koff(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, E, k_lo, k_hi)
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), device, contiguous    :: dx(nx-1), dy(ny-1), dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz), T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2)
    integer, intent(in), value                 :: k_lo, k_hi
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
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx, Hsgs
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt
    offset_yz = (jt-1)*sx + (kt-1)*sx*sy
    call load_smem_visc4_x(it, jt, kt, j, k, nx, ny, nz, dy, dz, Q, u, v, w, uy, vy, uz, wz)
    i  = (blockIdx%x-1)*blockDim%x + it
    idx = it + offset_yz
    if (nx-1 < i .or. ny-1 < j .or. k_hi < k) return
    if (3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8) :: mu3(3), mut3(3)
        mu3(:)  = 0.0625d0 * (9.d0 * ( mu(i-1:i+1,j,k)  + mu(i:i+2,j,k)) - ( mu(i-2:i,j,k) +  mu(i+1:i+3,j,k)))
        mut3(:) = 0.0625d0 * (9.d0 * (mut(i-1:i+1,j,k) + mut(i:i+2,j,k)) - (mut(i-2:i,j,k) + mut(i+1:i+3,j,k)))
        block ! dQdx
          real(8) :: kTx3(3)
          kTx3(:) = Cp_over_Pr * mu3(:) * &
                    (1.125d0 * (-T(i-1:i+1,j,k) + T(i:i+2,j,k)) - (-T(i-2:i,j,k) + T(i+1:i+3,j,k)) * one_24) * dx(i)
          kTx     = flux4(kTx3(:))
        end block
        call calc_tau_straight_LES(mu3, mut3, u(idx-2:idx+3), vy(idx-2:idx+3), wz(idx-2:idx+3), dx(i), txx, utxx)
        call calc_tau_cross_LES(mu3, mut3, v(idx-2:idx+3), uy(idx-2:idx+3), dx(i), txy, vtxy)
        call calc_tau_cross_LES(mu3, mut3, w(idx-2:idx+3), uz(idx-2:idx+3), dx(i), txz, wtxz)
        block
          real(8) :: H(4)
          H(:) = Cp * T(i-1:i+2,j,k) + 0.5d0 * (u(idx-1:idx+2)**2 + v(idx-1:idx+2)**2 + w(idx-1:idx+2)**2) + qc2(i-1:i+2,j,k)
          Hsgs = -flux4(mut3) * (1.125d0 * (-H(2) + H(3)) - (-H(1) + H(4)) * one_24) * dx(i) / Prt
        end block
      end block
    else
      block
        real(8), dimension(2) :: my, mysgs, mz, mzsgs
        real(8) mx, mxsgs, mux, muxsgs, mvx, mvxsgs, mwx, mwxsgs, muy, muysgs, mvy, mvysgs, muz, muzsgs, mwz, mwzsgs
        mxsgs    = 0.5d0 * (mut(i,j,k) + mut(i+1,j,k))
        mysgs(:) = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i+1,j-1,k) + mut(i+1,j,k)), &
                     0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i+1,j,k) + mut(i+1,j+1,k))/)
        mzsgs(:) = (/0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i+1,j,k-1) + mut(i+1,j,k)), &
                     0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i+1,j,k) + mut(i+1,j,k+1))/)
        mx  = 0.5d0 * (mu(i,j,k) + mu(i+1,j,k))
        kTx = Cp_over_Pr * mx * (-T(i,j,k) + T(i+1,j,k)) * dx(i)
        block
          real(8) :: my(2)
          my(1)  = 0.25d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i+1,j-1,k) + mu(i+1,j,  k))
          my(2)  = 0.25d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i+1,j,  k) + mu(i+1,j+1,k))
          muy    = 0.25d0 * (my(1)    * (-Q(i,2,j-1,k) + Q(i,2,j,k) - Q(i+1,2,j-1,k) + Q(i+1,2,j,k)) &
                           + my(2)    * (-Q(i,2,j,k) + Q(i,2,j+1,k) - Q(i+1,2,j,k) + Q(i+1,2,j+1,k))) * dy(j)
          muysgs = 0.25d0 * (mysgs(1) * (-Q(i,2,j-1,k) + Q(i,2,j,k) - Q(i+1,2,j-1,k) + Q(i+1,2,j,k)) &
                           + mysgs(2) * (-Q(i,2,j,k) + Q(i,2,j+1,k) - Q(i+1,2,j,k) + Q(i+1,2,j+1,k))) * dy(j)
          mvy    = 0.25d0 * (my(1)    * (-Q(i,3,j-1,k) + Q(i,3,j,k) - Q(i+1,3,j-1,k) + Q(i+1,3,j,k)) &
                           + my(2)    * (-Q(i,3,j,k) + Q(i,3,j+1,k) - Q(i+1,3,j,k) + Q(i+1,3,j+1,k))) * dy(j)
          mvysgs = 0.25d0 * (mysgs(1) * (-Q(i,3,j-1,k) + Q(i,3,j,k) - Q(i+1,3,j-1,k) + Q(i+1,3,j,k)) &
                           + mysgs(2) * (-Q(i,3,j,k) + Q(i,3,j+1,k) - Q(i+1,3,j,k) + Q(i+1,3,j+1,k))) * dy(j)
        end block
        block
          real(8) :: mz(2)
          mz(1)  = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i+1,j,k-1) + mu(i+1,j,k  ))
          mz(2)  = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i+1,j,k  ) + mu(i+1,j,k+1))
          muz    = 0.25d0 * (mz(1)    * (-Q(i,2,j,k-1) + Q(i,2,j,k) - Q(i+1,2,j,k-1) + Q(i+1,2,j,k)) &
                           + mz(2)    * (-Q(i,2,j,k) + Q(i,2,j,k+1) - Q(i+1,2,j,k) + Q(i+1,2,j,k+1))) * dz(k)
          muzsgs = 0.25d0 * (mzsgs(1) * (-Q(i,2,j,k-1) + Q(i,2,j,k) - Q(i+1,2,j,k-1) + Q(i+1,2,j,k)) &
                           + mzsgs(2) * (-Q(i,2,j,k) + Q(i,2,j,k+1) - Q(i+1,2,j,k) + Q(i+1,2,j,k+1))) * dz(k)
          mwz    = 0.25d0 * (mz(1)    * (-Q(i,4,j,k-1) + Q(i,4,j,k) - Q(i+1,4,j,k-1) + Q(i+1,4,j,k)) &
                           + mz(2)    * (-Q(i,4,j,k) + Q(i,4,j,k+1) - Q(i+1,4,j,k) + Q(i+1,4,j,k+1))) * dz(k)
          mwzsgs = 0.25d0 * (mzsgs(1) * (-Q(i,4,j,k-1) + Q(i,4,j,k) - Q(i+1,4,j,k-1) + Q(i+1,4,j,k)) &
                           + mzsgs(2) * (-Q(i,4,j,k) + Q(i,4,j,k+1) - Q(i+1,4,j,k) + Q(i+1,4,j,k+1))) * dz(k)
        end block
        mux    = mx    * (-u(idx) + u(idx+1)) * dx(i)
        muxsgs = mxsgs * (-u(idx) + u(idx+1)) * dx(i)
        mvx    = mx    * (-v(idx) + v(idx+1)) * dx(i)
        mvxsgs = mxsgs * (-v(idx) + v(idx+1)) * dx(i)
        mwx    = mx    * (-w(idx) + w(idx+1)) * dx(i)
        mwxsgs = mxsgs * (-w(idx) + w(idx+1)) * dx(i)
        txx    = two_third * (2.d0 * mux - mvy - mwz)
        txy    = muy + mvx
        txz    = mwx + muz
        utxx   = 0.5d0 * (u(idx) + u(idx+1)) * txx
        vtxy   = 0.5d0 * (v(idx) + v(idx+1)) * txy
        wtxz   = 0.5d0 * (w(idx) + w(idx+1)) * txz
        txx    = txx + two_third * (2.d0 * muxsgs - mvysgs - mwzsgs)
        txy    = txy + muysgs + mvxsgs
        txz    = txz + mwxsgs + muzsgs
        block
          real(8) :: H(2)
          H(:) = Cp * T(i:i+1,j,k) + 0.5d0 * (u(idx:idx+1)**2 + v(idx:idx+1)**2 + w(idx:idx+1)**2) + qc2(i:i+1,j,k)
          Hsgs = -mxsgs * (-H(1) + H(2)) * dx(i) / Prt
        end block
      end block
    endif
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx + Hsgs)
  end subroutine calc_Ev_LES4_koff


  attributes(global) subroutine calc_Fv_LES4_koff(nx, ny, nz, dy, dx, dz, Q, T, mu, mut, qc2, F, k_lo, k_hi)
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), device, contiguous    :: dy(ny-1), dx(nx-1), dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz), T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: F(5,nx-2,ny-1,nz-2)
    integer, intent(in), value                 :: k_lo, k_hi
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
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, Hsgs
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt
    offset_xz = (it-1)*sy + (kt-1)*sy*sx
    call load_smem_visc4_y(it, jt, kt, i, k, nx, ny, nz, dx, dz, Q, u, v, w, ux, vx, vz, wz)
    j  = (blockIdx%y-1)*blockDim%y + jt
    idx = jt + offset_xz
    if (nx-1 < i .or. ny-1 < j .or. k_hi < k) return
    if (3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8) :: mu3(3), mut3(3)
        mu3(:)  = 0.0625d0 * (9.d0 * ( mu(i,j-1:j+1,k) +  mu(i,j:j+2,k)) - ( mu(i,j-2:j,k) +  mu(i,j+1:j+3,k)))
        mut3(:) = 0.0625d0 * (9.d0 * (mut(i,j-1:j+1,k) + mut(i,j:j+2,k)) - (mut(i,j-2:j,k) + mut(i,j+1:j+3,k)))
        block ! dQdy
          real(8) :: kTy3(3)
          kTy3(:) = Cp_over_Pr * mu3(:) * &
                    (1.125d0 * (-T(i,j-1:j+1,k) + T(i,j:j+2,k)) - (-T(i,j-2:j,k) + T(i,j+1:j+3,k)) * one_24) * dy(j)
          kTy     = flux4(kTy3(:))
        end block
        call calc_tau_straight_LES(mu3, mut3, v(idx-2:idx+3), wz(idx-2:idx+3), ux(idx-2:idx+3), dy(j), tyy, vtyy)
        call calc_tau_cross_LES(mu3, mut3, u(idx-2:idx+3), vx(idx-2:idx+3), dy(j), tyx, utyx)
        call calc_tau_cross_LES(mu3, mut3, w(idx-2:idx+3), vz(idx-2:idx+3), dy(j), tyz, wtyz)
        block
          real(8) :: H(4)
          H(:) = Cp * T(i,j-1:j+2,k) + 0.5d0 * (u(idx-1:idx+2)**2 + v(idx-1:idx+2)**2 + w(idx-1:idx+2)**2) + qc2(i,j-1:j+2,k)
          Hsgs = -flux4(mut3) * (1.125d0 * (-H(2) + H(3)) - (-H(1) + H(4)) * one_24) * dy(j) / Prt
        end block
      end block
    else
      block
        real(8), dimension(2) :: u2, v2, w2, mz, mzsgs, mx, mxsgs
        real(8) my, mysgs, muy, muysgs, mvy, mvysgs, mwy, mwysgs, mvz, mvzsgs, mwz, mwzsgs, mux, muxsgs, mvx, mvxsgs
        mysgs    = 0.5d0 * (mut(i,j,k) + mut(i,j+1,k))
        mzsgs(:) = (/0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i,j+1,k-1) + mut(i,j+1,k)), &
                     0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i,j+1,k) + mut(i,j+1,k+1))/)
        mxsgs(:) = (/0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j+1,k) + mut(i,j+1,k)), &
                     0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j+1,k) + mut(i+1,j+1,k))/)
        block
          real(8) :: mx(2)
          mx(1)  = 0.25d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j+1,k) + mu(i,  j+1,k))
          mx(2)  = 0.25d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j+1,k) + mu(i+1,j+1,k))
          mux    = 0.25d0 * (mx(1)    * (-Q(i-1,2,j,k) + Q(i,2,j,k) - Q(i-1,2,j+1,k) + Q(i,2,j+1,k)) &
                           + mx(2)    * (-Q(i,2,j,k) + Q(i+1,2,j,k) - Q(i,2,j+1,k) + Q(i+1,2,j+1,k))) * dx(i)
          muxsgs = 0.25d0 * (mxsgs(1) * (-Q(i-1,2,j,k) + Q(i,2,j,k) - Q(i-1,2,j+1,k) + Q(i,2,j+1,k)) &
                           + mxsgs(2) * (-Q(i,2,j,k) + Q(i+1,2,j,k) - Q(i,2,j+1,k) + Q(i+1,2,j+1,k))) * dx(i)
          mvx    = 0.25d0 * (mx(1)    * (-Q(i-1,3,j,k) + Q(i,3,j,k) - Q(i-1,3,j+1,k) + Q(i,3,j+1,k)) &
                           + mx(2)    * (-Q(i,3,j,k) + Q(i+1,3,j,k) - Q(i,3,j+1,k) + Q(i+1,3,j+1,k))) * dx(i)
          mvxsgs = 0.25d0 * (mxsgs(1) * (-Q(i-1,3,j,k) + Q(i,3,j,k) - Q(i-1,3,j+1,k) + Q(i,3,j+1,k)) &
                           + mxsgs(2) * (-Q(i,3,j,k) + Q(i+1,3,j,k) - Q(i,3,j+1,k) + Q(i+1,3,j+1,k))) * dx(i)
        end block
        my  = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k))
        kTy = Cp_over_Pr * my * (-T(i,j,k) + T(i,j+1,k)) * dy(j)
        block
          real(8) :: mz(2)
          mz(1)  = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i,j+1,k-1) + mu(i,j+1,k  ))
          mz(2)  = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i,j+1,k  ) + mu(i,j+1,k+1))
          mvz    = 0.25d0 * (mz(1)    * (-Q(i,3,j,k-1) + Q(i,3,j,k) - Q(i,3,j+1,k-1) + Q(i,3,j+1,k)) &
                           + mz(2)    * (-Q(i,3,j,k) + Q(i,3,j,k+1) - Q(i,3,j+1,k) + Q(i,3,j+1,k+1))) * dz(k)
          mvzsgs = 0.25d0 * (mzsgs(1) * (-Q(i,3,j,k-1) + Q(i,3,j,k) - Q(i,3,j+1,k-1) + Q(i,3,j+1,k)) &
                           + mzsgs(2) * (-Q(i,3,j,k) + Q(i,3,j,k+1) - Q(i,3,j+1,k) + Q(i,3,j+1,k+1))) * dz(k)
          mwz    = 0.25d0 * (mz(1)    * (-Q(i,4,j,k-1) + Q(i,4,j,k) - Q(i,4,j+1,k-1) + Q(i,4,j+1,k)) &
                           + mz(2)    * (-Q(i,4,j,k) + Q(i,4,j,k+1) - Q(i,4,j+1,k) + Q(i,4,j+1,k+1))) * dz(k)
          mwzsgs = 0.25d0 * (mzsgs(1) * (-Q(i,4,j,k-1) + Q(i,4,j,k) - Q(i,4,j+1,k-1) + Q(i,4,j+1,k)) &
                           + mzsgs(2) * (-Q(i,4,j,k) + Q(i,4,j,k+1) - Q(i,4,j+1,k) + Q(i,4,j+1,k+1))) * dz(k)
        end block
        muy    = my    * (-u(idx) + u(idx+1)) * dy(j)
        muysgs = mysgs * (-u(idx) + u(idx+1)) * dy(j)
        mvy    = my    * (-v(idx) + v(idx+1)) * dy(j)
        mvysgs = mysgs * (-v(idx) + v(idx+1)) * dy(j)
        mwy    = my    * (-w(idx) + w(idx+1)) * dy(j)
        mwysgs = mysgs * (-w(idx) + w(idx+1)) * dy(j)
        tyx    = muy + mvx
        tyy    = two_third * (2.d0 * mvy - mwz - mux)
        tyz    = mvz + mwy
        utyx   = 0.5d0 * (u(idx) + u(idx+1)) * tyx
        vtyy   = 0.5d0 * (v(idx) + v(idx+1)) * tyy
        wtyz   = 0.5d0 * (w(idx) + w(idx+1)) * tyz
        tyx    = tyx + muysgs + mvxsgs
        tyy    = tyy + two_third * (2.d0 * mvysgs - mwzsgs - muxsgs)
        tyz    = tyz + mvzsgs + mwysgs
        block
          real(8) :: H(2)
          H(:) = Cp * T(i,j:j+1,k) + 0.5d0 * (u(idx:idx+1)**2 + v(idx:idx+1)**2 + w(idx:idx+1)**2) + qc2(i,j:j+1,k)
          Hsgs = -mysgs * (-H(1) + H(2)) * dy(j) / Prt
        end block
      end block
    endif
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy + Hsgs)
  end subroutine calc_Fv_LES4_koff


  !> Like calc_Gv_LES4 but restricted to z-range [k_lo, k_hi]; uses load_smem_visc4_z_koff.
  attributes(global) subroutine calc_Gv_LES4_koff(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, G, k_lo, k_hi)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), device, contiguous    :: dx(nx-1)
    real(8), intent(in), device, contiguous    :: dy(ny-1)
    real(8), intent(in), device, contiguous    :: dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    integer, intent(in), value                 :: k_lo
    integer, intent(in), value                 :: k_hi
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
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, Hsgs
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    offset_xy = (jt-1)*sz + (it-1)*sz*sy
    call load_smem_visc4_z_koff(it, jt, kt, i, j, nx, ny, nz, dx, dy, Q, u, v, w, ux, wx, vy, wy, k_lo)
    k  = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt
    idx = kt + offset_xy
    if (nx-1 < i .or. ny-1 < j .or. k_hi < k) return
    if (3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      block
        real(8) :: mu3(3), mut3(3)
        mu3(:)  = 0.0625d0 * (9.d0 * ( mu(i,j,k-1:k+1) +  mu(i,j,k:k+2)) - ( mu(i,j,k-2:k) +  mu(i,j,k+1:k+3)))
        mut3(:) = 0.0625d0 * (9.d0 * (mut(i,j,k-1:k+1) + mut(i,j,k:k+2)) - (mut(i,j,k-2:k) + mut(i,j,k+1:k+3)))
        block ! dQdz
          real(8) :: kTz3(3)
          kTz3(:) = Cp_over_Pr * mu3(:) * &
                    (1.125d0 * (-T(i,j,k-1:k+1) + T(i,j,k:k+2)) - (-T(i,j,k-2:k) + T(i,j,k+1:k+3)) * one_24) * dz(k)
          kTz     = flux4(kTz3(:))
        end block
        call calc_tau_straight_LES(mu3, mut3, w(idx-2:idx+3), ux(idx-2:idx+3), vy(idx-2:idx+3), dz(k), tzz, wtzz)
        call calc_tau_cross_LES(mu3, mut3, u(idx-2:idx+3), wx(idx-2:idx+3), dz(k), tzx, utzx)
        call calc_tau_cross_LES(mu3, mut3, v(idx-2:idx+3), wy(idx-2:idx+3), dz(k), tzy, vtzy)
        block
          real(8) :: H(4)
          H(:) = Cp * T(i,j,k-1:k+2) + 0.5d0 * (u(idx-1:idx+2)**2 + v(idx-1:idx+2)**2 + w(idx-1:idx+2)**2) + qc2(i,j,k-1:k+2)
          Hsgs = -flux4(mut3) * (1.125d0 * (-H(2) + H(3)) - (-H(1) + H(4)) * one_24) * dz(k) / Prt
        end block
      end block
    else
      block
        real(8), dimension(2) :: mx, mxsgs, my, mysgs
        real(8) mz, mzsgs, muz, muzsgs, mvz, mvzsgs, mwz, mwzsgs, mwx, mwxsgs, mux, muxsgs, mvy, mvysgs, mwy, mwysgs
        mzsgs    = 0.5d0 * (mut(i,j,k) + mut(i,j,k+1))
        mxsgs(:) = (/0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j,k+1) + mut(i,j,k+1)), &
                     0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j,k+1) + mut(i+1,j,k+1))/)
        mysgs(:) = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i,j-1,k+1) + mut(i,j,k+1)), &
                     0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i,j,k+1) + mut(i,j+1,k+1))/)
        block
          real(8) :: mx(2)
          mx(1)  = 0.25d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j,k+1) + mu(i,  j,k+1))
          mx(2)  = 0.25d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j,k+1) + mu(i+1,j,k+1))
          mux    = 0.25d0 * (mx(1)    * (-Q(i-1,2,j,k) + Q(i,2,j,k) - Q(i-1,2,j,k+1) + Q(i,2,j,k+1)) &
                           + mx(2)    * (-Q(i,2,j,k) + Q(i+1,2,j,k) - Q(i,2,j,k+1) + Q(i+1,2,j,k+1))) * dx(i)
          muxsgs = 0.25d0 * (mxsgs(1) * (-Q(i-1,2,j,k) + Q(i,2,j,k) - Q(i-1,2,j,k+1) + Q(i,2,j,k+1)) &
                           + mxsgs(2) * (-Q(i,2,j,k) + Q(i+1,2,j,k) - Q(i,2,j,k+1) + Q(i+1,2,j,k+1))) * dx(i)
          mwx    = 0.25d0 * (mx(1)    * (-Q(i-1,4,j,k) + Q(i,4,j,k) - Q(i-1,4,j,k+1) + Q(i,4,j,k+1)) &
                           + mx(2)    * (-Q(i,4,j,k) + Q(i+1,4,j,k) - Q(i,4,j,k+1) + Q(i+1,4,j,k+1))) * dx(i)
          mwxsgs = 0.25d0 * (mxsgs(1) * (-Q(i-1,4,j,k) + Q(i,4,j,k) - Q(i-1,4,j,k+1) + Q(i,4,j,k+1)) &
                           + mxsgs(2) * (-Q(i,4,j,k) + Q(i+1,4,j,k) - Q(i,4,j,k+1) + Q(i+1,4,j,k+1))) * dx(i)
        end block
        block
          real(8) :: my(2)
          my(1)  = 0.25d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i,j-1,k+1) + mu(i,j,  k+1))
          my(2)  = 0.25d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i,j,  k+1) + mu(i,j+1,k+1))
          mvy    = 0.25d0 * (my(1)    * (-Q(i,3,j-1,k) + Q(i,3,j,k) - Q(i,3,j-1,k+1) + Q(i,3,j,k+1)) &
                           + my(2)    * (-Q(i,3,j,k) + Q(i,3,j+1,k) - Q(i,3,j,k+1) + Q(i,3,j+1,k+1))) * dy(j)
          mvysgs = 0.25d0 * (mysgs(1) * (-Q(i,3,j-1,k) + Q(i,3,j,k) - Q(i,3,j-1,k+1) + Q(i,3,j,k+1)) &
                           + mysgs(2) * (-Q(i,3,j,k) + Q(i,3,j+1,k) - Q(i,3,j,k+1) + Q(i,3,j+1,k+1))) * dy(j)
          mwy    = 0.25d0 * (my(1)    * (-Q(i,4,j-1,k) + Q(i,4,j,k) - Q(i,4,j-1,k+1) + Q(i,4,j,k+1)) &
                           + my(2)    * (-Q(i,4,j,k) + Q(i,4,j+1,k) - Q(i,4,j,k+1) + Q(i,4,j+1,k+1))) * dy(j)
          mwysgs = 0.25d0 * (mysgs(1) * (-Q(i,4,j-1,k) + Q(i,4,j,k) - Q(i,4,j-1,k+1) + Q(i,4,j,k+1)) &
                           + mysgs(2) * (-Q(i,4,j,k) + Q(i,4,j+1,k) - Q(i,4,j,k+1) + Q(i,4,j+1,k+1))) * dy(j)
        end block
        mz     = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1))
        kTz    = Cp_over_Pr * mz * (-T(i,j,k) + T(i,j,k+1)) * dz(k)
        muz    = mz    * (-u(idx) + u(idx+1)) * dz(k)
        muzsgs = mzsgs * (-u(idx) + u(idx+1)) * dz(k)
        mvz    = mz    * (-v(idx) + v(idx+1)) * dz(k)
        mvzsgs = mzsgs * (-v(idx) + v(idx+1)) * dz(k)
        mwz    = mz    * (-w(idx) + w(idx+1)) * dz(k)
        mwzsgs = mzsgs * (-w(idx) + w(idx+1)) * dz(k)
        tzx    = mwx + muz
        tzy    = mvz + mwy
        tzz    = two_third * (2.d0 * mwz - mux - mvy)
        utzx   = 0.5d0 * (u(idx) + u(idx+1)) * tzx
        vtzy   = 0.5d0 * (v(idx) + v(idx+1)) * tzy
        wtzz   = 0.5d0 * (w(idx) + w(idx+1)) * tzz
        tzx    = tzx + mwxsgs + muzsgs
        tzy    = tzy + mvzsgs + mwysgs
        tzz    = tzz + two_third * (2.d0 * mwzsgs - muxsgs - mvysgs)
        block
          real(8) :: H(2)
          H(:) = Cp * T(i,j,k:k+1) + 0.5d0 * (u(idx:idx+1)**2 + v(idx:idx+1)**2 + w(idx:idx+1)**2) + qc2(i,j,k:k+1)
          Hsgs = -mzsgs * (-H(1) + H(2)) * dz(k) / Prt
        end block
      end block
    endif
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz + Hsgs)
  end subroutine calc_Gv_LES4_koff
end module calc_visc4
