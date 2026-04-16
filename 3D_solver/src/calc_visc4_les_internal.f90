!> Module containing 4th-order viscous flux computation kernels with LES for interior-only regions
!> Optimized variants without boundary condition handling for periodic domains with SGS modeling
!> Uses centered difference stencils; assumes all threads execute interior 4th-order path
module calc_visc4_les_internal
  use mod_globals, only : id_visc, gamma, R, Pr, Prt, dt, threadsEv, threadsFv, threadsGv
  use mod_constant, only : Cp, gamma_1, Cp_over_Pr, one_third, two_third, one_twelfth
  implicit none
  private
  public calc_Ev_LES4_in, calc_Fv_LES4_in, calc_Gv_LES4_in
  real(8), parameter :: one_24 = 1.d0 / 24.d0 !< coefficient for 4th-order flux (1/24)
contains
  include 'calc_visc_me4_base.f90'

  !> CUDA Fortran kernel for 4th-order viscous flux with LES SGS model in x direction (interior only)
  !> Combines molecular viscosity (Navier-Stokes) with turbulent viscosity (from Smagorinsky model)
  !> No boundary condition handling — all threads execute interior 4th-order stencil
  attributes(global) subroutine calc_Ev_LES4_in(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, E)
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
    real(8), intent(inout), device, contiguous :: E(nx-1,5,ny-2,nz-2) !< viscous + SGS flux in x direction
    real(8), shared ::  u(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)
    real(8), shared ::  v(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)
    real(8), shared ::  w(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)
    real(8), shared :: uy(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)
    real(8), shared :: vy(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)
    real(8), shared :: uz(-2:threadsEv%x+3,threadsEv%z,threadsEv%y)
    real(8), shared :: wz(-2:threadsEv%x+3,threadsEv%z,threadsEv%y)
    integer i, j, k, it, jt, kt, ii, i_base
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx, Hsgs
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    do ii = it-2, threadsEv%x+3, blockDim%x
      i = i_base + ii
      if (1 <= i .and. i <= nx .and. j <= ny .and. k <= nz) then
        u(ii,jt,kt) = Q(i,2,j,k)
        v(ii,jt,kt) = Q(i,3,j,k)
        w(ii,jt,kt) = Q(i,4,j,k)
      endif
      if (1 <= i .and. i <= nx .and. 3 <= j .and. j <= ny-2 .and. k <= nz) then
        uy(ii,jt,kt) = (two_third * (-Q(i,2,j-1,k) + Q(i,2,j+1,k)) - one_twelfth * (-Q(i,2,j-2,k) + Q(i,2,j+2,k))) * dy(j)
        vy(ii,jt,kt) = (two_third * (-Q(i,3,j-1,k) + Q(i,3,j+1,k)) - one_twelfth * (-Q(i,3,j-2,k) + Q(i,3,j+2,k))) * dy(j)
      endif
      if (1 <= i .and. i <= nx .and. j <= ny .and. 3 <= k .and. k <= nz-2) then
        uz(ii,jt,kt) = (two_third * (-Q(i,2,j,k-1) + Q(i,2,j,k+1)) - one_twelfth * (-Q(i,2,j,k-2) + Q(i,2,j,k+2))) * dz(k)
        wz(ii,jt,kt) = (two_third * (-Q(i,4,j,k-1) + Q(i,4,j,k+1)) - one_twelfth * (-Q(i,4,j,k-2) + Q(i,4,j,k+2))) * dz(k)
      endif
    enddo
    call syncthreads()
    i  = (blockIdx%x-1)*blockDim%x + it
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
        call calc_tau_straight_LES(mu3, mut3, u(it-2:it+3,jt,kt), vy(it-2:it+3,jt,kt), wz(it-2:it+3,jt,kt), dx(i), txx, utxx)
        call calc_tau_cross_LES(mu3, mut3, v(it-2:it+3,jt,kt), uy(it-2:it+3,jt,kt), dx(i), txy, vtxy)
        call calc_tau_cross_LES(mu3, mut3, w(it-2:it+3,jt,kt), uz(it-2:it+3,jt,kt), dx(i), txz, wtxz)
        block
          real(8) :: H(4)
          H(:) = Cp * T(i-1:i+2,j,k) + 0.5d0 * (u(it-1:it+2,jt,kt)**2 + v(it-1:it+2,jt,kt)**2 + w(it-1:it+2,jt,kt)**2) + qc2(i-1:i+2,j,k)
          Hsgs = -flux4(mut3) * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) * one_third) * dx(i) / Prt
        end block
      end block
      E(i,2,j-1,k-1) = E(i,2,j-1,k-1) - txx
      E(i,3,j-1,k-1) = E(i,3,j-1,k-1) - txy
      E(i,4,j-1,k-1) = E(i,4,j-1,k-1) - txz
      E(i,5,j-1,k-1) = E(i,5,j-1,k-1) - (utxx + vtxy + wtxz + kTx + Hsgs)
    endif
  end subroutine calc_Ev_LES4_in

  !> CUDA Fortran kernel for 4th-order viscous flux with LES SGS model in y direction (interior only)
  !> No boundary condition handling — all threads execute interior 4th-order stencil
  attributes(global) subroutine calc_Fv_LES4_in(nx, ny, nz, dy, dx, dz, Q, T, mu, mut, qc2, F)
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
    real(8), intent(inout), device, contiguous :: F(nx-2,5,ny-1,nz-2) !< viscous + SGS flux in y direction
    real(8), shared ::  u(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)
    real(8), shared ::  v(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)
    real(8), shared ::  w(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)
    real(8), shared :: ux(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)
    real(8), shared :: vx(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)
    real(8), shared :: vz(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)
    real(8), shared :: wz(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)
    integer i, j, k, it, jt, kt, jj, j_base
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, Hsgs
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    do jj = jt-2, threadsFv%y+3, blockDim%y
      j = j_base + jj
      if (i <= nx .and. 1 <= j .and. j <= ny .and. k <= nz) then
        u(jj,it,kt) = Q(i,2,j,k)
        v(jj,it,kt) = Q(i,3,j,k)
        w(jj,it,kt) = Q(i,4,j,k)
      endif
      if (3 <= i .and. i <= nx-2 .and. 1 <= j .and. j <= ny .and. k <= nz) then
        ux(jj,it,kt) = (two_third * (-Q(i-1,2,j,k) + Q(i+1,2,j,k)) - one_twelfth * (-Q(i-2,2,j,k) + Q(i+2,2,j,k))) * dx(i)
        vx(jj,it,kt) = (two_third * (-Q(i-1,3,j,k) + Q(i+1,3,j,k)) - one_twelfth * (-Q(i-2,3,j,k) + Q(i+2,3,j,k))) * dx(i)
      endif
      if (i <= nx .and. 1 <= j .and. j <= ny .and. 3 <= k .and. k <= nz-2) then
        vz(jj,it,kt) = (two_third * (-Q(i,3,j,k-1) + Q(i,3,j,k+1)) - one_twelfth * (-Q(i,3,j,k-2) + Q(i,3,j,k+2))) * dz(k)
        wz(jj,it,kt) = (two_third * (-Q(i,4,j,k-1) + Q(i,4,j,k+1)) - one_twelfth * (-Q(i,4,j,k-2) + Q(i,4,j,k+2))) * dz(k)
      endif
    enddo
    call syncthreads()
    j  = (blockIdx%y-1)*blockDim%y + jt
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
        call calc_tau_straight_LES(mu3, mut3, v(jt-2:jt+3,it,kt), wz(jt-2:jt+3,it,kt), ux(jt-2:jt+3,it,kt), dy(j), tyy, vtyy)
        call calc_tau_cross_LES(mu3, mut3, u(jt-2:jt+3,it,kt), vx(jt-2:jt+3,it,kt), dy(j), tyx, utyx)
        call calc_tau_cross_LES(mu3, mut3, w(jt-2:jt+3,it,kt), vz(jt-2:jt+3,it,kt), dy(j), tyz, wtyz)
        block
          real(8) :: H(4)
          H(:) = Cp * T(i,j-1:j+2,k) + 0.5d0 * (u(jt-1:jt+2,it,kt)**2 + v(jt-1:jt+2,it,kt)**2 + w(jt-1:jt+2,it,kt)**2) + qc2(i,j-1:j+2,k)
          Hsgs = -flux4(mut3) * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) * one_third) * dy(j) / Prt
        end block
      end block
    F(i-1,2,j,k-1) = F(i-1,2,j,k-1) - tyx
    F(i-1,3,j,k-1) = F(i-1,3,j,k-1) - tyy
    F(i-1,4,j,k-1) = F(i-1,4,j,k-1) - tyz
    F(i-1,5,j,k-1) = F(i-1,5,j,k-1) - (utyx + vtyy + wtyz + kTy + Hsgs)
    endif
  end subroutine calc_Fv_LES4_in

  !> CUDA Fortran kernel for 4th-order viscous flux with LES SGS model in z direction (interior only)
  !> No boundary condition handling — all threads execute interior 4th-order stencil
  attributes(global) subroutine calc_Gv_LES4_in(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, G)
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
    real(8), intent(inout), device, contiguous :: G(nx-2,5,ny-2,nz-1) !< viscous + SGS flux in z direction
    real(8), shared ::  u(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)
    real(8), shared ::  v(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)
    real(8), shared ::  w(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)
    real(8), shared :: wx(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)
    real(8), shared :: wy(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)
    real(8), shared :: ux(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)
    real(8), shared :: vy(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)
    integer i, j, k, it, jt, kt, kk, k_base
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, Hsgs
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    do kk = kt-2, threadsGv%z+3, blockDim%z
      k = k_base + kk
      if (i <= nx .and. j <= ny .and. 1 <= k .and. k <= nz) then
        u(kk,jt,it) = Q(i,2,j,k)
        v(kk,jt,it) = Q(i,3,j,k)
        w(kk,jt,it) = Q(i,4,j,k)
      endif
      if (3 <= i .and. i <= nx-2 .and. j <= ny .and. 1 <= k .and. k <= nz) then
        ux(kk,jt,it) = (two_third * (-Q(i-1,2,j,k) + Q(i+1,2,j,k)) - one_twelfth * (-Q(i-2,2,j,k) + Q(i+2,2,j,k))) * dx(i)
        wx(kk,jt,it) = (two_third * (-Q(i-1,4,j,k) + Q(i+1,4,j,k)) - one_twelfth * (-Q(i-2,4,j,k) + Q(i+2,4,j,k))) * dx(i)
      endif
      if (i <= nx .and. 3 <= j .and. j <= ny-2 .and. 1 <= k .and. k <= nz) then
        vy(kk,jt,it) = (two_third * (-Q(i,3,j-1,k) + Q(i,3,j+1,k)) - one_twelfth * (-Q(i,3,j-2,k) + Q(i,3,j+2,k))) * dy(j)
        wy(kk,jt,it) = (two_third * (-Q(i,4,j-1,k) + Q(i,4,j+1,k)) - one_twelfth * (-Q(i,4,j-2,k) + Q(i,4,j+2,k))) * dy(j)
      endif
    enddo
    call syncthreads()
    k  = (blockIdx%z-1)*blockDim%z + kt
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
        call calc_tau_straight_LES(mu3, mut3, w(kt-2:kt+3,jt,it), ux(kt-2:kt+3,jt,it), vy(kt-2:kt+3,jt,it), dz(k), tzz, wtzz)
        call calc_tau_cross_LES(mu3, mut3, u(kt-2:kt+3,jt,it), wx(kt-2:kt+3,jt,it), dz(k), tzx, utzx)
        call calc_tau_cross_LES(mu3, mut3, v(kt-2:kt+3,jt,it), wy(kt-2:kt+3,jt,it), dz(k), tzy, vtzy)
        block
          real(8) :: H(4)
          H(:) = Cp * T(i,j,k-1:k+2) + 0.5d0 * (u(kt-1:kt+2,jt,it)**2 + v(kt-1:kt+2,jt,it)**2 + w(kt-1:kt+2,jt,it)**2) + qc2(i,j,k-1:k+2)
          Hsgs = -flux4(mut3) * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) * one_third) * dz(k) / Prt
        end block
      end block
    G(i-1,2,j-1,k) = G(i-1,2,j-1,k) - tzx
    G(i-1,3,j-1,k) = G(i-1,3,j-1,k) - tzy
    G(i-1,4,j-1,k) = G(i-1,4,j-1,k) - tzz
    G(i-1,5,j-1,k) = G(i-1,5,j-1,k) - (utzx + vtzy + wtzz + kTz + Hsgs)
    endif
  end subroutine calc_Gv_LES4_in
end module calc_visc4_les_internal
