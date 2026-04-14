!> Module containing 2nd-order viscous flux computation kernels
!> Computes viscous stresses and heat flux for compressible Navier-Stokes equations
module calc_visc2
  use mod_globals, only : gamma, R, Pr, Prt, dt, threadsEv, threadsFv, threadsGv
  use mod_constant, only : Cp, gamma_1, Cp_over_Pr, one_third, two_third
  implicit none
contains

  !> CUDA Fortran kernel for 2nd-order viscous flux in x direction
  !> Computes stress tensor components and heat flux at cell faces
  !> Accounts for molecular viscosity and thermal conductivity
  attributes(global) subroutine calc_Ev2(nx, ny, nz, dx, dy, dz, Q, T, mu, E)
    integer, intent(in), value                 :: nx                  !< number of grid points in x direction
    integer, intent(in), value                 :: ny                  !< number of grid points in y direction
    integer, intent(in), value                 :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous    :: dx(nx-1)            !< inverse grid spacing in x (1/dx)
    real(8), intent(in), device, contiguous    :: dy(ny-1)            !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous    :: dz(nz-1)            !< inverse grid spacing in z (1/dz)
    real(8), intent(in), device, contiguous    :: Q(5,nx,ny,nz)       !< conservative variables: rho, rho*u, rho*v, rho*w, E
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)         !< temperature at grid points
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)        !< molecular viscosity coefficient
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2) !< viscous flux components in x direction
    real(8), shared :: u(threadsEv%x+1,0:threadsEv%y+1,0:threadsEv%z+1)
    real(8), shared :: v(threadsEv%x+1,0:threadsEv%y+1,threadsEv%z)
    real(8), shared :: w(threadsEv%x+1,0:threadsEv%z+1,threadsEv%y)
    integer i, j, k, it, jt, kt
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx
    real(8) m1, m2, mx, mux, mvx, mwx, muy, mvy, muz, mwz
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    u(it,jt-1:jt+1,kt-1:kt+1) = Q(2,i,j-1:j+1,k-1:k+1)
    v(it,jt-1:jt+1,kt)        = Q(3,i,j-1:j+1,k)
    w(it,kt-1:kt+1,jt)        = Q(4,i,j,k-1:k+1)
    if (it == blockDim%x) then
      u(it+1,jt-1:jt+1,kt-1:kt+1) = Q(2,i+1,j-1:j+1,k-1:k+1)
      v(it+1,jt-1:jt+1,kt)        = Q(3,i+1,j-1:j+1,k)
      w(it+1,kt-1:kt+1,jt)        = Q(4,i+1,j,k-1:k+1)
    endif 
    call syncthreads()

    ! Compute viscous stress tensor components using 2nd-order finite differences
    ! Algorithm: tau_ij = mu * (du_i/dx_j + du_j/dx_i) - (2/3)*mu*delta_ij*(div u)
    ! Viscosity is averaged at cell faces, velocities at grid points in local shared memory
    mx  = 0.5d0 * (mu(i,j,k) + mu(i+1,j,k))        ! Face-centered viscosity in x
    
    ! Heat flux from Fourier's law: q_x = -k * dT/dx where k = density * Cp / Pr * mu
    kTx = Cp_over_Pr * mx * (-T(i,j,k) + T(i+1,j,k)) * dx(i)
    
    ! Compute du/dy and dv/dy with viscosity weighted at cell edges (y-direction)
    ! Using centered differences: du/dy ≈ 0.5 * (u(j+1) - u(j-1)) * dy
    ! Two viscosities at y-edges averaged for interpolation to face centers
    block
      real(8) my1, my2
      my1 = 0.25d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i+1,j-1,k) + mu(i+1,j,  k))
      my2 = 0.25d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i+1,j,  k) + mu(i+1,j+1,k))
      ! du/dy with linear interpolation of mu between two y-faces
      muy = 0.25d0 * (my1 * (-u(it,jt-1,kt) - u(it+1,jt-1,kt)) + (my1 - my2) * (u(it,jt,kt) + u(it+1,jt,kt)) &
                    + my2 * ( u(it,jt+1,kt) + u(it+1,jt+1,kt))) * dy(j)
      mvy = 0.25d0 * (my1 * (-v(it,jt-1,kt) - v(it+1,jt-1,kt)) + (my1 - my2) * (v(it,jt,kt) + v(it+1,jt,kt)) &
                    + my2 * ( v(it,jt+1,kt) + v(it+1,jt+1,kt))) * dy(j)
    end block
    block
      real(8) mz1, mz2
      mz1 = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i+1,j,k-1) + mu(i+1,j,k  ))
      mz2 = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i+1,j,k  ) + mu(i+1,j,k+1))
      muz = 0.25d0 * (mz1 * (-u(it,jt,kt-1) - u(it+1,jt,kt-1)) + (mz1 - mz2) * (u(it,jt,kt) + u(it+1,jt,kt)) &
                    + mz2 * ( u(it,jt,kt+1) + u(it+1,jt,kt+1))) * dz(k)
      mwz = 0.25d0 * (mz1 * (-w(it,kt-1,jt) - w(it+1,kt-1,jt)) + (mz1 - mz2) * (w(it,kt,jt) + w(it+1,kt,jt)) &
                    + mz2 * ( w(it,kt+1,jt) + w(it+1,kt+1,jt))) * dz(k)
    end block
    ! Compute individual strain rate components: dv_j/dx_i
    mux = mx * (-u(it,jt,kt) + u(it+1,jt,kt)) * dx(i)  ! du/dx
    mvx = mx * (-v(it,jt,kt) + v(it+1,jt,kt)) * dx(i)  ! dv/dx
    mwx = mx * (-w(it,kt,jt) + w(it+1,kt,jt)) * dx(i)  ! dw/dx (note: from different shared array index)
    
    ! Assemble stress tensor: tau_ij = mu * S_ij where S is strain rate
    ! Diagonal: tau_ii = (2/3)*mu*(2*du_i/dx_i - du_j/dx_j - du_k/dx_k)
    ! Off-diagonal: tau_ij = mu*(du_i/dx_j + du_j/dx_i)
    txx = two_third * (2.d0 * mux - mvy - mwz)  ! tau_xx with bulk viscosity correction
    txy = muy + mvx                              ! tau_xy = tau_yx (symmetric)
    txz = mwx + muz                              ! tau_xz = tau_zx (symmetric)
    
    ! Compute work terms u*tau for energy equation: (u_i * tau_ij)
    utxx = 0.5d0 * (u(it,jt,kt) + u(it+1,jt,kt)) * txx  ! Average u at cell face × stress
    vtxy = 0.5d0 * (v(it,jt,kt) + v(it+1,jt,kt)) * txy
    wtxz = 0.5d0 * (w(it,kt,jt) + w(it+1,kt,jt)) * txz
    
    ! Accumulate viscous flux components (note subtraction: fluxes point outward)
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx            ! Momentum: -tau_xx in x-flux
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy            ! Momentum: -tau_xy in x-flux  
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz            ! Momentum: -tau_xz in x-flux
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx)  ! Energy: -(u*tau + q)
  end subroutine calc_Ev2
 

  !> CUDA Fortran kernel for 2nd-order viscous flux with LES SGS model in x direction
  !> Computes molecular + subgrid-scale viscous stresses and heat flux
  attributes(global) subroutine calc_Ev_LES2(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, E)
    integer, intent(in), value                 :: nx                  !< number of grid points in x direction
    integer, intent(in), value                 :: ny                  !< number of grid points in y direction
    integer, intent(in), value                 :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous    :: dx(nx-1)            !< inverse grid spacing in x (1/dx)
    real(8), intent(in), device, contiguous    :: dy(ny-1)            !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous    :: dz(nz-1)            !< inverse grid spacing in z (1/dz)
    real(8), intent(in), device, contiguous    :: Q(5,nx,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)         !< temperature at grid points
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)        !< molecular viscosity coefficient
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)       !< turbulent eddy viscosity (LES model)
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)       !< quadratic constitutive relation correction
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2) !< viscous + SGS flux in x direction
    integer i, j, k
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx, Hsgs
    real(8), dimension(2) :: my, mysgs, mz, mzsgs
    real(8) mx, mxsgs, mux, muxsgs, mvx, mvxsgs, mwx, mwxsgs, muy, muysgs, mvy, mvysgs, muz, muzsgs, mwz, mwzsgs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
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
      muy    = 0.25d0 * (my(1)    * (-Q(2,i,j-1,k) + Q(2,i,j,k) - Q(2,i+1,j-1,k) + Q(2,i+1,j,k)) &
                       + my(2)    * (-Q(2,i,j,k) + Q(2,i,j+1,k) - Q(2,i+1,j,k) + Q(2,i+1,j+1,k))) * dy(j)
      muysgs = 0.25d0 * (mysgs(1) * (-Q(2,i,j-1,k) + Q(2,i,j,k) - Q(2,i+1,j-1,k) + Q(2,i+1,j,k)) &
                       + mysgs(2) * (-Q(2,i,j,k) + Q(2,i,j+1,k) - Q(2,i+1,j,k) + Q(2,i+1,j+1,k))) * dy(j)
      mvy    = 0.25d0 * (my(1)    * (-Q(3,i,j-1,k) + Q(3,i,j,k) - Q(3,i+1,j-1,k) + Q(3,i+1,j,k)) &
                       + my(2)    * (-Q(3,i,j,k) + Q(3,i,j+1,k) - Q(3,i+1,j,k) + Q(3,i+1,j+1,k))) * dy(j)
      mvysgs = 0.25d0 * (mysgs(1) * (-Q(3,i,j-1,k) + Q(3,i,j,k) - Q(3,i+1,j-1,k) + Q(3,i+1,j,k)) &
                       + mysgs(2) * (-Q(3,i,j,k) + Q(3,i,j+1,k) - Q(3,i+1,j,k) + Q(3,i+1,j+1,k))) * dy(j)
    end block
    block
      real(8) :: mz(2)
      mz(1)  = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i+1,j,k-1) + mu(i+1,j,k  ))
      mz(2)  = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i+1,j,k  ) + mu(i+1,j,k+1))
      muz    = 0.25d0 * (mz(1)    * (-Q(2,i,j,k-1) + Q(2,i,j,k) - Q(2,i+1,j,k-1) + Q(2,i+1,j,k)) &
                       + mz(2)    * (-Q(2,i,j,k) + Q(2,i,j,k+1) - Q(2,i+1,j,k) + Q(2,i+1,j,k+1))) * dz(k)
      muzsgs = 0.25d0 * (mzsgs(1) * (-Q(2,i,j,k-1) + Q(2,i,j,k) - Q(2,i+1,j,k-1) + Q(2,i+1,j,k)) &
                       + mzsgs(2) * (-Q(2,i,j,k) + Q(2,i,j,k+1) - Q(2,i+1,j,k) + Q(2,i+1,j,k+1))) * dz(k)
      mwz    = 0.25d0 * (mz(1)    * (-Q(4,i,j,k-1) + Q(4,i,j,k) - Q(4,i+1,j,k-1) + Q(4,i+1,j,k)) &
                       + mz(2)    * (-Q(4,i,j,k) + Q(4,i,j,k+1) - Q(4,i+1,j,k) + Q(4,i+1,j,k+1))) * dz(k)
      mwzsgs = 0.25d0 * (mzsgs(1) * (-Q(4,i,j,k-1) + Q(4,i,j,k) - Q(4,i+1,j,k-1) + Q(4,i+1,j,k)) &
                       + mzsgs(2) * (-Q(4,i,j,k) + Q(4,i,j,k+1) - Q(4,i+1,j,k) + Q(4,i+1,j,k+1))) * dz(k)
    end block
    mux    = mx    * (-Q(2,i,j,k) + Q(2,i+1,j,k)) * dx(i)
    muxsgs = mxsgs * (-Q(2,i,j,k) + Q(2,i+1,j,k)) * dx(i)
    mvx    = mx    * (-Q(3,i,j,k) + Q(3,i+1,j,k)) * dx(i)
    mvxsgs = mxsgs * (-Q(3,i,j,k) + Q(3,i+1,j,k)) * dx(i)
    mwx    = mx    * (-Q(4,i,j,k) + Q(4,i+1,j,k)) * dx(i)
    mwxsgs = mxsgs * (-Q(4,i,j,k) + Q(4,i+1,j,k)) * dx(i)
    txx    = two_third * (2.d0 * mux - mvy - mwz)
    txy    = muy + mvx
    txz    = mwx + muz
    utxx   = 0.5d0 * (Q(2,i,j,k) + Q(2,i+1,j,k)) * txx
    vtxy   = 0.5d0 * (Q(3,i,j,k) + Q(3,i+1,j,k)) * txy
    wtxz   = 0.5d0 * (Q(4,i,j,k) + Q(4,i+1,j,k)) * txz
    txx    = txx + two_third * (2.d0 * muxsgs - mvysgs - mwzsgs)
    txy    = txy + muysgs + mvxsgs
    txz    = txz + mwxsgs + muzsgs
    block
      real(8) :: H(2)
      H(:) = Cp * T(i:i+1,j,k) + 0.5d0 * (Q(2,i:i+1,j,k)**2 + Q(3,i:i+1,j,k)**2 + Q(4,i:i+1,j,k)**2) + qc2(i:i+1,j,k)
      Hsgs = -mxsgs * (-H(1) + H(2)) * dx(i) / Prt
    end block
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx + Hsgs)
  end subroutine calc_Ev_LES2
 

  !> CUDA Fortran kernel for 2nd-order viscous flux in y direction
  !> Computes stress tensor components and heat flux at cell faces
  attributes(global) subroutine calc_Fv2(nx, ny, nz, dy, dx, dz, Q, T, mu, F)
    integer, intent(in), value                 :: nx                  !< number of grid points in x direction
    integer, intent(in), value                 :: ny                  !< number of grid points in y direction
    integer, intent(in), value                 :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous    :: dy(ny-1)            !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous    :: dx(nx-1)            !< inverse grid spacing in x (1/dx)
    real(8), intent(in), device, contiguous    :: dz(nz-1)            !< inverse grid spacing in z (1/dz)
    real(8), intent(in), device, contiguous    :: Q(5,nx,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)         !< temperature at grid points
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)        !< molecular viscosity coefficient
    real(8), intent(inout), device, contiguous :: F(5,nx-2,ny-1,nz-2) !< viscous flux components in y direction
    real(8), shared :: u(threadsFv%y+1,0:threadsFv%x+1,threadsFv%z)
    real(8), shared :: v(threadsFv%y+1,0:threadsFv%x+1,0:threadsFv%z+1)
    real(8), shared :: w(threadsFv%y+1,0:threadsFv%z+1,threadsFv%x)
    integer i, j, k, it, jt, kt
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy
    real(8) m1, m2, my, muy, mvy, mwy, mvz, mwz, mux, mvx
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    u(jt,it-1:it+1,kt)        = Q(2,i-1:i+1,j,k)
    v(jt,it-1:it+1,kt-1:kt+1) = Q(3,i-1:i+1,j,k-1:k+1)
    w(jt,kt-1:kt+1,it)        = Q(4,i,j,k-1:k+1)
    if (jt == blockDim%y) then
      u(jt+1,it-1:it+1,kt)        = Q(2,i-1:i+1,j+1,k)
      v(jt+1,it-1:it+1,kt-1:kt+1) = Q(3,i-1:i+1,j+1,k-1:k+1)
      w(jt+1,kt-1:kt+1,it)        = Q(4,i,j+1,k-1:k+1)
    endif 
    call syncthreads()

    block
      real(8) mx1, mx2
      mx1 = 0.25d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j+1,k) + mu(i,  j+1,k))
      mx2 = 0.25d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j+1,k) + mu(i+1,j+1,k))
      mux = 0.25d0 * (mx1 * (-u(jt,it-1,kt) - u(jt+1,it-1,kt)) + (mx1 - mx2) * (u(jt,it,kt) + u(jt+1,it,kt)) &
                    + mx2 * ( u(jt,it+1,kt) + u(jt+1,it+1,kt))) * dx(i)
      mvx = 0.25d0 * (mx1 * (-v(jt,it-1,kt) - v(jt+1,it-1,kt)) + (mx1 - mx2) * (v(jt,it,kt) + v(jt+1,it,kt)) &
                    + mx2 * ( v(jt,it+1,kt) + v(jt+1,it+1,kt))) * dx(i)
    end block
    my  = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k))
    kTy = Cp_over_Pr * my * (-T(i,j,k) + T(i,j+1,k)) * dy(j)
    block
      real(8) mz1, mz2
      mz1 = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i,j+1,k-1) + mu(i,j+1,k  ))
      mz2 = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i,j+1,k  ) + mu(i,j+1,k+1))
      mvz = 0.25d0 * (mz1 * (-v(jt,it,kt-1) - v(jt+1,it,kt-1)) + (mz1 - mz2) * (v(jt,it,kt) + v(jt+1,it,kt)) &
                    + mz2 * ( v(jt,it,kt+1) + v(jt+1,it,kt+1))) * dz(k)
      mwz = 0.25d0 * (mz1 * (-w(jt,kt-1,it) - w(jt+1,kt-1,it)) + (mz1 - mz2) * (w(jt,kt,it) + w(jt+1,kt,it)) &
                    + mz2 * ( w(jt,kt+1,it) + w(jt+1,kt+1,it))) * dz(k)
    end block
    muy = my * (-u(jt,it,kt) + u(jt+1,it,kt)) * dy(j)
    mvy = my * (-v(jt,it,kt) + v(jt+1,it,kt)) * dy(j)
    mwy = my * (-w(jt,kt,it) + w(jt+1,kt,it)) * dy(j)
    tyx = muy + mvx
    tyy = two_third * (2.d0 * mvy - mwz - mux)
    tyz = mvz + mwy
    utyx = 0.5d0 * (u(jt,it,kt) + u(jt+1,it,kt)) * tyx
    vtyy = 0.5d0 * (v(jt,it,kt) + v(jt+1,it,kt)) * tyy
    wtyz = 0.5d0 * (w(jt,kt,it) + w(jt+1,kt,it)) * tyz
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy)
  end subroutine calc_Fv2
 

  !> CUDA Fortran kernel for 2nd-order viscous flux with LES SGS model in y direction
  !> Computes molecular + subgrid-scale viscous stresses and heat flux
  attributes(global) subroutine calc_Fv_LES2(nx, ny, nz, dy, dx, dz, Q, T, mu, mut, qc2, F)
    integer, intent(in), value                 :: nx                  !< number of grid points in x direction
    integer, intent(in), value                 :: ny                  !< number of grid points in y direction
    integer, intent(in), value                 :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous    :: dy(ny-1)            !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous    :: dx(nx-1)            !< inverse grid spacing in x (1/dx)
    real(8), intent(in), device, contiguous    :: dz(nz-1)            !< inverse grid spacing in z (1/dz)
    real(8), intent(in), device, contiguous    :: Q(5,nx,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)         !< temperature at grid points
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)        !< molecular viscosity coefficient
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)       !< turbulent eddy viscosity (LES model)
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)       !< quadratic constitutive relation correction
    real(8), intent(inout), device, contiguous :: F(5,nx-2,ny-1,nz-2) !< viscous + SGS flux in y direction
    integer i, j, k
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, Hsgs
    real(8), dimension(2) :: u2, v2, w2, mz, mzsgs, mx, mxsgs
    real(8) my, mysgs, muy, muysgs, mvy, mvysgs, mwy, mwysgs, mvz, mvzsgs, mwz, mwzsgs, mux, muxsgs, mvx, mvxsgs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
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
      mux    = 0.25d0 * (mx(1)    * (-Q(2,i-1,j,k) + Q(2,i,j,k) - Q(2,i-1,j+1,k) + Q(2,i,j+1,k)) &
                       + mx(2)    * (-Q(2,i,j,k) + Q(2,i+1,j,k) - Q(2,i,j+1,k) + Q(2,i+1,j+1,k))) * dx(i)
      muxsgs = 0.25d0 * (mxsgs(1) * (-Q(2,i-1,j,k) + Q(2,i,j,k) - Q(2,i-1,j+1,k) + Q(2,i,j+1,k)) &
                       + mxsgs(2) * (-Q(2,i,j,k) + Q(2,i+1,j,k) - Q(2,i,j+1,k) + Q(2,i+1,j+1,k))) * dx(i)
      mvx    = 0.25d0 * (mx(1)    * (-Q(3,i-1,j,k) + Q(3,i,j,k) - Q(3,i-1,j+1,k) + Q(3,i,j+1,k)) &
                       + mx(2)    * (-Q(3,i,j,k) + Q(3,i+1,j,k) - Q(3,i,j+1,k) + Q(3,i+1,j+1,k))) * dx(i)
      mvxsgs = 0.25d0 * (mxsgs(1) * (-Q(3,i-1,j,k) + Q(3,i,j,k) - Q(3,i-1,j+1,k) + Q(3,i,j+1,k)) &
                       + mxsgs(2) * (-Q(3,i,j,k) + Q(3,i+1,j,k) - Q(3,i,j+1,k) + Q(3,i+1,j+1,k))) * dx(i)
    end block
    my  = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k))
    kTy = Cp_over_Pr * my * (-T(i,j,k) + T(i,j+1,k)) * dy(j)
    block
      real(8) :: mz(2)
      mz(1)  = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i,j+1,k-1) + mu(i,j+1,k  ))
      mz(2)  = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i,j+1,k  ) + mu(i,j+1,k+1))
      mvz    = 0.25d0 * (mz(1)    * (-Q(3,i,j,k-1) + Q(3,i,j,k) - Q(3,i,j+1,k-1) + Q(3,i,j+1,k)) &
                       + mz(2)    * (-Q(3,i,j,k) + Q(3,i,j,k+1) - Q(3,i,j+1,k) + Q(3,i,j+1,k+1))) * dz(k)
      mvzsgs = 0.25d0 * (mzsgs(1) * (-Q(3,i,j,k-1) + Q(3,i,j,k) - Q(3,i,j+1,k-1) + Q(3,i,j+1,k)) &
                       + mzsgs(2) * (-Q(3,i,j,k) + Q(3,i,j,k+1) - Q(3,i,j+1,k) + Q(3,i,j+1,k+1))) * dz(k)
      mwz    = 0.25d0 * (mz(1)    * (-Q(4,i,j,k-1) + Q(4,i,j,k) - Q(4,i,j+1,k-1) + Q(4,i,j+1,k)) &
                       + mz(2)    * (-Q(4,i,j,k) + Q(4,i,j,k+1) - Q(4,i,j+1,k) + Q(4,i,j+1,k+1))) * dz(k)
      mwzsgs = 0.25d0 * (mzsgs(1) * (-Q(4,i,j,k-1) + Q(4,i,j,k) - Q(4,i,j+1,k-1) + Q(4,i,j+1,k)) &
                       + mzsgs(2) * (-Q(4,i,j,k) + Q(4,i,j,k+1) - Q(4,i,j+1,k) + Q(4,i,j+1,k+1))) * dz(k)
    end block
    muy    = my    * (-Q(2,i,j,k) + Q(2,i,j+1,k)) * dy(j)
    muysgs = mysgs * (-Q(2,i,j,k) + Q(2,i,j+1,k)) * dy(j)
    mvy    = my    * (-Q(3,i,j,k) + Q(3,i,j+1,k)) * dy(j)
    mvysgs = mysgs * (-Q(3,i,j,k) + Q(3,i,j+1,k)) * dy(j)
    mwy    = my    * (-Q(4,i,j,k) + Q(4,i,j+1,k)) * dy(j)
    mwysgs = mysgs * (-Q(4,i,j,k) + Q(4,i,j+1,k)) * dy(j)
    tyx    = muy + mvx
    tyy    = two_third * (2.d0 * mvy - mwz - mux)
    tyz    = mvz + mwy
    utyx   = 0.5d0 * (Q(2,i,j,k) + Q(2,i,j+1,k)) * tyx
    vtyy   = 0.5d0 * (Q(3,i,j,k) + Q(3,i,j+1,k)) * tyy
    wtyz   = 0.5d0 * (Q(4,i,j,k) + Q(4,i,j+1,k)) * tyz
    tyx    = tyx + muysgs + mvxsgs
    tyy    = tyy + two_third * (2.d0 * mvysgs - mwzsgs - muxsgs)
    tyz    = tyz + mvzsgs + mwysgs
    block
      real(8) :: H(2)
      H(:) = Cp * T(i,j:j+1,k) + 0.5d0 * (Q(2,i,j:j+1,k)**2 + Q(3,i,j:j+1,k)**2 + Q(4,i,j:j+1,k)**2) + qc2(i,j:j+1,k)
      Hsgs = -mysgs * (-H(1) + H(2)) * dy(j) / Prt
    end block
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy + Hsgs)
  end subroutine calc_Fv_LES2
 

  !> CUDA Fortran kernel for 2nd-order viscous flux in z direction
  !> Computes stress tensor components and heat flux at cell faces
  attributes(global) subroutine calc_Gv2(nx, ny, nz, dx, dy, dz, Q, T, mu, G)
    integer, intent(in), value                 :: nx                  !< number of grid points in x direction
    integer, intent(in), value                 :: ny                  !< number of grid points in y direction
    integer, intent(in), value                 :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous    :: dx(nx-1)            !< inverse grid spacing in x (1/dx)
    real(8), intent(in), device, contiguous    :: dy(ny-1)            !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous    :: dz(nz-1)            !< inverse grid spacing in z (1/dz)
    real(8), intent(in), device, contiguous    :: Q(5,nx,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)         !< temperature at grid points
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)        !< molecular viscosity coefficient
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1) !< viscous flux components in z direction
    real(8), shared :: u(threadsGv%z+1,0:threadsGv%x+1,threadsGv%y)
    real(8), shared :: v(threadsGv%z+1,0:threadsGv%y+1,threadsGv%x)
    real(8), shared :: w(threadsGv%z+1,0:threadsGv%x+1,0:threadsGv%y+1)
    integer i, j, k, it, jt, kt
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz
    real(8) m1, m2, mz, muz, mvz, mwz, mwx, mux, mvy, mwy
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    u(kt,it-1:it+1,jt)        = Q(2,i-1:i+1,j,k)
    v(kt,jt-1:jt+1,it)        = Q(3,i,j-1:j+1,k)
    w(kt,it-1:it+1,jt-1:jt+1) = Q(4,i-1:i+1,j-1:j+1,k)
    if (kt == blockDim%z) then
      u(kt+1,it-1:it+1,jt)        = Q(2,i-1:i+1,j,k+1)
      v(kt+1,jt-1:jt+1,it)        = Q(3,i,j-1:j+1,k+1)
      w(kt+1,it-1:it+1,jt-1:jt+1) = Q(4,i-1:i+1,j-1:j+1,k+1)
    endif 
    call syncthreads()

    block
      real(8) mx1, mx2
      mx1 = 0.25d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j,k+1) + mu(i,  j,k+1))
      mx2 = 0.25d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j,k+1) + mu(i+1,j,k+1))
      mux = 0.25d0 * (mx1 * (-u(kt,it-1,jt) - u(kt+1,it-1,jt)) + (mx1 - mx2) * (u(kt,it,jt) + u(kt+1,it,jt)) &
                    + mx2 * ( u(kt,it+1,jt) + u(kt+1,it+1,jt))) * dx(i)
      mwx = 0.25d0 * (mx1 * (-w(kt,it-1,jt) - w(kt+1,it-1,jt)) + (mx1 - mx2) * (w(kt,it,jt) + w(kt+1,it,jt)) &
                    + mx2 * ( w(kt,it+1,jt) + w(kt+1,it+1,jt))) * dx(i)
    end block
    block
      real(8) my1, my2
      my1 = 0.25d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i,j-1,k+1) + mu(i,j,  k+1))
      my2 = 0.25d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i,j,  k+1) + mu(i,j+1,k+1))
      mvy = 0.25d0 * (my1 * (-v(kt,jt-1,it) - v(kt+1,jt-1,it)) + (my1 - my2) * (v(kt,jt,it) + v(kt+1,jt,it)) &
                    + my2 * ( v(kt,jt+1,it) + v(kt+1,jt+1,it))) * dy(j)
      mwy = 0.25d0 * (my1 * (-w(kt,it,jt-1) - w(kt+1,it,jt-1)) + (my1 - my2) * (w(kt,it,jt) + w(kt+1,it,jt)) &
                    + my2 * ( w(kt,it,jt+1) + w(kt+1,it,jt+1))) * dy(j)
    end block
    mz  = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1))
    kTz = Cp_over_Pr * mz * (-T(i,j,k) + T(i,j,k+1)) * dz(k)
    muz = mz * (-u(kt,it,jt) + u(kt+1,it,jt)) * dz(k)
    mvz = mz * (-v(kt,jt,it) + v(kt+1,jt,it)) * dz(k)
    mwz = mz * (-w(kt,it,jt) + w(kt+1,it,jt)) * dz(k)
    tzx = mwx + muz
    tzy = mvz + mwy
    tzz = two_third * (2.d0 * mwz - mux - mvy)
    utzx = 0.5d0 * (u(kt,it,jt) + u(kt+1,it,jt)) * tzx
    vtzy = 0.5d0 * (v(kt,jt,it) + v(kt+1,jt,it)) * tzy
    wtzz = 0.5d0 * (w(kt,it,jt) + w(kt+1,it,jt)) * tzz
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz)
  end subroutine calc_Gv2


  attributes(global) subroutine calc_Gv_LES2(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, G)
    integer, intent(in), value                 :: nx                  !< number of grid points in x direction
    integer, intent(in), value                 :: ny                  !< number of grid points in y direction
    integer, intent(in), value                 :: nz                  !< number of grid points in z direction
    real(8), intent(in), device, contiguous    :: dx(nx-1)            !< inverse grid spacing in x (1/dx)
    real(8), intent(in), device, contiguous    :: dy(ny-1)            !< inverse grid spacing in y (1/dy)
    real(8), intent(in), device, contiguous    :: dz(nz-1)            !< inverse grid spacing in z (1/dz)
    real(8), intent(in), device, contiguous    :: Q(5,nx,ny,nz)       !< conservative variables
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)         !< temperature at grid points
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)        !< molecular viscosity coefficient
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)       !< turbulent eddy viscosity (LES model)
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)       !< quadratic constitutive relation correction
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1) !< viscous + SGS flux in z direction
    integer i, j, k
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, Hsgs
    real(8), dimension(2) :: mx, mxsgs, my, mysgs
    real(8) mz, mzsgs, muz, muzsgs, mvz, mvzsgs, mwz, mwzsgs, mwx, mwxsgs, mux, muxsgs, mvy, mvysgs, mwy, mwysgs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
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
      mux    = 0.25d0 * (mx(1)    * (-Q(2,i-1,j,k) + Q(2,i,j,k) - Q(2,i-1,j,k+1) + Q(2,i,j,k+1)) &
                       + mx(2)    * (-Q(2,i,j,k) + Q(2,i+1,j,k) - Q(2,i,j,k+1) + Q(2,i+1,j,k+1))) * dx(i)
      muxsgs = 0.25d0 * (mxsgs(1) * (-Q(2,i-1,j,k) + Q(2,i,j,k) - Q(2,i-1,j,k+1) + Q(2,i,j,k+1)) &
                       + mxsgs(2) * (-Q(2,i,j,k) + Q(2,i+1,j,k) - Q(2,i,j,k+1) + Q(2,i+1,j,k+1))) * dx(i)
      mwx    = 0.25d0 * (mx(1)    * (-Q(4,i-1,j,k) + Q(4,i,j,k) - Q(4,i-1,j,k+1) + Q(4,i,j,k+1)) &
                       + mx(2)    * (-Q(4,i,j,k) + Q(4,i+1,j,k) - Q(4,i,j,k+1) + Q(4,i+1,j,k+1))) * dx(i)
      mwxsgs = 0.25d0 * (mxsgs(1) * (-Q(4,i-1,j,k) + Q(4,i,j,k) - Q(4,i-1,j,k+1) + Q(4,i,j,k+1)) &
                       + mxsgs(2) * (-Q(4,i,j,k) + Q(4,i+1,j,k) - Q(4,i,j,k+1) + Q(4,i+1,j,k+1))) * dx(i)
    end block
    block
      real(8) :: my(2)
      my(1)  = 0.25d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i,j-1,k+1) + mu(i,j,  k+1))
      my(2)  = 0.25d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i,j,  k+1) + mu(i,j+1,k+1))
      mvy    = 0.25d0 * (my(1)    * (-Q(3,i,j-1,k) + Q(3,i,j,k) - Q(3,i,j-1,k+1) + Q(3,i,j,k+1)) &
                       + my(2)    * (-Q(3,i,j,k) + Q(3,i,j+1,k) - Q(3,i,j,k+1) + Q(3,i,j+1,k+1))) * dy(j)
      mvysgs = 0.25d0 * (mysgs(1) * (-Q(3,i,j-1,k) + Q(3,i,j,k) - Q(3,i,j-1,k+1) + Q(3,i,j,k+1)) &
                       + mysgs(2) * (-Q(3,i,j,k) + Q(3,i,j+1,k) - Q(3,i,j,k+1) + Q(3,i,j+1,k+1))) * dy(j)
      mwy    = 0.25d0 * (my(1)    * (-Q(4,i,j-1,k) + Q(4,i,j,k) - Q(4,i,j-1,k+1) + Q(4,i,j,k+1)) &
                       + my(2)    * (-Q(4,i,j,k) + Q(4,i,j+1,k) - Q(4,i,j,k+1) + Q(4,i,j+1,k+1))) * dy(j)
      mwysgs = 0.25d0 * (mysgs(1) * (-Q(4,i,j-1,k) + Q(4,i,j,k) - Q(4,i,j-1,k+1) + Q(4,i,j,k+1)) &
                       + mysgs(2) * (-Q(4,i,j,k) + Q(4,i,j+1,k) - Q(4,i,j,k+1) + Q(4,i,j+1,k+1))) * dy(j)
    end block
    mz     = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1))
    kTz    = Cp_over_Pr * mz * (-T(i,j,k) + T(i,j,k+1)) * dz(k)
    muz    = mz    * (-Q(2,i,j,k) + Q(2,i,j,k+1)) * dz(k)
    muzsgs = mzsgs * (-Q(2,i,j,k) + Q(2,i,j,k+1)) * dz(k)
    mvz    = mz    * (-Q(3,i,j,k) + Q(3,i,j,k+1)) * dz(k)
    mvzsgs = mzsgs * (-Q(3,i,j,k) + Q(3,i,j,k+1)) * dz(k)
    mwz    = mz    * (-Q(4,i,j,k) + Q(4,i,j,k+1)) * dz(k)
    mwzsgs = mzsgs * (-Q(4,i,j,k) + Q(4,i,j,k+1)) * dz(k)
    tzx    = mwx + muz
    tzy    = mvz + mwy
    tzz    = two_third * (2.d0 * mwz - mux - mvy)
    utzx   = 0.5d0 * (Q(2,i,j,k) + Q(2,i,j,k+1)) * tzx
    vtzy   = 0.5d0 * (Q(3,i,j,k) + Q(3,i,j,k+1)) * tzy
    wtzz   = 0.5d0 * (Q(4,i,j,k) + Q(4,i,j,k+1)) * tzz
    tzx    = tzx + mwx + muz
    tzy    = tzy + mvz + mwy
    tzz    = tzz + two_third * (2.d0 * mwz - mux - mvy)
    block
      real(8) :: H(2)
      H(:) = Cp * T(i,j,k:k+1) + 0.5d0 * (Q(2,i,j,k:k+1)**2 + Q(3,i,j,k:k+1)**2 + Q(4,i,j,k:k+1)**2) + qc2(i,j,k:k+1)
      Hsgs = -mzsgs * (-H(1) + H(2)) * dz(k) / Prt
    end block
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz + Hsgs)
  end subroutine calc_Gv_LES2
end module calc_visc2

