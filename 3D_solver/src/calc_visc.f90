module calc_visc
  use mod_globals, only : id_visc, id_turbulence, id_av, gamma, R, Pr, Prt, threadsEv, threadsFv, threadsGv
  use mod_constant, only : Cp, gamma_1
  use calc_me4_base
  use calc_visc_visbal
  implicit none
contains
  attributes(global) subroutine calc_Ev(nx, ny, nz, dx, dy, dz, Q, E)
    use calc_sutherland, only : mu6, mu2, mu23
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(inout), device :: E(5,nx-1,ny-2,nz-2)
    integer i, j, k, it, jt, kt
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx
    real(8), shared :: u(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), shared :: v(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), shared :: w(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    ! 4th-order 2
    real(8), device :: T6(6), mu(3)
    ! 2nd-order 2
    real(8) :: T233(2,3,3), Ty(2,3), Tz(2,3), Tx(2), my(2), mz(2), mx
    it = threadIdx%x
    jt = threadIdx%y + 1
    kt = threadIdx%z + 1
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      u(it-2:it+3,jt-2:jt+2,kt-2:kt+2) = Q(2,i-2:i+3,j-2:j+2,k-2:k+2)
      v(it-2:it+3,jt-2:jt+2,kt) = Q(3,i-2:i+3,j-2:j+2,k)
      w(it-2:it+3,jt,kt-2:kt+2) = Q(4,i-2:i+3,j,k-2:k+2)
    endif
    call syncthreads()
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      T6(:) = Q(5,i-2:i+3,j,k) / (R * Q(1,i-2:i+3,j,k))
      mu(:) = mu6(T6(:))
      kTx   = heat_conduction6(mu(:), T6(:), dx(i))
      call calc_derivative_x4(nx, ny, nz, it, jt, kt, dx(i), dy(j), dz(k), mu, u, v, w, txx, txy, txz, utxx, vtxy, wtxz)
    else
      T233(:,:,:) = Q(5,i:i+1,j-1:j+1,k-1:k+1) / (R * Q(1,i:i+1,j-1:j+1,k-1:k+1))
      Tx(:)       = T233(:,2,2)
      Ty(:,:)     = T233(:,:,2)
      Tz(:,:)     = T233(:,2,:)
      mx          = mu2(Tx(:))
      my(:)       = mu23(Ty(:,:))
      mz(:)       = mu23(Tz(:,:))
      kTx         = Cp * mx * (-Tx(1) + Tx(2)) * dx(i) / Pr
      call calc_derivative_x2(nx, ny, nz, it, jt, kt, dx(i), dy(j), dz(k), mx, my, mz, u, v, w, txx, txy, txz, utxx, vtxy, wtxz)
    endif
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx)
  end subroutine calc_Ev
  
  attributes(global) subroutine calc_Ev_LES(nx, ny, nz, dx, dy, dz, Q, mut, qc2, E)
    use calc_sutherland, only : mu6, mu2, mu23
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(in), device    :: mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device :: E(5,nx-1,ny-2,nz-2)
    integer i, j, k, it, jt, kt
    real(8) :: txx_txxsgs, txy_txysgs, txz_txzsgs, utxx, vtxy, wtxz, kTx, mutx, H(4), Hsgs
    real(8), shared :: u(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), shared :: v(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), shared :: w(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    ! 4th-order 2
    real(8), device :: T6(6), mu(3), mt(3)
    ! 2nd-order 2
    real(8), device :: T233(2,3,3), Ty(2,3), Tz(2,3), Tx(2), my(2), mty(2), mz(2), mtz(2)
    real(8) mx, mtx
    it = threadIdx%x
    jt = threadIdx%y + 1
    kt = threadIdx%z + 1
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      u(it-2:it+3,jt-2:jt+2,kt-2:kt+2) = Q(2,i-2:i+3,j-2:j+2,k-2:k+2)
      v(it-2:it+3,jt-2:jt+2,kt) = Q(3,i-2:i+3,j-2:j+2,k)
      w(it-2:it+3,jt,kt-2:kt+2) = Q(4,i-2:i+3,j,k-2:k+2)
    endif
    call syncthreads()
    if (id_visc ==2 .and. 3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      T6(:) = Q(5,i-2:i+3,j,k) / (R * Q(1,i-2:i+3,j,k))
      mu(:) = mu6(T6(:))
      kTx   = heat_conduction6(mu(:), T6(:), dx(i))
      mt(:) = mut(i-1:i+1,j,k)
      H(:)  = (gamma * Q(5,i-1:i+2,j,k) / (Q(1,i-1:i+2,j,k) * gamma_1)) &
              + 0.5d0 * (u(it-1:it+2,jt,kt)**2 + v(it-1:it+2,jt,kt)**2 + w(it-1:it+2,jt,kt)**2) + qc2(i-1:i+2,j,k)
      Hsgs  = -mutx * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) / 3.d0) * dx(i) / Prt
      call calc_derivative_sgs_x4(nx, ny, nz, it, jt, kt, dx(i), dy(j), dz(k), mu, mt, u, v, w, txx_txxsgs, txy_txysgs, txz_txzsgs, utxx, vtxy, wtxz)
    else
      T233(:,:,:) = Q(5,i:i+1,j-1:j+1,k-1:k+1) / (R * Q(1,i:i+1,j-1:j+1,k-1:k+1))
      Tx(:)       = T233(:,2,2)
      Ty(:,:)     = T233(:,:,2)
      Tz(:,:)     = T233(:,2,:)
      mx          = mu2(Tx(:))
      my(:)       = mu23(Ty(:,:))
      mz(:)       = mu23(Tz(:,:))
      kTx         = Cp * mx * (-Tx(1) + Tx(2)) * dx(i) / Pr
      mtx    = 0.5d0 * (mut(i,j,k) + mut(i+1,j,k))
      mty(:) = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i+1,j-1,k) + mut(i+1,j,k)), &
                 0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i+1,j,k) + mut(i+1,j+1,k))/)
      mtz(:) = (/0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i+1,j,k-1) + mut(i+1,j,k)), &
                 0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i+1,j,k) + mut(i+1,j,k+1))/)
      H(2:3) = (gamma * Q(5,i:i+1,j,k) / (Q(1,i:i+1,j,k) * gamma_1)) &
               + 0.5d0 * (Q(2,i:i+1,j,k)**2 + Q(3,i,j:j+1,k)**2 + Q(4,i:i+1,j,k)**2) + qc2(i:i+1,j,k)
      Hsgs   = -mx * (-H(2) + H(3)) * dx(i) / Prt
      call calc_derivative_sgs_x2(nx, ny, nz, it, jt, kt, dx(i), dy(j), dz(k), mx, mtx, my, mty, mz, mtz, u, v, w, txx_txxsgs, txy_txysgs, txz_txzsgs, utxx, vtxy, wtxz)
    endif
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx_txxsgs
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy_txysgs
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz_txzsgs
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx + Hsgs)
  end subroutine calc_Ev_LES
  
  attributes(global) subroutine calc_Fv(nx, ny, nz, dy, dx, dz, Q, F)
    use calc_sutherland, only : mu6, mu2, mu23, mu32
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(inout), device :: F(5,nx-2,ny-1,nz-2)
    integer i, j, k, it, jt, kt
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, muty
    real(8), shared :: u(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), shared :: v(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), shared :: w(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    ! 4th-order 2
    real(8), device :: T6(6), mu(3)
    ! 2nd-order 2
    real(8) :: T323(3,2,3), Tx(3,2), Tz(2,3), Ty(2), mx(2), mz(2), my
    it = threadIdx%x + 1
    jt = threadIdx%y
    kt = threadIdx%z + 1
    i = (blockIdx%x-1)*blockDim%x + it
    j = (blockIdx%y-1)*blockDim%y + jt
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      v(it-2:it+2,jt-2:jt+3,kt-2:kt+2) = Q(3,i-2:i+2,j-2:j+3,k-2:k+2)
      u(it-2:it+2,jt-2:jt+3,kt) = Q(2,i-2:i+2,j-2:j+3,k)
      w(it,jt-2:jt+3,kt-2:kt+2) = Q(4,i,j-2:j+3,k-2:k+2)
    endif
    call syncthreads()
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      T6(:) = Q(5,i,j-2:j+3,k) / (R * Q(1,i,j-2:j+3,k))
      mu(:) = mu6(T6(:))
      kTy   = heat_conduction6(mu(:), T6(:), dy(j))
      call calc_derivative_y4(nx, ny, nz, it, jt, kt, dx(i), dy(j), dz(k), mu, u, v, w, tyy, tyx, tyz, vtyy, utyx, wtyz)
    else
      T323(:,:,:) = Q(5,i-1:i+1,j:j+1,k-1:k+1) / (R * Q(1,i-1:i+1,j:j+1,k-1:k+1))
      Tx(:,:)     = T323(:,:,2)
      Ty(:)       = T323(2,:,2)
      Tz(:,:)     = T323(2,:,:)
      mx(:)       = mu23(Tx(:,:))
      my          = mu2(Ty(:))
      mz(:)       = mu32(Tz(:,:))
      kTy         = Cp * my * (-Ty(1) + Ty(2)) * dy(j) / Pr
      call calc_derivative_y2(nx, ny, nz, it, jt, kt, dx(i), dy(j), dz(k), my, mx, mz, u, v, w, tyy, tyx, tyz, vtyy, utyx, wtyz)
    endif
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy)
  end subroutine calc_Fv
  
  attributes(global) subroutine calc_Fv_LES(nx, ny, nz, dy, dx, dz, Q, mut, qc2, F)
    use calc_sutherland, only : mu6, mu2, mu23, mu32
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(in), device    :: mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device :: F(5,nx-2,ny-1,nz-2)
    integer i, j, k, it, jt, kt
    real(8) :: tyx_tyxsgs, tyy_tyysgs, tyz_tyzsgs, utyx, vtyy, wtyz, kTy, muty, H(4), Hsgs 
    real(8), shared :: u(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), shared :: v(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), shared :: w(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    ! 4th-order 2
    real(8), device :: T6(6), mu(3), mt(3)
    ! 2nd-order 2
    real(8), dimension(3,2,3), device :: T323(3,2,3), Tx(3,2), Tz(2,3), Ty(2), mz(2), mtz(2), mx(2), mtx(2)
    real(8) my, mty
    it = threadIdx%x + 1
    jt = threadIdx%y
    kt = threadIdx%z + 1
    i = (blockIdx%x-1)*blockDim%x + it
    j = (blockIdx%y-1)*blockDim%y + jt
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      v(it-2:it+2,jt-2:jt+3,kt-2:kt+2) = Q(3,i-2:i+2,j-2:j+3,k-2:k+2)
      u(it-2:it+2,jt-2:jt+3,kt) = Q(2,i-2:i+2,j-2:j+3,k)
      w(it,jt-2:jt+3,kt-2:kt+2) = Q(4,i,j-2:j+3,k-2:k+2)
    endif
    call syncthreads()
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      T6(:) = Q(5,i,j-2:j+3,k) / (R * Q(1,i,j-2:j+3,k))
      mu(:) = mu6(T6(:))
      kTy   = heat_conduction6(mu(:), T6(:), dy(j))
      mt(:) = mut(i,j-1:j+1,k)
      H(:)  = (gamma * Q(5,i,j-1:j+2,k) / (Q(1,i,j-1:j+2,k) * gamma_1)) &
              + 0.5d0 * (u(it,jt-1:jt+2,kt)**2 + v(it,jt-1:jt+2,kt)**2 + w(it,jt-1:jt+2,kt)**2) + qc2(i,j-1:j+2,k)
      Hsgs  = -muty * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) / 3.d0) * dy(j) / Prt
      call calc_derivative_sgs_y4(nx, ny, nz, it, jt, kt, dx(i), dy(j), dz(k), mu, mt, u, v, w, tyy_tyysgs, tyx_tyxsgs, tyz_tyzsgs, vtyy, utyx, wtyz)
    else
      T323(:,:,:) = Q(5,i-1:i+1,j:j+1,k-1:k+1) / (R * Q(1,i-1:i+1,j:j+1,k-1:k+1))
      Tx(:,:)     = T323(:,:,2)
      Ty(:)       = T323(2,:,2)
      Tz(:,:)     = T323(2,:,:)
      mx(:)       = mu23(Tx(:,:))
      my          = mu2(Ty(:))
      mz(:)       = mu32(Tz(:,:))
      kTy         = Cp * my * (-Ty(1) + Ty(2)) * dy(j) / Pr
      mty    = 0.5d0 * (mut(i,j,k) + mut(i,j+1,k))
      mtz(:) = (/0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i,j+1,k-1) + mut(i,j+1,k)), &
                0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i,j+1,k) + mut(i,j+1,k+1))/)
      mtx(:) = (/0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j+1,k) + mut(i,j+1,k)), &
                0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j+1,k) + mut(i+1,j+1,k))/)
      H(2:3) = (gamma * Q(5,i,j:j+1,k) / (Q(1,i,j:j+1,k) * gamma_1)) &
               + 0.5d0 * (Q(2,i,j:j+1,k)**2 + Q(3,i,j:j+1,k)**2 + Q(4,i,j:j+1,k)**2) + qc2(i,j:j+1,k)
      Hsgs   = -my * (-H(2) + H(3)) * dy(j) / Prt
      call calc_derivative_sgs_y2(nx, ny, nz, it, jt, kt, dx(i), dy(j), dz(k), my, mty, mx, mtx, mtz, mtz, u, v, w, tyy_tyysgs, tyx_tyxsgs, tyz_tyzsgs, vtyy, utyx, wtyz)
    endif
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx_tyxsgs
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy_tyysgs
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz_tyzsgs
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy + Hsgs)
  end subroutine calc_Fv_LES
  
  attributes(global) subroutine calc_Gv(nx, ny, nz, dx, dy, dz, Q, G)
    use calc_sutherland, only : mu6, mu2, mu32
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(inout), device :: G(5,nx-2,ny-2,nz-1)
    integer i, j, k, it, jt, kt
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz
    real(8), shared :: u(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), shared :: v(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), shared :: w(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    ! 4th-order 2
    real(8), device :: T6(6), mu(3)
    ! 2nd-order 2
    real(8) :: T332(3,3,2), Tx(3,2), Ty(3,2), Tz(2), mx(2), my(2), mz
    it = threadIdx%x + 1
    jt = threadIdx%y + 1
    kt = threadIdx%z 
    i = (blockIdx%x-1)*blockDim%x + it
    j = (blockIdx%y-1)*blockDim%y + jt
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      u(it-2:it+2,jt,kt-2:kt+3) = Q(2,i-2:i+2,j,k-2:k+3)
      v(it,jt-2:jt+2,kt-2:kt+3) = Q(3,i,j-2:j+2,k-2:k+3)
      w(it-2:it+2,jt-2:jt+2,kt-2:kt+3) = Q(4,i-2:i+2,j-2:j+2,k-2:k+3)
    endif
    call syncthreads()
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      T6(:) = Q(5,i,j,k-2:k+3) / (R * Q(1,i,j,k-2:k+3))
      mu(:) = mu6(T6(:))
      kTz   = heat_conduction6(mu(:), T6(:), dz(k))
      call calc_derivative_z4(nx, ny, nz, it, jt, kt, dx(i), dy(j), dz(k), mu, u, v, w, tzz, tzx, tzy, wtzz, utzx, vtzy)
    else
      T332(:,:,:) = Q(5,i-1:i+1,j-1:j+1,k:k+1) / (R * Q(1,i-1:i+1,j-1:j+1,k:k+1))
      Tx(:,:)     = T332(:,2,:)
      Ty(:,:)     = T332(2,:,:)
      Tz(:)       = T332(2,2,:)
      mx(:)       = mu32(Tx(:,:))
      my(:)       = mu32(Ty(:,:))
      mz          = mu2(Tz(:))
      kTz         = Cp * mz * (-Tz(1) + Tz(2)) * dz(k) / Pr
      call calc_derivative_z2(nx, ny, nz, it, jt, kt, dx(i), dy(j), dz(k), mz, mx, my, u, v, w, tzz, tzx, tzy, wtzz, utzx, vtzy)
    endif
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz)
  end subroutine calc_Gv

  attributes(global) subroutine calc_Gv_LES(nx, ny, nz, dx, dy, dz, Q, mut, qc2, G)
    use calc_sutherland, only : mu6, mu2, mu32
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(in), device    :: mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device :: G(5,nx-2,ny-2,nz-1)
    integer i, j, k, it, jt, kt
    real(8) :: tzx_tzxsgs, tzy_tzysgs, tzz_tzzsgs, utzx, vtzy, wtzz, kTz, mutz, H(4), Hsgs
    real(8), shared :: u(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), shared :: v(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), shared :: w(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    ! 4th-order 2
    real(8), device :: T6(6), mu(3), mt(3)
    ! 2nd-order 2
    real(8), device :: T332(3,3,2), Tx(3,2), Ty(3,2), Tz(2), mx(2), mtx(2), my(2), mty(2)
    real(8) mz, mtz
    it = threadIdx%x + 1
    jt = threadIdx%y + 1
    kt = threadIdx%z 
    i = (blockIdx%x-1)*blockDim%x + it
    j = (blockIdx%y-1)*blockDim%y + jt
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      u(it-2:it+2,jt,kt-2:kt+3) = Q(2,i-2:i+2,j,k-2:k+3)
      v(it,jt-2:jt+2,kt-2:kt+3) = Q(3,i,j-2:j+2,k-2:k+3)
      w(it-2:it+2,jt-2:jt+2,kt-2:kt+3) = Q(4,i-2:i+2,j-2:j+2,k-2:k+3)
    endif
    call syncthreads()
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      T6(:) = Q(5,i,j,k-2:k+3) / (R * Q(1,i,j,k-2:k+3))
      mu(:) = mu6(T6(:))
      kTz   = heat_conduction6(mu(:), T6(:), dz(k))
      mt(:) = 0.0625d0 * (-mut(i,j,k-1) + 9.d0 * (mut(i,j,k) + mut(i,j,k+1)) -mut(i,j,k+2))
      H(:)  = (gamma * Q(5,i,j,k-1:k+2) / (Q(1,i,j,k-1:k+2) * gamma_1)) &
              + 0.5d0 * (u(it,jt,kt-1:kt+2)**2 + v(it,jt,kt-1:kt+2)**2 + w(it,jt,kt-1:kt+2)**2) + qc2(i,j,k-1:k+2)
      Hsgs  = -mutz * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) / 3.d0) * dz(k) / Prt
      call calc_derivative_sgs_z4(nx, ny, nz, it, jt, kt, dx(i), dy(j), dz(k), mu, mt, u, v, w, tzz_tzzsgs, tzx_tzxsgs, tzy_tzysgs, wtzz, utzx, vtzy)
    else
      T332(:,:,:) = Q(5,i-1:i+1,j-1:j+1,k:k+1) / (R * Q(1,i-1:i+1,j-1:j+1,k:k+1))
      Tx(:,:)     = T332(:,2,:)
      Ty(:,:)     = T332(2,:,:)
      Tz(:)       = T332(2,2,:)
      mx(:)       = mu32(Tx(:,:))
      my(:)       = mu32(Ty(:,:))
      mz          = mu2(Tz(:))
      kTz         = Cp * mz * (-Tz(1) + Tz(2)) * dz(k) / Pr
      mtz    = 0.5d0 * (mut(i,j,k) + mut(i,j,k+1))
      mtx(:) = (/0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j,k+1) + mut(i,j,k+1)), &
                 0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j,k+1) + mut(i+1,j,k+1))/)
      mty(:) = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i,j-1,k+1) + mut(i,j,k+1)), &
                 0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i,j,k+1) + mut(i,j+1,k+1))/)
      H(2:3) = (gamma * Q(5,i,j,k:k+1) / (Q(1,i,j,k:k+1) * gamma_1)) &
               + 0.5d0 * (Q(2,i,j,k:k+1)**2 + Q(3,i,j,k:k+1)**2 + Q(4,i,j,k:k+1)**2) + qc2(i,j,k:k+1)
      Hsgs   = -mz * (-H(2) + H(3)) * dz(k) / Prt
      call calc_derivative_sgs_z2(nx, ny, nz, it, jt, kt, dx(i), dy(j), dz(k), mz, mtz, mx, mtx, my, mty, u, v, w, tzz_tzzsgs, tzx_tzxsgs, tzy_tzysgs, wtzz, utzx, vtzy)
    endif
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx_tzxsgs
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy_tzysgs
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz_tzzsgs
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz + Hsgs)
  end subroutine calc_Gv_LES
end module calc_visc

