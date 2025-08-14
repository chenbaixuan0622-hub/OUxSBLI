module calc_visc_visbal
  use mod_globals, only : R, Pr, Prt, gamma, threadsEv, threadsFv, threadsGv
  use mod_constant, only : Cp, gamma_1
  use calc_sutherland, only : mu2, mu23, mu32
  use calc_visc_base
  implicit none
contains
  attributes(device) subroutine calc_derivative_x2(nx, ny, nz, i, j, k, dx, dy, dz, Q, txx, txy, txz, utxx, vtxy, wtxz, kTx)
    integer, intent(in), value  :: nx, ny, nz, i, j, k
    real(8), intent(in), device :: dx(nx-1), dy(ny-1), dz(nz-1), Q(5,nx,ny,nz)
    real(8), intent(out)        :: txx, txy, txz, utxx, vtxy, wtxz, kTx
    real(8), dimension(2,3,3), device :: T233
    real(8), dimension(2,3), device   :: tmp2
    real(8), dimension(2), device     :: Tx, my, mz
    real(8) u, v, w, mx, mux, mvx, mwx, muy, mvy, muz, mwz
    ! dTdx & mu
    T233(:,:,:) = Q(5,i:i+1,j-1:j+1,k-1:k+1) / (R * Q(1,i:i+1,j-1:j+1,k-1:k+1))
    Tx(:)       = T233(:,2,2)
    mx          = mu2(Tx(:))
    kTx         = Cp * mx * (-Tx(1) + Tx(2)) * dx(i) / Pr
    tmp2(:,:)   = T233(:,:,2)
    my(:)       = mu23(tmp2(:,:))
    tmp2(:,:)   = T233(:,2,:)
    mz(:)       = mu23(tmp2(:,:))
    ! dudx & dudy 
    tmp2(:,:)   = Q(2,i:i+1,j-1:j+1,k)
    muy         = dy23(my(:), tmp2(:,:), dy(j))
    mux         = mx * (-tmp2(1,2) + tmp2(2,2)) * dx(i)
    u           = 0.5d0 * (tmp2(1,2) + tmp2(2,2))
    ! dudz
    tmp2(:,:)   = Q(2,i:i+1,j,k-1:k+1)
    muz         = dy23(mz(:), tmp2(:,:), dz(k))
    ! dvdx & dvdy
    tmp2(:,:)   = Q(3,i:i+1,j-1:j+1,k)
    mvy         = dy23(my(:), tmp2(:,:), dy(j))
    mvx         = mx * (-tmp2(1,2) + tmp2(2,2)) * dx(i)
    v           = 0.5d0 * (tmp2(1,2) + tmp2(2,2))
    ! dwdx & dwdz
    tmp2(:,:)   = Q(4,i:i+1,j,k-1:k+1)
    mwz         = dy23(mz(:), tmp2(:,:), dz(k))
    mwx         = mx * (-tmp2(1,2) + tmp2(2,2)) * dx(i)
    w           = 0.5d0 * (tmp2(1,2) + tmp2(2,2))
    txx         = 2.d0 * (2.d0 * mux - mvy - mwz) / 3.d0
    txy         = muy + mvx
    txz         = mwx + muz
    utxx        = u * txx
    vtxy        = v * txy
    wtxz        = w * txz
  end subroutine calc_derivative_x2
      
  attributes(device) subroutine calc_derivative_les_x2(nx, ny, nz, i, j, k, dx, dy, dz, Q, mut, qc2, txx, txy, txz, utxx, vtxy, wtxz, kTx, Hsgs)
    integer, intent(in), value  :: nx, ny, nz, i, j, k
    real(8), intent(in), device :: dx(nx-1), dy(ny-1), dz(nz-1), Q(5,nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(out)        :: txx, txy, txz, utxx, vtxy, wtxz, kTx, Hsgs
    real(8), dimension(2,3,3), device :: T233
    real(8), dimension(2,3), device   :: tmp2
    real(8), dimension(4), device     :: H
    real(8), dimension(2), device     :: Tx, my, mysgs, mz, mzsgs
    real(8) u, v, w, mx, mxsgs, mux, muxsgs, mvx, mvxsgs, mwx, mwxsgs, muy, muysgs, mvy, mvysgs, muz, muzsgs, mwz, mwzsgs
    ! dTdx & mu
    T233(:,:,:) = Q(5,i:i+1,j-1:j+1,k-1:k+1) / (R * Q(1,i:i+1,j-1:j+1,k-1:k+1))
    Tx(:)       = T233(:,2,2)
    mx          = mu2(Tx(:))
    kTx         = Cp * mx * (-Tx(1) + Tx(2)) * dx(i) / Pr
    tmp2(:,:)   = T233(:,:,2)
    my(:)       = mu23(tmp2(:,:))
    tmp2(:,:)   = T233(:,2,:)
    mz(:)       = mu23(tmp2(:,:))
    ! SGS
    mxsgs    = 0.5d0 * (mut(i,j,k) + mut(i+1,j,k))
    mysgs(:) = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i+1,j-1,k) + mut(i+1,j,k)), &
                 0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i+1,j,k) + mut(i+1,j+1,k))/)
    mzsgs(:) = (/0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i+1,j,k-1) + mut(i+1,j,k)), &
                 0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i+1,j,k) + mut(i+1,j,k+1))/)
    ! dudx & dudy 
    tmp2(:,:) = Q(2,i:i+1,j-1:j+1,k)
    muy       = dy23(my(:), tmp2(:,:), dy(j))
    muysgs    = dy23(mysgs(:), tmp2(:,:), dy(j))
    mux       = mx * (-tmp2(1,2) + tmp2(2,2)) * dx(i)
    muxsgs    = mxsgs * (-tmp2(1,2) + tmp2(2,2)) * dx(i)
    u         = 0.5d0 * (tmp2(1,2) + tmp2(2,2))
    ! dudz
    tmp2(:,:) = Q(2,i:i+1,j,k-1:k+1)
    muz       = dy23(mz(:), tmp2(:,:), dz(k))
    muzsgs    = dy23(mzsgs(:), tmp2(:,:), dz(k))
    ! dvdx & dvdy
    tmp2(:,:) = Q(3,i:i+1,j-1:j+1,k)
    mvy       = dy23(my(:), tmp2(:,:), dy(j))
    mvysgs    = dy23(mysgs(:), tmp2(:,:), dy(j))
    mvx       = mx * (-tmp2(1,2) + tmp2(2,2)) * dx(i)
    mvxsgs    = mxsgs * (-tmp2(1,2) + tmp2(2,2)) * dx(i)
    v         = 0.5d0 * (tmp2(1,2) + tmp2(2,2))
    ! dwdx & dwdz
    tmp2(:,:) = Q(4,i:i+1,j,k-1:k+1)
    mwz       = dy23(mz(:), tmp2(:,:), dz(k))
    mwx       = mx * (-tmp2(1,2) + tmp2(2,2)) * dx(i)
    mwxsgs    = mxsgs * (-tmp2(1,2) + tmp2(2,2)) * dx(i)
    mwzsgs    = dy23(mzsgs(:), tmp2(:,:), dz(k))
    w         = 0.5d0 * (tmp2(1,2) + tmp2(2,2))
    txx       = 2.d0 * (2.d0 * mux - mvy - mwz) / 3.d0
    txy       = muy + mvx
    txz       = mwx + muz
    utxx      = u * txx
    vtxy      = v * txy
    wtxz      = w 
    txx = txx + 2.d0 * (2.d0 * muxsgs - mvysgs - mwzsgs) / 3.d0
    txy = txy + muysgs + mvxsgs
    txz = txz + mwxsgs + muzsgs
    H(2:3) = (gamma * Q(5,i:i+1,j,k) / (Q(1,i:i+1,j,k) * gamma_1)) &
             + 0.5d0 * (Q(2,i:i+1,j,k)**2 + Q(3,i:i+1,j,k)**2 + Q(4,i:i+1,j,k)**2) + qc2(i:i+1,j,k)
    Hsgs   = -mx * (-H(2) + H(3)) * dx(i) / Prt
  end subroutine calc_derivative_les_x2
   
  attributes(device) subroutine calc_derivative_y2(nx, ny, nz, i, j, k, dx, dy, dz, Q, tyx, tyy, tyz, utyx, vtyy, wtyz, kTy)
    integer, intent(in), value  :: nx, ny, nz, i, j, k
    real(8), intent(in), device :: dx(nx-1), dy(ny-1), dz(nz-1), Q(5,nx,ny,nz)
    real(8), intent(out)        :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy
    real(8), dimension(3,2,3), device :: T323
    real(8), dimension(3,2), device   :: tmp2
    real(8), dimension(2,3), device   :: tmp3
    real(8), dimension(2), device     :: Ty, u2, v2, w2, mz, mx
    real(8) my, muy, mvy, mwy, mvz, mwz, mux, mvx
    ! dTdy
    T323(:,:,:) = Q(5,i-1:i+1,j:j+1,k-1:k+1) / (R * Q(1,i-1:i+1,j:j+1,k-1:k+1))
    tmp2(:,:)   = T323(:,:,2)
    mx(:)       = mu23(tmp2(:,:))
    Ty(:)       = T323(2,:,2)
    my          = mu2(Ty(:))
    kTy         = Cp * my * (-Ty(1) + Ty(2)) * dy(j) / Pr
    tmp3(:,:)   = T323(2,:,:)
    mz(:)       = mu32(tmp3(:,:))
    ! dudx & dudy
    tmp2(:,:)   = Q(2,i-1:i+1,j:j+1,k)
    mux         = dy32(mx(:), tmp2(:,:), dx(i))
    u2(:)       = tmp2(2,:)
    muy         = my * (-u2(1) + u2(2)) * dy(j)
    ! dvdx & dvdy
    tmp2(:,:)   = Q(3,i-1:i+1,j:j+1,k)
    mvx         = dy32(mx(:), tmp2(:,:), dx(i))
    v2(:)       = tmp2(2,:)
    mvy         = my * (-v2(1) + v2(2)) * dy(j)
    ! dvdz
    tmp3(:,:)   = Q(3,i,j:j+1,k-1:k+1)
    mvz         = dy23(mz(:), tmp3(:,:), dz(k))
    ! dwdz
    tmp3(:,:)   = Q(4,i,j:j+1,k-1:k+1)
    mwz         = dy23(mz(:), tmp3(:,:), dz(k))
    w2(:)       = tmp3(:,2)
    mwy         = my * (-w2(1) + w2(2)) * dy(j)
    tyx         = muy + mvx
    tyy         = 2.d0 * (2.d0 * mvy - mwz - mux) / 3.d0
    tyz         = mvz + mwy
    utyx        = 0.5d0 * (u2(1) + u2(2)) * tyx
    vtyy        = 0.5d0 * (v2(1) + v2(2)) * tyy
    wtyz        = 0.5d0 * (w2(1) + w2(2)) * tyz
  end subroutine calc_derivative_y2

  attributes(device) subroutine calc_derivative_les_y2(nx, ny, nz, i, j, k, dx, dy, dz, Q, mut, qc2, tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, Hsgs)
    integer, intent(in), value  :: nx, ny, nz, i, j, k
    real(8), intent(in), device :: dx(nx-1), dy(ny-1), dz(nz-1), Q(5,nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(out)        :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, Hsgs
    real(8), dimension(3,2,3), device :: T323
    real(8), dimension(3,2), device   :: tmp2
    real(8), dimension(2,3), device   :: tmp3
    real(8), dimension(4), device     :: H
    real(8), dimension(2), device     :: Ty, u2, v2, w2, mz, mzsgs, mx, mxsgs
    real(8) my, mysgs, muy, muysgs, mvy, mvysgs, mwy, mwysgs, mvz, mvzsgs, mwz, mwzsgs, mux, muxsgs, mvx, mvxsgs
    ! dTdy
    T323(:,:,:) = Q(5,i-1:i+1,j:j+1,k-1:k+1) / (R * Q(1,i-1:i+1,j:j+1,k-1:k+1))
    tmp2(:,:)   = T323(:,:,2)
    mx(:)       = mu23(tmp2(:,:))
    Ty(:)       = T323(2,:,2)
    my          = mu2(Ty(:))
    kTy         = Cp * my * (-Ty(1) + Ty(2)) * dy(j) / Pr
    tmp3(:,:)   = T323(2,:,:)
    mz(:)       = mu32(tmp3(:,:))
    ! SGS
    mysgs    = 0.5d0 * (mut(i,j,k) + mut(i,j+1,k))
    mzsgs(:) = (/0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i,j+1,k-1) + mut(i,j+1,k)), &
                 0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i,j+1,k) + mut(i,j+1,k+1))/)
    mxsgs(:) = (/0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j+1,k) + mut(i,j+1,k)), &
                 0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j+1,k) + mut(i+1,j+1,k))/)
    ! dudx & dudy
    tmp2(:,:) = Q(2,i-1:i+1,j:j+1,k)
    mux       = dy32(mx(:), tmp2(:,:), dx(i))
    muxsgs    = dy32(mxsgs(:), tmp2(:,:), dx(i))
    u2(:)     = tmp2(2,:)
    muy       = my * (-u2(1) + u2(2)) * dy(j)
    muysgs    = my * (-u2(1) + u2(2)) * dy(j)
    ! dvdx & dvdy
    tmp2(:,:) = Q(3,i-1:i+1,j:j+1,k)
    mvx       = dy32(mx(:), tmp2(:,:), dx(i))
    mvxsgs    = dy32(mxsgs(:), tmp2(:,:), dx(i))
    v2(:)     = tmp2(2,:)
    mvy       = my * (-v2(1) + v2(2)) * dy(j)
    mvysgs    = my * (-v2(1) + v2(2)) * dy(j)
    ! dvdz
    tmp3(:,:) = Q(3,i,j:j+1,k-1:k+1)
    mvz       = dy23(mz(:), tmp3(:,:), dz(k))
    mvzsgs    = dy23(mzsgs(:), tmp3(:,:), dz(k))
    ! dwdz
    tmp3(:,:) = Q(4,i,j:j+1,k-1:k+1)
    mwz       = dy23(mz(:), tmp3(:,:), dz(k))
    mwzsgs    = dy23(mzsgs(:), tmp3(:,:), dz(k))
    w2(:)     = tmp3(:,2)
    mwy       = my * (-w2(1) + w2(2)) * dy(j)
    mwysgs    = my * (-w2(1) + w2(2)) * dy(j)
    tyx       = muy + mvx
    tyy       = 2.d0 * (2.d0 * mvy - mwz - mux) / 3.d0
    tyz       = mvz + mwy
    utyx      = 0.5d0 * (u2(1) + u2(2)) * tyx
    vtyy      = 0.5d0 * (v2(1) + v2(2)) * tyy
    wtyz      = 0.5d0 * (w2(1) + w2(2)) * tyz
    tyx = tyx + muy + mvx
    tyy = tyy + 2.d0 * (2.d0 * mvy - mwz - mux) / 3.d0
    tyz = tyz + mvz + mwy
    H(2:3) = (gamma * Q(5,i,j:j+1,k) / (Q(1,i,j:j+1,k) * gamma_1)) &
             + 0.5d0 * (Q(2,i,j:j+1,k)**2 + Q(3,i,j:j+1,k)**2 + Q(4,i,j:j+1,k)**2) + qc2(i,j:j+1,k)
    Hsgs   = -my * (-H(2) + H(3)) * dy(j) / Prt
  end subroutine calc_derivative_les_y2

  attributes(device) subroutine calc_derivative_z2(nx, ny, nz, i, j, k, dx, dy, dz, Q, tzx, tzy, tzz, utzx, vtzy, wtzz, kTz)
    integer, intent(in), value  :: nx, ny, nz, i, j, k
    real(8), intent(in), device :: dx(nx-1), dy(ny-1), dz(nz-1), Q(5,nx,ny,nz)
    real(8), intent(out)        :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz
    real(8), dimension(3,3,2), device :: T332
    real(8), dimension(3,2), device   :: tmp2
    real(8), dimension(2), device     :: Tz, u2, v2, w2, mx, my
    real(8) mz, muz, mvz, mwz, mwx, mux, mvy, mwy
    ! dTdz & mu
    T332(:,:,:) = Q(5,i-1:i+1,j-1:j+1,k:k+1) / (R * Q(1,i-1:i+1,j-1:j+1,k:k+1))
    tmp2(:,:)   = T332(:,2,:)
    mx(:)       = mu32(tmp2(:,:))
    tmp2(:,:)   = T332(2,:,:)
    my(:)       = mu32(tmp2(:,:))
    Tz(:)       = T332(2,2,:)
    mz          = mu2(Tz(:))
    !dudx & dudz
    tmp2(:,:)   = Q(2,i-1:i+1,j,k:k+1)
    mux         = dy32(mx(:), tmp2(:,:), dx(i))
    u2(:)       = tmp2(2,:)
    muz         = mz * (-u2(1) + u2(2)) * dz(k)
    ! dvdy & dvdz
    tmp2(:,:)   = Q(3,i,j-1:j+1,k:k+1)
    mvy         = dy32(my(:), tmp2(:,:), dy(j))
    v2(:)       = tmp2(2,:)
    mvz         = mz * (-v2(1) + v2(2)) * dz(k)
    ! dwdx
    tmp2(:,:)   = Q(4,i-1:i+1,j,k:k+1)
    mwx         = dy32(mx(:), tmp2(:,:), dx(i))
    ! dwdz & dwdy
    tmp2(:,:)   = Q(4,i,j-1:j+1,k:k+1)
    w2(:)       = tmp2(2,:)
    mwz         = mz * (-w2(1) + w2(2)) * dz(k)
    mwy         = dy32(my(:), tmp2(:,:), dy(j))
    tzx         = mwx + muz
    tzy         = mvz + mwy
    tzz         = 2.d0 * (2.d0 * mwz - mux - mvy) / 3.d0
    utzx        = 0.5d0 * (u2(1) + u2(2)) * tzx
    vtzy        = 0.5d0 * (v2(1) + v2(2)) * tzy
    wtzz        = 0.5d0 * (w2(1) + w2(2)) * tzz
    kTz         = Cp * mz * (-Tz(1) + Tz(2)) * dz(k) / Pr
  end subroutine calc_derivative_z2

  attributes(device) subroutine calc_derivative_les_z2(nx, ny, nz, i, j, k, dx, dy, dz, Q, mut, qc2, tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, Hsgs)
    integer, intent(in), value  :: nx, ny, nz, i, j, k
    real(8), intent(in), device :: dx(nx-1), dy(ny-1), dz(nz-1), Q(5,nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(out)        :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, Hsgs
    real(8), dimension(3,3,2), device :: T332
    real(8), dimension(3,2), device   :: tmp2
    real(8), dimension(4), device     :: H
    real(8), dimension(2), device     :: Tz, u2, v2, w2, mx, mxsgs, my, mysgs
    real(8) mz, mzsgs, muz, muzsgs, mvz, mvzsgs, mwz, mwzsgs, mwx, mwxsgs, mux, muxsgs, mvy, mvysgs, mwy, mwysgs
    T332(:,:,:) = Q(5,i-1:i+1,j-1:j+1,k:k+1) / (R * Q(1,i-1:i+1,j-1:j+1,k:k+1))
    tmp2(:,:)   = T332(:,2,:)
    mx(:)       = mu32(tmp2(:,:))
    tmp2(:,:)   = T332(2,:,:)
    my(:)       = mu32(tmp2(:,:))
    Tz(:)       = T332(2,2,:)
    mz          = mu2(Tz(:))
    ! SGS
    mzsgs    = 0.5d0 * (mut(i,j,k) + mut(i,j,k+1))
    mxsgs(:) = (/0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j,k+1) + mut(i,j,k+1)), &
                 0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j,k+1) + mut(i+1,j,k+1))/)
    mysgs(:) = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i,j-1,k+1) + mut(i,j,k+1)), &
                 0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i,j,k+1) + mut(i,j+1,k+1))/)
    ! dudx & dudz
    tmp2(:,:) = Q(2,i-1:i+1,j,k:k+1)
    mux       = dy32(mx(:), tmp2(:,:), dx(i))
    muxsgs    = dy32(mxsgs(:), tmp2, dx(i))
    u2(:)     = tmp2(2,:)
    muz       = mz * (-u2(1) + u2(2)) * dz(k)
    muzsgs    = mzsgs * (-u2(1) + u2(2)) * dz(k)
    ! dvdy & dvdz
    tmp2(:,:) = Q(3,i,j-1:j+1,k:k+1)
    mvy       = dy32(my(:), tmp2(:,:), dy(j))
    mvysgs    = dy32(mysgs(:), tmp2(:,:), dy(j))
    v2(:)     = tmp2(2,:)
    mvz       = mz * (-v2(1) + v2(2)) * dz(k)
    mvzsgs    = mzsgs * (-v2(1) + v2(2)) * dz(k)
    ! dwdx
    tmp2(:,:) = Q(4,i-1:i+1,j,k:k+1)
    mwx       = dy32(mx(:), tmp2(:,:), dx(i))
    mwxsgs    = dy32(mxsgs(:), tmp2(:,:), dx(i))
    ! dwdz & dwdy
    tmp2(:,:) = Q(4,i,j-1:j+1,k:k+1)
    w2(:)     = tmp2(2,:)
    mwz       = mz * (-w2(1) + w2(2)) * dz(k)
    mwzsgs    = mzsgs * (-w2(1) + w2(2)) * dz(k)
    mwy       = dy32(my(:), tmp2(:,:), dy(j))
    mwysgs    = dy32(mysgs(:), tmp2(:,:), dy(j))
    tzx       = mwx + muz
    tzy       = mvz + mwy
    tzz       = 2.d0 * (2.d0 * mwz - mux - mvy) / 3.d0
    utzx      = 0.5d0 * (u2(1) + u2(2)) * tzx
    vtzy      = 0.5d0 * (v2(1) + v2(2)) * tzy
    wtzz      = 0.5d0 * (w2(1) + w2(2)) * tzz
    kTz       = Cp * mz * (-Tz(1) + Tz(2)) * dz(k) / Pr
    tzx = tzx + mwxsgs + muzsgs
    tzy = tzy + mvzsgs + mwysgs
    tzz = tzz + 2.d0 * (2.d0 * mwzsgs - muxsgs - mvysgs) / 3.d0
    H(2:3) = (gamma * Q(5,i,j,k:k+1) / (Q(1,i,j,k:k+1) * gamma_1)) &
             + 0.5d0 * (Q(2,i,j,k:k+1)**2 + Q(3,i,j,k:k+1)**2 + Q(4,i,j,k:k+1)**2) + qc2(i,j,k:k+1)
    Hsgs   = -mz * (-H(2) + H(3)) * dz(k) / Prt
  end subroutine calc_derivative_les_z2
end module calc_visc_visbal

