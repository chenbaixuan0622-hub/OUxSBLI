module calc_visc
  use mod_globals, only : id_visc, gamma, R, Pr, Prt
  use mod_constant, only : Cp, gamma_1
  use calc_visc_common
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
    integer i, j, k
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc ==2 .and. 3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8), dimension(6), device :: u6, v6, w6
        real(8), dimension(3), device :: ux3, vx3, wx3, uy3, vy3, uz3, wz3, mu
        block
          real(8), device :: u651(6,5)
          u651(:,:) = Q(2,i-2:i+3,j-2:j+2,k)
          uy3(:)    = dy65(u651(:,:), dy(j))
          u6(:)     = u651(:,3) 
          ux3(:)    = dx6(u6(:), dx(i))
        end block
        block
          real(8), device :: u615(6,5)
          u615(:,:) = Q(2,i-2:i+3,j,k-2:k+2)
          uz3(:)    = dy65(u615(:,:), dz(k))
        end block
        block
          real(8), device :: v651(6,5)
          v651(:,:) = Q(3,i-2:i+3,j-2:j+2,k)
          vy3(:)    = dy65(v651(:,:), dy(j))
          v6(:)     = v651(:,3) 
          vx3(:)    = dx6(v6(:), dx(i))
        end block
        block
          real(8), device :: w615(6,5)
          w615(:,:) = Q(4,i-2:i+3,j,k-2:k+2)
          wz3(:)    = dy65(w615(:,:), dz(k))
          w6(:)     = w615(:,3) 
          wx3(:)    = dx6(w6(:), dx(i))
        end block
        block
          real(8), device :: T6(6)
          T6(:) = Q(5,i-2:i+3,j,k) / (R * Q(1,i-2:i+3,j,k))
          mu(:) = mu6(T6(:))
          kTx   = heat_conduction6(mu(:), T6(:), dx(i))
        end block
        call tauxx4(mu(:), ux3(:), vy3(:), wz3(:), u6(:), txx, utxx)
        call tauxy4(mu(:), uy3(:), vx3(:), v6(:), txy, vtxy)
        call tauxy4(mu(:), wx3(:), uz3(:), w6(:), txz, wtxz)
      end block
    else
      block
        real(8), dimension(2), device :: u2, v2, w2, my, mz
        real(8) mx, mux, mvx, mwx, muy, mvy, muz, mwz
        block
          real(8), device :: T233(2,3,3)
          T233(:,:,:) = Q(5,i:i+1,j-1:j+1,k-1:k+1) / (R * Q(1,i:i+1,j-1:j+1,k-1:k+1))
          block
            real(8), device :: Tx(2)
            Tx(:) = T233(:,2,2)
            mx    = mu2(Tx(:))
            kTx   = Cp * mx * (-Tx(1) + Tx(2)) * dx(i) / Pr
          end block
          block
            real(8), device :: Ty(2,3)
            Ty(:,:) = T233(:,:,2)
            my(:)   = mu23(Ty(:,:))
          end block
          block
            real(8), device :: Tz(2,3)
            Tz(:,:) = T233(:,2,:)
            mz(:)   = mu23(Tz(:,:))
          end block
        end block
        block
          real(8), device :: u231(2,3)
          u231(:,:) = Q(2,i:i+1,j-1:j+1,k)
          muy       = dy23(my(:), u231(:,:), dy(j))
          u2(:)     = u231(:,2)
          mux       = mx * (-u2(1) + u2(2)) * dx(i)
        end block
        block
          real(8), device :: u213(2,3)
          u213(:,:) = Q(2,i:i+1,j,k-1:k+1)
          muz       = dy23(mz(:), u213(:,:), dz(k))
        end block
        block
          real(8), device :: v231(2,3)
          v231(:,:) = Q(3,i:i+1,j-1:j+1,k)
          mvy       = dy23(my(:), v231(:,:), dy(j))
          v2(:)     = v231(:,2)
          mvx       = mx * (-v2(1) + v2(2)) * dx(i)
        end block
        block
          real(8), device :: w213(2,3)
          w213(:,:) = Q(4,i:i+1,j,k-1:k+1)
          mwz       = dy23(mz(:), w213(:,:), dz(k))
          w2(:)     = w213(:,2)
          mwx       = mx * (-w2(1) + w2(2)) * dx(i)
        end block
        txx  = 2.d0 * (2.d0 * mux - mvy - mwz) / 3.d0
        txy  = muy + mvx
        txz  = mwx + muz
        utxx = 0.5d0 * (u2(1) + u2(2)) * txx
        vtxy = 0.5d0 * (v2(1) + v2(2)) * txy
        wtxz = 0.5d0 * (w2(1) + w2(2)) * txz
      end block
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
    integer i, j, k
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx, Hsgs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc ==2 .and. 3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8), dimension(6), device :: u6, v6, w6
        real(8), dimension(3), device :: ux3, vx3, wx3, uy3, vy3, uz3, wz3, mu
        block
          real(8), device :: u651(6,5)
          u651(:,:) = Q(2,i-2:i+3,j-2:j+2,k)
          uy3(:)    = dy65(u651(:,:), dy(j))
          u6(:)     = u651(:,3) 
          ux3(:)    = dx6(u6(:), dx(i))
        end block
        block
          real(8), device :: u615(6,5)
          u615(:,:) = Q(2,i-2:i+3,j,k-2:k+2)
          uz3(:)    = dy65(u615(:,:), dz(k))
        end block
        block
          real(8), device :: v651(6,5)
          v651(:,:) = Q(3,i-2:i+3,j-2:j+2,k)
          vy3(:)    = dy65(v651(:,:), dy(j))
          v6(:)     = v651(:,3) 
          vx3(:)    = dx6(v6(:), dx(i))
        end block
        block
          real(8), device :: w615(6,5)
          w615(:,:) = Q(4,i-2:i+3,j,k-2:k+2)
          wz3(:)    = dy65(w615(:,:), dz(k))
          w6(:)     = w615(:,3) 
          wx3(:)    = dx6(w6(:), dx(i))
        end block
        block
          real(8), device :: T6(6)
          T6(:) = Q(5,i-2:i+3,j,k) / (R * Q(1,i-2:i+3,j,k))
          mu(:) = mu6(T6(:))
          kTx   = heat_conduction6(mu(:), T6(:), dx(i))
        end block
        call tauxx4(mu(:), ux3(:), vy3(:), wz3(:), u6(:), txx, utxx)
        call tauxy4(mu(:), uy3(:), vx3(:), v6(:), txy, vtxy)
        call tauxy4(mu(:), wx3(:), uz3(:), w6(:), txz, wtxz)
        block
          real(8) :: mutx, H(4)
          mutx = 0.0625d0 * (-mut(i-1,j,k) + 9.d0 * (mut(i,j,k) + mut(i+1,j,k)) -mut(i+2,j,k))
          txx  = txx + 2.d0 * mutx * (2.d0 * ux3(2) - vy3(2) - wz3(2)) / 3.d0
          txy  = txy + mutx * (uy3(2) + vx3(2))
          txz  = txz + mutx * (wx3(2) + uz3(2))
          H(:) = (gamma * Q(5,i-1:i+2,j,k) / (Q(1,i-1:i+2,j,k) * gamma_1)) &
                 + 0.5d0 * (Q(2,i-1:i+2,j,k)**2 + Q(3,i-1:i+2,j,k)**2 + Q(4,i-1:i+2,j,k)**2) + qc2(i-1:i+2,j,k)
          Hsgs = -mutx * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) / 3.d0) * dx(i) / Prt
        end block
      end block
    else
      block
        real(8), dimension(2), device :: u2, v2, w2, my, mysgs, mz, mzsgs
        real(8) mx, mxsgs, mux, muxsgs, mvx, mvxsgs, mwx, mwxsgs, muy, muysgs, mvy, mvysgs, muz, muzsgs, mwz, mwzsgs
        ! SGS
        mx    = 0.5d0 * (mut(i,j,k) + mut(i+1,j,k))
        my(:) = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i+1,j-1,k) + mut(i+1,j,k)), &
                  0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i+1,j,k) + mut(i+1,j+1,k))/)
        mz(:) = (/0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i+1,j,k-1) + mut(i+1,j,k)), &
                  0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i+1,j,k) + mut(i+1,j,k+1))/)
        block
          real(8), device :: T233(2,3,3)
          T233(:,:,:) = Q(5,i:i+1,j-1:j+1,k-1:k+1) / (R * Q(1,i:i+1,j-1:j+1,k-1:k+1))
          block
            real(8), device :: Tx(2)
            Tx(:) = T233(:,2,2)
            mx    = mu2(Tx(:))
            kTx   = Cp * mx * (-Tx(1) + Tx(2)) * dx(i) / Pr
          end block
          block
            real(8), device :: Ty(2,3)
            Ty(:,:) = T233(:,:,2)
            my(:)   = mu23(Ty(:,:))
          end block
          block
            real(8), device :: Tz(2,3)
            Tz(:,:) = T233(:,2,:)
            mz(:)   = mu23(Tz(:,:))
          end block
        end block
        block
          real(8), device :: u231(2,3)
          u231(:,:) = Q(2,i:i+1,j-1:j+1,k)
          muy       = dy23(my(:), u231(:,:), dy(j))
          muysgs    = dy23(mysgs(:), u231(:,:), dy(j))
          u2(:)     = u231(:,2)
          mux       = mx * (-u2(1) + u2(2)) * dx(i)
          muxsgs    = mxsgs * (-u2(1) + u2(2)) * dx(i)
        end block
        block
          real(8), device :: u213(2,3)
          u213(:,:) = Q(2,i:i+1,j,k-1:k+1)
          muz       = dy23(mz(:), u213(:,:), dz(k))
          muzsgs    = dy23(mzsgs(:), u213(:,:), dz(k))
        end block
        block
          real(8), device :: v231(2,3)
          v231(:,:) = Q(3,i:i+1,j-1:j+1,k)
          mvy       = dy23(my(:), v231(:,:), dy(j))
          mvysgs    = dy23(mysgs(:), v231(:,:), dy(j))
          v2(:)     = v231(:,2)
          mvx       = mx * (-v2(1) + v2(2)) * dx(i)
          mvxsgs    = mxsgs * (-v2(1) + v2(2)) * dx(i)
        end block
        block
          real(8), device :: w213(2,3)
          w213(:,:) = Q(4,i:i+1,j,k-1:k+1)
          mwz       = dy23(mz(:), w213(:,:), dz(k))
          mwzsgs    = dy23(mzsgs(:), w213(:,:), dz(k))
          w2(:)     = w213(:,2)
          mwx       = mx * (-w2(1) + w2(2)) * dx(i)
          mwxsgs    = mxsgs * (-w2(1) + w2(2)) * dx(i)
        end block
        txx  = 2.d0 * (2.d0 * mux - mvy - mwz) / 3.d0
        txy  = muy + mvx
        txz  = mwx + muz
        utxx = 0.5d0 * (u2(1) + u2(2)) * txx
        vtxy = 0.5d0 * (v2(1) + v2(2)) * txy
        txx = txx + 2.d0 * (2.d0 * muxsgs - mvysgs - mwzsgs) / 3.d0
        txy = txy + muysgs + mvxsgs
        txz = txz + mwxsgs + muzsgs
        block
          real(8), device :: H(2)
          H(:) = (gamma * Q(5,i:i+1,j,k) / (Q(1,i:i+1,j,k) * gamma_1)) &
                 + 0.5d0 * (Q(2,i:i+1,j,k)**2 + Q(3,i:i+1,j,k)**2 + Q(4,i:i+1,j,k)**2) + qc2(i:i+1,j,k)
          Hsgs = -mx * (-H(1) + H(2)) * dx(i) / Prt
        end block
      end block
    endif
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
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
    integer i, j, k
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8), dimension(6), device :: u6, v6, w6
        real(8), dimension(3), device :: uy3, vy3, wy3, vz3, wz3, ux3, vx3, mu
        block
          real(8), device :: u561(5,6)
          u561(:,:) = Q(2,i-2:i+2,j-2:j+3,k)
          ux3(:)    = dy56(u561(:,:), dx(i))
          u6(:)     = u561(3,:)
          uy3(:)    = dx6(u6(:), dy(j))
        end block
        block
          real(8), device :: v561(5,6)
          v561(:,:) = Q(3,i-2:i+2,j-2:j+3,k)
          vx3(:)    = dy56(v561(:,:), dx(i))
          v6(:)     = v561(3,:)
          vy3(:)    = dx6(v6(:), dy(j))
        end block
        block
          real(8), device :: v165(6,5)
          v165(:,:) = Q(3,i,j-2:j+3,k-2:k+2)
          vz3(:)    = dy65(v165(:,:), dz(k))
        end block
        block
          real(8), device :: w165(6,5)
          w165(:,:) = Q(4,i,j-2:j+3,k-2:k+2)
          w6(:)     = w165(:,3)
          wz3(:)    = dy65(w165(:,:), dz(k))
          wy3(:)    = dx6(w6(:), dy(j))
        end block
        block
          real(8), device :: T6(6)
          T6(:) = Q(5,i,j-2:j+3,k) / (R * Q(1,i,j-2:j+3,k))
          mu(:) = mu6(T6(:))
          kTy   = heat_conduction6(mu(:), T6(:), dy(j))
        end block
        call tauxy4(mu(:), uy3(:), vx3(:), u6(:), tyx, utyx)
        call tauxx4(mu(:), vy3(:), wz3(:), ux3(:), v6(:), tyy, vtyy)
        call tauxy4(mu(:), vz3(:), wy3(:), w6(:), tyz, wtyz)
      end block
    else
      block
        real(8), dimension(2), device :: u2, v2, w2, mz, mx
        real(8) my, muy, mvy, mwy, mvz, mwz, mux, mvx
        block
          real(8), device :: T323(3,2,3)
          T323(:,:,:) = Q(5,i-1:i+1,j:j+1,k-1:k+1) / (R * Q(1,i-1:i+1,j:j+1,k-1:k+1))
          block
            real(8), device :: Tx(3,2)
            Tx(:,:) = T323(:,:,2)
            mx(:)   = mu23(Tx(:,:))
          end block
          block
            real(8), device :: Ty(2)
            Ty(:) = T323(2,:,2)
            my    = mu2(Ty(:))
            kTy   = Cp * my * (-Ty(1) + Ty(2)) * dy(j) / Pr
          end block
          block
            real(8), device :: Tz(2,3)
            Tz(:,:) = T323(2,:,:)
            mz(:)   = mu32(Tz(:,:))
          end block
        end block
        block
          real(8), device :: u321(3,2)
          u321(:,:) = Q(2,i-1:i+1,j:j+1,k)
          mux       = dy32(mx(:), u321(:,:), dx(i))
          u2(:)     = u321(2,:)
          muy       = my * (-u2(1) + u2(2)) * dy(j)
        end block
        block
          real(8), device :: v321(3,2)
          v321(:,:) = Q(3,i-1:i+1,j:j+1,k)
          mvx       = dy32(mx(:), v321(:,:), dx(i))
          v2(:)     = v321(2,:)
          mvy       = my * (-v2(1) + v2(2)) * dy(j)
        end block
        block
          real(8), device :: v123(2,3)
          v123(:,:) = Q(3,i,j:j+1,k-1:k+1)
          mvz       = dy23(mz(:), v123(:,:), dz(k))
        end block
        block
          real(8), device :: w123(2,3)
          w123(:,:) = Q(4,i,j:j+1,k-1:k+1)
          mwz       = dy23(mz(:), w123(:,:), dz(k))
          w2(:)     = w123(:,2)
          mwy       = my * (-w2(1) + w2(2)) * dy(j)
        end block
        tyx  = muy + mvx
        tyy  = 2.d0 * (2.d0 * mvy - mwz - mux) / 3.d0
        tyz  = mvz + mwy
        utyx = 0.5d0 * (u2(1) + u2(2)) * tyx
        vtyy = 0.5d0 * (v2(1) + v2(2)) * tyy
        wtyz = 0.5d0 * (w2(1) + w2(2)) * tyz
      end block
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
    integer i, j, k
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, Hsgs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8), dimension(6), device   :: u6, v6, w6
        real(8), dimension(3), device   :: uy3, vy3, wy3, vz3, wz3, ux3, vx3, mu
        block
          real(8), device :: u561(5,6)
          u561(:,:) = Q(2,i-2:i+2,j-2:j+3,k)
          ux3(:)    = dy56(u561(:,:), dx(i))
          u6(:)     = u561(3,:)
          uy3(:)    = dx6(u6(:), dy(j))
        end block
        block
          real(8), device :: v561(5,6)
          v561(:,:) = Q(3,i-2:i+2,j-2:j+3,k)
          vx3(:)    = dy56(v561(:,:), dx(i))
          v6(:)     = v561(3,:)
          vy3(:)    = dx6(v6(:), dy(j))
        end block
        block
          real(8), device :: v165(6,5)
          v165(:,:) = Q(3,i,j-2:j+3,k-2:k+2)
          vz3(:)    = dy65(v165(:,:), dz(k))
        end block
        block
          real(8), device :: w165(6,5)
          w165(:,:) = Q(4,i,j-2:j+3,k-2:k+2)
          w6(:)     = w165(:,3)
          wz3(:)    = dy65(w165(:,:), dz(k))
          wy3(:)    = dx6(w6(:), dy(j))
        end block
        block
          real(8), device :: T6(6)
          T6(:) = Q(5,i,j-2:j+3,k) / (R * Q(1,i,j-2:j+3,k))
          mu(:) = mu6(T6(:))
          kTy   = heat_conduction6(mu(:), T6(:), dy(j))
        end block
        call tauxy4(mu(:), uy3(:), vx3(:), u6(:), tyx, utyx)
        call tauxx4(mu(:), vy3(:), wz3(:), ux3(:), v6(:), tyy, vtyy)
        call tauxy4(mu(:), vz3(:), wy3(:), w6(:), tyz, wtyz)
        block
          real(8), device :: H(4)
          real(8) muty
          muty   = 0.0625d0 * (-mut(i,j-1,k) + 9.d0 * (mut(i,j,k) + mut(i,j+1,k)) -mut(i,j+2,k))
          tyx = tyx + muty * (uy3(2) + vx3(2))
          tyy = tyy + 2.d0 * muty * (2.d0 * vy3(2) - ux3(2) - wz3(2)) / 3.d0
          tyz = tyz + muty * (vz3(2) + wy3(2))
          H(:)   = (gamma * Q(5,i,j-1:j+2,k) / (Q(1,i,j-1:j+2,k) * gamma_1)) &
                   + 0.5d0 * (Q(2,i,j-1:j+2,k)**2 + Q(3,i,j-1:j+2,k)**2 + Q(4,i,j-1:j+2,k)**2) + qc2(i,j-1:j+2,k)
          Hsgs   = -muty * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) / 3.d0) * dy(j) / Prt
        end block
      end block
    else
      block
        real(8), dimension(2), device :: u2, v2, w2, mz, mzsgs, mx, mxsgs
        real(8) my, mysgs, muy, muysgs, mvy, mvysgs, mwy, mwysgs, mvz, mvzsgs, mwz, mwzsgs, mux, muxsgs, mvx, mvxsgs
        ! SGS
        my     = 0.5d0 * (mut(i,j,k) + mut(i,j+1,k))
        mz(:)  = (/0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i,j+1,k-1) + mut(i,j+1,k)), &
                   0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i,j+1,k) + mut(i,j+1,k+1))/)
        mx(:)  = (/0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j+1,k) + mut(i,j+1,k)), &
                   0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j+1,k) + mut(i+1,j+1,k))/)
        block
          real(8), device :: T323(3,2,3)
          T323(:,:,:) = Q(5,i-1:i+1,j:j+1,k-1:k+1) / (R * Q(1,i-1:i+1,j:j+1,k-1:k+1))
          block
            real(8), device :: Tx(3,2)
            Tx(:,:) = T323(:,:,2)
            mx(:)   = mu23(Tx(:,:))
          end block
          block
            real(8), device :: Ty(2)
            Ty(:) = T323(2,:,2)
            my    = mu2(Ty(:))
            kTy   = Cp * my * (-Ty(1) + Ty(2)) * dy(j) / Pr
          end block
          block
            real(8), device :: Tz(2,3)
            Tz(:,:) = T323(2,:,:)
            mz(:)   = mu32(Tz(:,:))
          end block
        end block
        block
          real(8), device :: u321(3,2)
          u321(:,:) = Q(2,i-1:i+1,j:j+1,k)
          mux       = dy32(mx(:), u321(:,:), dx(i))
          muxsgs    = dy32(mxsgs(:), u321(:,:), dx(i))
          u2(:)     = u321(2,:)
          muy       = my * (-u2(1) + u2(2)) * dy(j)
          muysgs    = mysgs * (-u2(1) + u2(2)) * dy(j)
        end block
        block
          real(8), device :: v321(3,2)
          v321(:,:) = Q(3,i-1:i+1,j:j+1,k)
          mvx       = dy32(mx(:), v321(:,:), dx(i))
          mvxsgs    = dy32(mxsgs(:), v321(:,:), dx(i))
          v2(:)     = v321(2,:)
          mvy       = my * (-v2(1) + v2(2)) * dy(j)
          mvysgs    = mysgs * (-v2(1) + v2(2)) * dy(j)
        end block
        block
          real(8), device :: v123(2,3)
          v123(:,:) = Q(3,i,j:j+1,k-1:k+1)
          mvz       = dy23(mz(:), v123(:,:), dz(k))
          mvzsgs    = dy23(mzsgs(:), v123(:,:), dz(k))
        end block
        block
          real(8), device :: w123(2,3)
          w123(:,:) = Q(4,i,j:j+1,k-1:k+1)
          mwz       = dy23(mz(:), w123(:,:), dz(k))
          mwzsgs    = dy23(mzsgs(:), w123(:,:), dz(k))
          w2(:)     = w123(:,2)
          mwy       = my * (-w2(1) + w2(2)) * dy(j)
          mwysgs    = mysgs * (-w2(1) + w2(2)) * dy(j)
        end block
          tyx  = muy + mvx
          tyy  = 2.d0 * (2.d0 * mvy - mwz - mux) / 3.d0
          tyz  = mvz + mwy
          utyx = 0.5d0 * (u2(1) + u2(2)) * tyx
          vtyy = 0.5d0 * (v2(1) + v2(2)) * tyy
          wtyz = 0.5d0 * (w2(1) + w2(2)) * tyz
          tyx  = tyx + muysgs + mvxsgs
          tyy  = tyy + 2.d0 * (2.d0 * mvysgs - mwzsgs - muxsgs) / 3.d0
          tyz  = tyz + mvzsgs + mwysgs
        block
          real(8), device :: H(2)
          H(:) = (gamma * Q(5,i,j:j+1,k) / (Q(1,i,j:j+1,k) * gamma_1)) &
                 + 0.5d0 * (u2(:)**2 + v2(:)**2 + w2(:)**2) + qc2(i,j:j+1,k)
          Hsgs = -my * (-H(1) + H(2)) * dy(j) / Prt
        end block
      end block
    endif
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
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
    integer i, j, k
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      block
        real(8), dimension(6), device :: u6, v6, w6
        real(8), dimension(3), device :: uz3, vz3, wz3, wx3, ux3, vy3, wy3, mu
        block
          real(8), device :: u516(5,6)
          u516(:,:) = Q(2,i-2:i+2,j,k-2:k+3)
          ux3(:)    = dy56(u516(:,:), dx(i))
          u6(:)     = u516(3,:)
          uz3(:)    = dx6(u6(:), dz(k))
        end block
        block
          real(8), device :: v156(5,6)
          v156(:,:) = Q(3,i,j-2:j+2,k-2:k+3)
          vy3(:)    = dy56(v156(:,:), dy(j))
          v6(:)     = v156(3,:)
          vz3(:)    = dx6(v6(:), dz(k))
        end block
        block
          real(8), device :: w516(5,6)
          w516(:,:) = Q(4,i-2:i+2,j,k-2:k+3)
          wx3(:)    = dy56(w516(:,:), dx(i))
        end block
        block
          real(8), device :: w156(5,6)
          w156(:,:) = Q(4,i,j-2:j+2,k-2:k+3)
          wy3(:)    = dy56(w156(:,:), dy(j))
          w6(:)     = w156(3,:)
          wz3(:)    = dx6(w6(:), dz(k))
        end block
        block
          real(8), device :: T332(3,3,2), T6(6)
          T6(:) = Q(5,i,j,k-2:k+3) / (R * Q(1,i,j,k-2:k+3))
          mu(:) = mu6(T6(:))
          kTz   = heat_conduction6(mu(:), T6(:), dz(k))
        end block
        call tauxy4(mu(:), wx3(:), uz3(:), u6(:), tzx, utzx)
        call tauxy4(mu(:), vz3(:), wy3(:), v6(:), tzy, vtzy)
        call tauxx4(mu(:), wz3(:), ux3(:), vy3(:), w6(:), tzz, wtzz)
      end block
    else
      block
        real(8), dimension(2), device :: u2, v2, w2, mx, my
        real(8) mz, muz, mvz, mwz, mwx, mux, mvy, mwy
        block
          real(8), device :: T332(3,3,2)
          T332(:,:,:) = Q(5,i-1:i+1,j-1:j+1,k:k+1) / (R * Q(1,i-1:i+1,j-1:j+1,k:k+1))
          block
            real(8), device :: Tx(3,2)
            Tx(:,:) = T332(:,2,:)
            mx(:)   = mu32(Tx(:,:))
          end block
          block
            real(8), device :: Ty(3,2)
            Ty(:,:) = T332(2,:,:)
            my(:)   = mu32(Ty(:,:))
          end block
          block
            real(8), device :: Tz(2)
            Tz(:)   = T332(2,2,:)
            mz      = mu2(Tz(:))
            kTz     = Cp * mz * (-Tz(1) + Tz(2)) * dz(k) / Pr
          end block
        end block
        block
          real(8), device :: u312(3,2)
          u312(:,:) = Q(2,i-1:i+1,j,k:k+1)
          mux       = dy32(mx(:), u312(:,:), dx(i))
          u2(:)     = u312(2,:)
          muz       = mz * (-u2(1) + u2(2)) * dz(k)
        end block
        block
          real(8), device :: v132(3,2)
          v132(:,:) = Q(3,i,j-1:j+1,k:k+1)
          mvy       = dy32(my(:), v132(:,:), dy(j))
          v2(:)     = v132(2,:)
          mvz       = mz * (-v2(1) + v2(2)) * dz(k)
        end block
        block
          real(8), device :: w312(3,2)
          w312(:,:) = Q(4,i-1:i+1,j,k:k+1)
          mwx       = dy32(mx(:), w312(:,:), dx(i))
        end block
        block
          real(8), device :: w132(3,2)
          w132(:,:) = Q(4,i,j-1:j+1,k:k+1)
          mwy       = dy32(my(:), w132(:,:), dy(j))
          w2(:)     = w132(2,:)
          mwz       = mz * (-w2(1) + w2(2)) * dz(k)
        end block
        tzx  = mwx + muz
        tzy  = mvz + mwy
        tzz  = 2.d0 * (2.d0 * mwz - mux - mvy) / 3.d0
        utzx = 0.5d0 * (u2(1) + u2(2)) * tzx
        vtzy = 0.5d0 * (v2(1) + v2(2)) * tzy
        wtzz = 0.5d0 * (w2(1) + w2(2)) * tzz
      end block
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
    integer i, j, k
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, Hsgs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      block
        real(8), dimension(6), device :: u6, v6, w6
        real(8), dimension(3), device :: uz3, vz3, wz3, wx3, ux3, vy3, wy3, mu
        block
          real(8), device :: u516(5,6)
          u516(:,:) = Q(2,i-2:i+2,j,k-2:k+3)
          ux3(:)    = dy56(u516(:,:), dx(i))
          u6(:)     = u516(3,:)
          uz3(:)    = dx6(u6(:), dz(k))
        end block
        block
          real(8), device :: v156(5,6)
          v156(:,:) = Q(3,i,j-2:j+2,k-2:k+3)
          vy3(:)    = dy56(v156(:,:), dy(j))
          v6(:)     = v156(3,:)
          vz3(:)    = dx6(v6(:), dz(k))
        end block
        block
          real(8), device :: w516(5,6)
          w516(:,:) = Q(4,i-2:i+2,j,k-2:k+3)
          wx3(:)    = dy56(w516(:,:), dx(i))
        end block
        block
          real(8), device :: w156(5,6)
          w156(:,:) = Q(4,i,j-2:j+2,k-2:k+3)
          wy3(:)    = dy56(w156(:,:), dy(j))
          w6(:)     = w156(3,:)
          wz3(:)    = dx6(w6(:), dz(k))
        end block
        block
          real(8), device :: T332(3,3,2), T6(6)
          T6(:) = Q(5,i,j,k-2:k+3) / (R * Q(1,i,j,k-2:k+3))
          mu(:) = mu6(T6(:))
          kTz   = heat_conduction6(mu(:), T6(:), dz(k))
        end block
        call tauxy4(mu(:), wx3(:), uz3(:), u6(:), tzx, utzx)
        call tauxy4(mu(:), vz3(:), wy3(:), v6(:), tzy, vtzy)
        call tauxx4(mu(:), wz3(:), ux3(:), vy3(:), w6(:), tzz, wtzz)
        block
          real(8), device :: H(4)
          real(8) mutz
          mutz = 0.0625d0 * (-mut(i,j,k-1) + 9.d0 * (mut(i,j,k) + mut(i,j,k+1)) -mut(i,j,k+2))
          tzx  = tzx + mutz * (wx3(2) + uz3(2))
          tzy  = tzy + mutz * (vz3(2) + wy3(2))
          tzz  = tzz + 2.d0 * mutz * (2.d0 * wz3(2) - ux3(2) - vy3(2)) / 3.d0
          H(:) = (gamma * Q(5,i,j,k-1:k+2) / (Q(1,i,j,k-1:k+2) * gamma_1)) &
                 + 0.5d0 * (Q(2,i,j,k-1:k+2)**2 + Q(3,i,j,k-1:k+2)**2 + Q(4,i,j,k-1:k+2)**2) + qc2(i,j,k-1:k+2)
          Hsgs = -mutz * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) / 3.d0) * dz(k) / Prt
        end block
      end block
    else
      block
        real(8), dimension(2), device :: u2, v2, w2, mx, mxsgs, my, mysgs
        real(8) mz, mzsgs, muz, muzsgs, mvz, mvzsgs, mwz, mwzsgs, mwx, mwxsgs, mux, muxsgs, mvy, mvysgs, mwy, mwysgs
        ! SGS
        mz    = 0.5d0 * (mut(i,j,k) + mut(i,j,k+1))
        mx(:) = (/0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j,k+1) + mut(i,j,k+1)), &
                  0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j,k+1) + mut(i+1,j,k+1))/)
        my(:) = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i,j-1,k+1) + mut(i,j,k+1)), &
                  0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i,j,k+1) + mut(i,j+1,k+1))/)
        block
          real(8), device :: T332(3,3,2)
          T332(:,:,:) = Q(5,i-1:i+1,j-1:j+1,k:k+1) / (R * Q(1,i-1:i+1,j-1:j+1,k:k+1))
          block
            real(8), device :: Tx(3,2)
            Tx(:,:) = T332(:,2,:)
            mx(:)   = mu32(Tx(:,:))
          end block
          block
            real(8), device :: Ty(3,2)
            Ty(:,:) = T332(2,:,:)
            my(:)   = mu32(Ty(:,:))
          end block
          block
            real(8), device :: Tz(2)
            Tz(:)   = T332(2,2,:)
            mz      = mu2(Tz(:))
            kTz     = Cp * mz * (-Tz(1) + Tz(2)) * dz(k) / Pr
          end block
        end block
        block
          real(8), device :: u312(3,2)
          u312(:,:) = Q(2,i-1:i+1,j,k:k+1)
          mux       = dy32(mx(:), u312(:,:), dx(i))
          muxsgs    = dy32(mxsgs(:), u312(:,:), dx(i))
          u2(:)     = u312(2,:)
          muz       = mz * (-u2(1) + u2(2)) * dz(k)
          muzsgs    = mzsgs * (-u2(1) + u2(2)) * dz(k)
        end block
        block
          real(8), device :: v132(3,2)
          v132(:,:) = Q(3,i,j-1:j+1,k:k+1)
          mvy       = dy32(my(:), v132(:,:), dy(j))
          mvysgs    = dy32(mysgs(:), v132(:,:), dy(j))
          v2(:)     = v132(2,:)
          mvz       = mz * (-v2(1) + v2(2)) * dz(k)
          mvzsgs    = mzsgs * (-v2(1) + v2(2)) * dz(k)
        end block
        block
          real(8), device :: w312(3,2)
          w312(:,:) = Q(4,i-1:i+1,j,k:k+1)
          mwx       = dy32(mx(:), w312(:,:), dx(i))
          mwxsgs    = dy32(mxsgs(:), w312(:,:), dx(i))
        end block
        block
          real(8), device :: w132(3,2)
          w132(:,:) = Q(4,i,j-1:j+1,k:k+1)
          mwy       = dy32(my(:), w132(:,:), dy(j))
          mwysgs    = dy32(mysgs(:), w132(:,:), dy(j))
          w2(:)     = w132(2,:)
          mwz       = mz * (-w2(1) + w2(2)) * dz(k)
          mwzsgs    = mzsgs * (-w2(1) + w2(2)) * dz(k)
        end block
        tzx  = mwx + muz
        tzy  = mvz + mwy
        tzz  = 2.d0 * (2.d0 * mwz - mux - mvy) / 3.d0
        utzx = 0.5d0 * (u2(1) + u2(2)) * tzx
        vtzy = 0.5d0 * (v2(1) + v2(2)) * tzy
        wtzz = 0.5d0 * (w2(1) + w2(2)) * tzz
        tzx  = tzx + mwx + muz
        tzy  = tzy + mvz + mwy
        tzz  = tzz + 2.d0 * (2.d0 * mwz - mux - mvy) / 3.d0
        block
          real(8), device :: H(2)
          H(:) = (gamma * Q(5,i,j,k:k+1) / (Q(1,i,j,k:k+1) * gamma_1)) &
                 + 0.5d0 * (Q(2,i,j,k:k+1)**2 + Q(3,i,j,k:k+1)**2 + Q(4,i,j,k:k+1)**2) + qc2(i,j,k:k+1)
          Hsgs = -mz * (-H(1) + H(2)) * dz(k) / Prt
        end block
      end block
    endif
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz + Hsgs)
  end subroutine calc_Gv_LES
end module calc_visc

