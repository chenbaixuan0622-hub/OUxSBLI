module calc_visc4
  use mod_globals, only : id_visc, gamma, R, Pr, Prt, dt
  use mod_constant, only : Cp, gamma_1, Cp_over_Pr, one_third
  use calc_visc_common
  use calc_rand
  implicit none
contains
  attributes(global) subroutine calc_Ev4(nx, ny, nz, dx, dy, dz, Q, E, seed)
    use calc_sutherland, only : mu6, mu2, mu23
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(inout), device :: E(5,nx-1,ny-2,nz-2)
    integer(8), intent(inout), device, optional :: seed(nx,ny,nz)
    integer i, j, k
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8), device :: mu(3)
        block ! dQdx
          real(8), device :: T6(6)
          T6(:) = Q(5,i-2:i+3,j,k) / (R * Q(1,i-2:i+3,j,k))
          mu(:) = mu6(T6(:))
          kTx   = heat_conduction6(mu(:), T6(:), dx(i))
        end block


        block ! txx
          real(8), device :: ux3(3), vy3(3), wz3(3)
          ux3(:) = 0.125d0 * (9.d0 * (-Q(2,i-1:i+1,j,k) + Q(2,i:i+2,j,k)) &
                                   - (-Q(2,i-2:i,j,k)   + Q(2,i+1:i+3,j,k)) * one_third) * dx(i)
          block
            real(8) vy_1, vy_2, vy_3, vy_4, vy_5, vy_6
            vy_1 = dy5(Q(3,i-2,j-2,k), Q(3,i-2,j-1,k), Q(3,i-2,j+1,k), Q(3,i-2,j+2,k), dy(j))
            vy_2 = dy5(Q(3,i-1,j-2,k), Q(3,i-1,j-1,k), Q(3,i-1,j+1,k), Q(3,i-1,j+2,k), dy(j))
            vy_3 = dy5(Q(3,i  ,j-2,k), Q(3,i  ,j-1,k), Q(3,i  ,j+1,k), Q(3,i  ,j+2,k), dy(j))
            vy_4 = dy5(Q(3,i+1,j-2,k), Q(3,i+1,j-1,k), Q(3,i+1,j+1,k), Q(3,i+1,j+2,k), dy(j))
            vy_5 = dy5(Q(3,i+2,j-2,k), Q(3,i+2,j-1,k), Q(3,i+2,j+1,k), Q(3,i+2,j+2,k), dy(j))
            vy_6 = dy5(Q(3,i+3,j-2,k), Q(3,i+3,j-1,k), Q(3,i+3,j+1,k), Q(3,i+3,j+2,k), dy(j))
            vy3(:) = interpolation6_scalar(vy_1, vy_2, vy_3, vy_4, vy_5, vy_6)
          end block
          block
            real(8) wz_1, wz_2, wz_3, wz_4, wz_5, wz_6
            wz_1 = dy5(Q(4,i-2,j,k-2), Q(4,i-2,j,k-1), Q(4,i-2,j,k+1), Q(4,i-2,j,k+2), dz(k))
            wz_2 = dy5(Q(4,i-1,j,k-2), Q(4,i-1,j,k-1), Q(4,i-1,j,k+1), Q(4,i-1,j,k+2), dz(k))
            wz_3 = dy5(Q(4,i  ,j,k-2), Q(4,i  ,j,k-1), Q(4,i  ,j,k+1), Q(4,i  ,j,k+2), dz(k))
            wz_4 = dy5(Q(4,i+1,j,k-2), Q(4,i+1,j,k-1), Q(4,i+1,j,k+1), Q(4,i+1,j,k+2), dz(k))
            wz_5 = dy5(Q(4,i+2,j,k-2), Q(4,i+2,j,k-1), Q(4,i+2,j,k+1), Q(4,i+2,j,k+2), dz(k))
            wz_6 = dy5(Q(4,i+3,j,k-2), Q(4,i+3,j,k-1), Q(4,i+3,j,k+1), Q(4,i+3,j,k+2), dz(k))
            wz3(:) = interpolation6_scalar(wz_1, wz_2, wz_3, wz_4, wz_5, wz_6)
          end block
          call tauxx_4(mu(:), ux3(:), vy3(:), wz3(:), &
                      Q(2,i-2,j,k), Q(2,i-1,j,k), Q(2,i,j,k), Q(2,i+1,j,k), Q(2,i+2,j,k), Q(2,i+3,j,k), txx, utxx)
        end block


        block ! txy
          real(8), device :: vx3(3), uy3(3)
          vx3(:) = 0.125d0 * (9.d0 * (-Q(3,i-1:i+1,j,k) + Q(3,i:i+2,j,k)) &
                                   - (-Q(3,i-2:i,j,k)   + Q(3,i+1:i+3,j,k)) * one_third) * dx(i)
          block
            real(8) uy_1, uy_2, uy_3, uy_4, uy_5, uy_6
            uy_1 = dy5(Q(2,i-2,j-2,k), Q(2,i-2,j-1,k), Q(2,i-2,j+1,k), Q(2,i-2,j+2,k), dy(j))
            uy_2 = dy5(Q(2,i-1,j-2,k), Q(2,i-1,j-1,k), Q(2,i-1,j+1,k), Q(2,i-1,j+2,k), dy(j))
            uy_3 = dy5(Q(2,i  ,j-2,k), Q(2,i  ,j-1,k), Q(2,i  ,j+1,k), Q(2,i  ,j+2,k), dy(j))
            uy_4 = dy5(Q(2,i+1,j-2,k), Q(2,i+1,j-1,k), Q(2,i+1,j+1,k), Q(2,i+1,j+2,k), dy(j))
            uy_5 = dy5(Q(2,i+2,j-2,k), Q(2,i+2,j-1,k), Q(2,i+2,j+1,k), Q(2,i+2,j+2,k), dy(j))
            uy_6 = dy5(Q(2,i+3,j-2,k), Q(2,i+3,j-1,k), Q(2,i+3,j+1,k), Q(2,i+3,j+2,k), dy(j))
            uy3(:) = interpolation6_scalar(uy_1, uy_2, uy_3, uy_4, uy_5, uy_6)
          end block
          call tauxy_4(mu(:), uy3(:), vx3(:), &
                      Q(3,i-2,j,k), Q(3,i-1,j,k), Q(3,i,j,k), Q(3,i+1,j,k), Q(3,i+2,j,k), Q(3,i+3,j,k), txy, vtxy)
        end block


        block ! txz
          real(8), device :: wx3(3), uz3(3)
          wx3(:) = 0.125d0 * (9.d0 * (-Q(4,i-1:i+1,j,k) + Q(4,i:i+2,j,k)) &
                                   - (-Q(4,i-2:i,j,k)   + Q(4,i+1:i+3,j,k)) * one_third) * dx(i)
          block
            real(8) uz_1, uz_2, uz_3, uz_4, uz_5, uz_6
            uz_1 = dy5(Q(2,i-2,j,k-2), Q(2,i-2,j,k-1), Q(2,i-2,j,k+1), Q(2,i-2,j,k+2), dz(k))
            uz_2 = dy5(Q(2,i-1,j,k-2), Q(2,i-1,j,k-1), Q(2,i-1,j,k+1), Q(2,i-1,j,k+2), dz(k))
            uz_3 = dy5(Q(2,i  ,j,k-2), Q(2,i  ,j,k-1), Q(2,i  ,j,k+1), Q(2,i  ,j,k+2), dz(k))
            uz_4 = dy5(Q(2,i+1,j,k-2), Q(2,i+1,j,k-1), Q(2,i+1,j,k+1), Q(2,i+1,j,k+2), dz(k))
            uz_5 = dy5(Q(2,i+2,j,k-2), Q(2,i+2,j,k-1), Q(2,i+2,j,k+1), Q(2,i+2,j,k+2), dz(k))
            uz_6 = dy5(Q(2,i+3,j,k-2), Q(2,i+3,j,k-1), Q(2,i+3,j,k+1), Q(2,i+3,j,k+2), dz(k))
            uz3(:) = interpolation6_scalar(uz_1, uz_2, uz_3, uz_4, uz_5, uz_6)
          end block
          call tauxy_4(mu(:), wx3(:), uz3(:), &
                      Q(4,i-2,j,k), Q(4,i-1,j,k), Q(4,i,j,k), Q(4,i+1,j,k), Q(4,i+2,j,k), Q(4,i+3,j,k), txz, wtxz)
        end block
      end block
    else
      block
        real(8) mx, mux, mvx, mwx, muy, mvy, muz, mwz
        block
          real(8), device :: Tx(2)
          Tx(:) = Q(5,i:i+1,j,k) / (R * Q(1,i:i+1,j,k))
          mx    = mu2(Tx(:))
          kTx   = Cp_over_Pr * mx * (-Tx(1) + Tx(2)) * dx(i)
        end block
        block
          real(8), device :: Ty(2,3), my(2)
          Ty(:,:) = Q(5,i:i+1,j-1:j+1,k) / (R * Q(1,i:i+1,j-1:j+1,k))
          my(:)   = mu23(Ty(:,:))
          muy = 0.25d0 * (my(1) * (-Q(2,i,j-1,k) + Q(2,i,j,k) - Q(2,i+1,j-1,k) + Q(2,i+1,j,k)) &
                        + my(2) * (-Q(2,i,j,k) + Q(2,i,j+1,k) - Q(2,i+1,j,k) + Q(2,i+1,j+1,k))) * dy(j)
          mvy = 0.25d0 * (my(1) * (-Q(3,i,j-1,k) + Q(3,i,j,k) - Q(3,i+1,j-1,k) + Q(3,i+1,j,k)) &
                        + my(2) * (-Q(3,i,j,k) + Q(3,i,j+1,k) - Q(3,i+1,j,k) + Q(3,i+1,j+1,k))) * dy(j)
        end block
        block
          real(8), device :: Tz(2,3), mz(2)
          Tz(:,:) = Q(5,i:i+1,j,k-1:k+1) / (R * Q(1,i:i+1,j,k-1:k+1))
          mz(:)   = mu23(Tz(:,:))
          muz = 0.25d0 * (mz(1) * (-Q(2,i,j,k-1) + Q(2,i,j,k) - Q(2,i+1,j,k-1) + Q(2,i+1,j,k)) &
                        + mz(2) * (-Q(2,i,j,k) + Q(2,i,j,k+1) - Q(2,i+1,j,k) + Q(2,i+1,j,k+1))) * dz(k)
          mwz = 0.25d0 * (mz(1) * (-Q(4,i,j,k-1) + Q(4,i,j,k) - Q(4,i+1,j,k-1) + Q(4,i+1,j,k)) &
                        + mz(2) * (-Q(4,i,j,k) + Q(4,i,j,k+1) - Q(4,i+1,j,k) + Q(4,i+1,j,k+1))) * dz(k)
        end block
        mux = mx * (-Q(2,i,j,k) + Q(2,i+1,j,k)) * dx(i)
        mvx = mx * (-Q(3,i,j,k) + Q(3,i+1,j,k)) * dx(i)
        mwx = mx * (-Q(4,i,j,k) + Q(4,i+1,j,k)) * dx(i)
        txx  = 2.d0 * (2.d0 * mux - mvy - mwz) * one_third
        txy  = muy + mvx
        txz  = mwx + muz
        utxx = 0.5d0 * (Q(2,i,j,k) + Q(2,i+1,j,k)) * txx
        vtxy = 0.5d0 * (Q(3,i,j,k) + Q(3,i+1,j,k)) * txy
        wtxz = 0.5d0 * (Q(4,i,j,k) + Q(4,i+1,j,k)) * txz
      end block
    endif
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx)
  end subroutine calc_Ev4
 

  attributes(global) subroutine calc_Ev_LES4(nx, ny, nz, dx, dy, dz, Q, mut, qc2, E)
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
    if (3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8), device :: mu(3)
        real(8) mutx
        mutx = 0.0625d0 * (-mut(i-1,j,k) + 9.d0 * (mut(i,j,k) + mut(i+1,j,k)) -mut(i+2,j,k))
        block ! dQdx
          real(8), device :: T6(6)
          T6(:) = Q(5,i-2:i+3,j,k) / (R * Q(1,i-2:i+3,j,k))
          mu(:) = mu6(T6(:))
          kTx   = heat_conduction6(mu(:), T6(:), dx(i))
        end block


        block ! txx
          real(8), device :: ux3(3), vy3(3), wz3(3)
          ux3(:) = 0.125d0 * (9.d0 * (-Q(2,i-1:i+1,j,k) + Q(2,i:i+2,j,k)) &
                                   - (-Q(2,i-2:i,j,k)   + Q(2,i+1:i+3,j,k)) * one_third) * dx(i)
          block
            real(8) vy_1, vy_2, vy_3, vy_4, vy_5, vy_6
            vy_1 = dy5(Q(3,i-2,j-2,k), Q(3,i-2,j-1,k), Q(3,i-2,j+1,k), Q(3,i-2,j+2,k), dy(j))
            vy_2 = dy5(Q(3,i-1,j-2,k), Q(3,i-1,j-1,k), Q(3,i-1,j+1,k), Q(3,i-1,j+2,k), dy(j))
            vy_3 = dy5(Q(3,i  ,j-2,k), Q(3,i  ,j-1,k), Q(3,i  ,j+1,k), Q(3,i  ,j+2,k), dy(j))
            vy_4 = dy5(Q(3,i+1,j-2,k), Q(3,i+1,j-1,k), Q(3,i+1,j+1,k), Q(3,i+1,j+2,k), dy(j))
            vy_5 = dy5(Q(3,i+2,j-2,k), Q(3,i+2,j-1,k), Q(3,i+2,j+1,k), Q(3,i+2,j+2,k), dy(j))
            vy_6 = dy5(Q(3,i+3,j-2,k), Q(3,i+3,j-1,k), Q(3,i+3,j+1,k), Q(3,i+3,j+2,k), dy(j))
            vy3(:) = interpolation6_scalar(vy_1, vy_2, vy_3, vy_4, vy_5, vy_6)
          end block
          block
            real(8) wz_1, wz_2, wz_3, wz_4, wz_5, wz_6
            wz_1 = dy5(Q(4,i-2,j,k-2), Q(4,i-2,j,k-1), Q(4,i-2,j,k+1), Q(4,i-2,j,k+2), dz(k))
            wz_2 = dy5(Q(4,i-1,j,k-2), Q(4,i-1,j,k-1), Q(4,i-1,j,k+1), Q(4,i-1,j,k+2), dz(k))
            wz_3 = dy5(Q(4,i  ,j,k-2), Q(4,i  ,j,k-1), Q(4,i  ,j,k+1), Q(4,i  ,j,k+2), dz(k))
            wz_4 = dy5(Q(4,i+1,j,k-2), Q(4,i+1,j,k-1), Q(4,i+1,j,k+1), Q(4,i+1,j,k+2), dz(k))
            wz_5 = dy5(Q(4,i+2,j,k-2), Q(4,i+2,j,k-1), Q(4,i+2,j,k+1), Q(4,i+2,j,k+2), dz(k))
            wz_6 = dy5(Q(4,i+3,j,k-2), Q(4,i+3,j,k-1), Q(4,i+3,j,k+1), Q(4,i+3,j,k+2), dz(k))
            wz3(:) = interpolation6_scalar(wz_1, wz_2, wz_3, wz_4, wz_5, wz_6)
          end block
          call tauxx_4(mu(:), ux3(:), vy3(:), wz3(:), &
                      Q(2,i-2,j,k), Q(2,i-1,j,k), Q(2,i,j,k), Q(2,i+1,j,k), Q(2,i+2,j,k), Q(2,i+3,j,k), txx, utxx)
          txx  = txx + 2.d0 * mutx * (2.d0 * ux3(2) - vy3(2) - wz3(2)) * one_third
        end block


        block ! txy
          real(8), device :: vx3(3), uy3(3)
          vx3(:) = 0.125d0 * (9.d0 * (-Q(3,i-1:i+1,j,k) + Q(3,i:i+2,j,k)) &
                                   - (-Q(3,i-2:i,j,k)   + Q(3,i+1:i+3,j,k)) * one_third) * dx(i)
          block
            real(8) uy_1, uy_2, uy_3, uy_4, uy_5, uy_6
            uy_1 = dy5(Q(2,i-2,j-2,k), Q(2,i-2,j-1,k), Q(2,i-2,j+1,k), Q(2,i-2,j+2,k), dy(j))
            uy_2 = dy5(Q(2,i-1,j-2,k), Q(2,i-1,j-1,k), Q(2,i-1,j+1,k), Q(2,i-1,j+2,k), dy(j))
            uy_3 = dy5(Q(2,i  ,j-2,k), Q(2,i  ,j-1,k), Q(2,i  ,j+1,k), Q(2,i  ,j+2,k), dy(j))
            uy_4 = dy5(Q(2,i+1,j-2,k), Q(2,i+1,j-1,k), Q(2,i+1,j+1,k), Q(2,i+1,j+2,k), dy(j))
            uy_5 = dy5(Q(2,i+2,j-2,k), Q(2,i+2,j-1,k), Q(2,i+2,j+1,k), Q(2,i+2,j+2,k), dy(j))
            uy_6 = dy5(Q(2,i+3,j-2,k), Q(2,i+3,j-1,k), Q(2,i+3,j+1,k), Q(2,i+3,j+2,k), dy(j))
            uy3(:) = interpolation6_scalar(uy_1, uy_2, uy_3, uy_4, uy_5, uy_6)
          end block
          call tauxy_4(mu(:), uy3(:), vx3(:), &
                      Q(3,i-2,j,k), Q(3,i-1,j,k), Q(3,i,j,k), Q(3,i+1,j,k), Q(3,i+2,j,k), Q(3,i+3,j,k), txy, vtxy)
          txy  = txy + mutx * (uy3(2) + vx3(2))
        end block


        block ! txz
          real(8), device :: wx3(3), uz3(3)
          wx3(:) = 0.125d0 * (9.d0 * (-Q(4,i-1:i+1,j,k) + Q(4,i:i+2,j,k)) &
                                   - (-Q(4,i-2:i,j,k)   + Q(4,i+1:i+3,j,k)) * one_third) * dx(i)
          block
            real(8) uz_1, uz_2, uz_3, uz_4, uz_5, uz_6
            uz_1 = dy5(Q(2,i-2,j,k-2), Q(2,i-2,j,k-1), Q(2,i-2,j,k+1), Q(2,i-2,j,k+2), dz(k))
            uz_2 = dy5(Q(2,i-1,j,k-2), Q(2,i-1,j,k-1), Q(2,i-1,j,k+1), Q(2,i-1,j,k+2), dz(k))
            uz_3 = dy5(Q(2,i  ,j,k-2), Q(2,i  ,j,k-1), Q(2,i  ,j,k+1), Q(2,i  ,j,k+2), dz(k))
            uz_4 = dy5(Q(2,i+1,j,k-2), Q(2,i+1,j,k-1), Q(2,i+1,j,k+1), Q(2,i+1,j,k+2), dz(k))
            uz_5 = dy5(Q(2,i+2,j,k-2), Q(2,i+2,j,k-1), Q(2,i+2,j,k+1), Q(2,i+2,j,k+2), dz(k))
            uz_6 = dy5(Q(2,i+3,j,k-2), Q(2,i+3,j,k-1), Q(2,i+3,j,k+1), Q(2,i+3,j,k+2), dz(k))
            uz3(:) = interpolation6_scalar(uz_1, uz_2, uz_3, uz_4, uz_5, uz_6)
          end block
          call tauxy_4(mu(:), wx3(:), uz3(:), &
                      Q(4,i-2,j,k), Q(4,i-1,j,k), Q(4,i,j,k), Q(4,i+1,j,k), Q(4,i+2,j,k), Q(4,i+3,j,k), txz, wtxz)
          txz  = txz + mutx * (wx3(2) + uz3(2))
        end block


        block
          real(8) :: H(4)
          H(:) = (gamma * Q(5,i-1:i+2,j,k) / (Q(1,i-1:i+2,j,k) * gamma_1)) &
                 + 0.5d0 * (Q(2,i-1:i+2,j,k)**2 + Q(3,i-1:i+2,j,k)**2 + Q(4,i-1:i+2,j,k)**2) + qc2(i-1:i+2,j,k)
          Hsgs = -mutx * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) * one_third) * dx(i) / Prt
        end block
      end block
    else
      block
        real(8), dimension(2), device :: my, mysgs, mz, mzsgs
        real(8) mx, mxsgs, mux, muxsgs, mvx, mvxsgs, mwx, mwxsgs, muy, muysgs, mvy, mvysgs, muz, muzsgs, mwz, mwzsgs
        ! SGS
        mx    = 0.5d0 * (mut(i,j,k) + mut(i+1,j,k))
        my(:) = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i+1,j-1,k) + mut(i+1,j,k)), &
                  0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i+1,j,k) + mut(i+1,j+1,k))/)
        mz(:) = (/0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i+1,j,k-1) + mut(i+1,j,k)), &
                  0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i+1,j,k) + mut(i+1,j,k+1))/)
        block
          real(8), device :: Tx(2)
          Tx(:) = Q(5,i:i+1,j,k) / (R * Q(1,i:i+1,j,k))
          mx    = mu2(Tx(:))
          kTx   = Cp_over_Pr * mx * (-Tx(1) + Tx(2)) * dx(i)
        end block
        block
          real(8), device :: Ty(2,3)
          Ty(:,:) = Q(5,i:i+1,j-1:j+1,k) / (R * Q(1,i:i+1,j-1:j+1,k))
          my(:)   = mu23(Ty(:,:))
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
          real(8), device :: Tz(2,3)
          Tz(:,:) = Q(5,i:i+1,j,k-1:k+1) / (R * Q(1,i:i+1,j,k-1:k+1))
          mz(:)   = mu23(Tz(:,:))
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
        txx  = 2.d0 * (2.d0 * mux - mvy - mwz) * one_third
        txy  = muy + mvx
        txz  = mwx + muz
        utxx = 0.5d0 * (Q(2,i,j,k) + Q(2,i+1,j,k)) * txx
        vtxy = 0.5d0 * (Q(3,i,j,k) + Q(3,i+1,j,k)) * txy
        wtxz = 0.5d0 * (Q(4,i,j,k) + Q(4,i+1,j,k)) * txz
        txx = txx + 2.d0 * (2.d0 * muxsgs - mvysgs - mwzsgs) * one_third
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
  end subroutine calc_Ev_LES4
 

  attributes(global) subroutine calc_Fv4(nx, ny, nz, dy, dx, dz, Q, F, seed)
    use calc_sutherland, only : mu6, mu2, mu23, mu32
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(inout), device :: F(5,nx-2,ny-1,nz-2)
    integer(8), intent(inout), device, optional :: seed(nx,ny,nz)
    integer i, j, k
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8), device :: mu(3)
        block ! dQdy
          real(8), device :: T6(6)
          T6(:) = Q(5,i,j-2:j+3,k) / (R * Q(1,i,j-2:j+3,k))
          mu(:) = mu6(T6(:))
          kTy   = heat_conduction6(mu(:), T6(:), dy(j))
        end block
        

        block ! tyx
          real(8), device :: uy3(3), vx3(3)
          uy3(:) = 0.125d0 * (9.d0 * (-Q(2,i,j-1:j+1,k) + Q(2,i,j:j+1,k)) &
                                   - (-Q(2,i,j-2:j,k)   + Q(2,i,j+1:j+3,k)) * one_third) * dy(j)
          block
            real(8) vx_1, vx_2, vx_3, vx_4, vx_5, vx_6
            vx_1 = dy5(Q(3,i-2,j-2,k), Q(3,i-1,j-2,k), Q(3,i+1,j-2,k), Q(3,i+2,j-2,k), dx(i))
            vx_2 = dy5(Q(3,i-2,j-1,k), Q(3,i-1,j-1,k), Q(3,i+1,j-1,k), Q(3,i+2,j-1,k), dx(i))
            vx_3 = dy5(Q(3,i-2,j  ,k), Q(3,i-1,j  ,k), Q(3,i+1,j  ,k), Q(3,i+2,j  ,k), dx(i))
            vx_4 = dy5(Q(3,i-2,j+1,k), Q(3,i-1,j+1,k), Q(3,i+1,j+1,k), Q(3,i+2,j+1,k), dx(i))
            vx_5 = dy5(Q(3,i-2,j+2,k), Q(3,i-1,j+2,k), Q(3,i+1,j+2,k), Q(3,i+2,j+2,k), dx(i))
            vx_6 = dy5(Q(3,i-2,j+3,k), Q(3,i-1,j+3,k), Q(3,i+1,j+3,k), Q(3,i+2,j+3,k), dx(i))
            vx3(:) = interpolation6_scalar(vx_1, vx_2, vx_3, vx_4, vx_5, vx_6)
          end block
          call tauxy_4(mu(:), uy3(:), vx3(:), &
                      Q(2,i,j-2,k), Q(2,i,j-1,k), Q(2,i,j,k), Q(2,i,j+1,k), Q(2,i,j+2,k), Q(2,i,j+3,k), tyx, utyx)
        end block


        block ! tyy
          real(8), device :: vy3(3), wz3(3), ux3(3)
          vy3(:) = 0.125d0 * (9.d0 * (-Q(3,i,j-1:j+1,k) + Q(3,i,j:j+1,k)) &
                                   - (-Q(3,i,j-2:j,k)   + Q(3,i,j+1:j+3,k)) * one_third) * dy(j)
          block
            real(8) wz_1, wz_2, wz_3, wz_4, wz_5, wz_6
            wz_1 = dy5(Q(4,i,j-2,k-2), Q(4,i,j-2,k-1), Q(4,i,j-2,k+1), Q(4,i,j-2,k+2), dz(k))
            wz_2 = dy5(Q(4,i,j-1,k-2), Q(4,i,j-1,k-1), Q(4,i,j-1,k+1), Q(4,i,j-1,k+2), dz(k))
            wz_3 = dy5(Q(4,i,j  ,k-2), Q(4,i,j  ,k-1), Q(4,i,j  ,k+1), Q(4,i,j  ,k+2), dz(k))
            wz_4 = dy5(Q(4,i,j+1,k-2), Q(4,i,j+1,k-1), Q(4,i,j+1,k+1), Q(4,i,j+1,k+2), dz(k))
            wz_5 = dy5(Q(4,i,j+2,k-2), Q(4,i,j+2,k-1), Q(4,i,j+2,k+1), Q(4,i,j+2,k+2), dz(k))
            wz_6 = dy5(Q(4,i,j+3,k-2), Q(4,i,j+3,k-1), Q(4,i,j+3,k+1), Q(4,i,j+3,k+2), dz(k))
            wz3(:) = interpolation6_scalar(wz_1, wz_2, wz_3, wz_4, wz_5, wz_6)
          end block
          block
            real(8) ux_1, ux_2, ux_3, ux_4, ux_5, ux_6
            ux_1 = dy5(Q(2,i-2,j-2,k), Q(2,i-1,j-2,k), Q(2,i+1,j-2,k), Q(2,i+2,j-2,k), dx(i))
            ux_2 = dy5(Q(2,i-2,j-1,k), Q(2,i-1,j-1,k), Q(2,i+1,j-1,k), Q(2,i+2,j-1,k), dx(i))
            ux_3 = dy5(Q(2,i-2,j  ,k), Q(2,i-1,j  ,k), Q(2,i+1,j  ,k), Q(2,i+2,j  ,k), dx(i))
            ux_4 = dy5(Q(2,i-2,j+1,k), Q(2,i-1,j+1,k), Q(2,i+1,j+1,k), Q(2,i+2,j+1,k), dx(i))
            ux_5 = dy5(Q(2,i-2,j+2,k), Q(2,i-1,j+2,k), Q(2,i+1,j+2,k), Q(2,i+2,j+2,k), dx(i))
            ux_6 = dy5(Q(2,i-2,j+3,k), Q(2,i-1,j+3,k), Q(2,i+1,j+3,k), Q(2,i+2,j+3,k), dx(i))
            ux3(:) = interpolation6_scalar(ux_1, ux_2, ux_3, ux_4, ux_5, ux_6)
          end block
          call tauxx_4(mu(:), vy3(:), wz3(:), ux3(:), &
                      Q(3,i,j-2,k), Q(3,i,j-1,k), Q(3,i,j,k), Q(3,i,j+1,k), Q(3,i,j+2,k), Q(3,i,j+3,k), tyy, vtyy)
        end block


        block ! tyz
          real(8), device :: wy3(3), vz3(3)
          wy3(:) = 0.125d0 * (9.d0 * (-Q(4,i,j-1:j+1,k) + Q(4,i,j:j+1,k)) &
                                   - (-Q(4,i,j-2:j,k)   + Q(4,i,j+1:j+3,k)) * one_third) * dy(j)
          block
            real(8) vz_1, vz_2, vz_3, vz_4, vz_5, vz_6
            vz_1 = dy5(Q(3,i,j-2,k-2), Q(3,i,j-2,k-1), Q(3,i,j-2,k+1), Q(3,i,j-2,k+2), dz(k))
            vz_2 = dy5(Q(3,i,j-1,k-2), Q(3,i,j-1,k-1), Q(3,i,j-1,k+1), Q(3,i,j-1,k+2), dz(k))
            vz_3 = dy5(Q(3,i,j,  k-2), Q(3,i,j  ,k-1), Q(3,i,j  ,k+1), Q(3,i,j  ,k+2), dz(k))
            vz_4 = dy5(Q(3,i,j+1,k-2), Q(3,i,j+1,k-1), Q(3,i,j+1,k+1), Q(3,i,j+1,k+2), dz(k))
            vz_5 = dy5(Q(3,i,j+2,k-2), Q(3,i,j+2,k-1), Q(3,i,j+2,k+1), Q(3,i,j+2,k+2), dz(k))
            vz_6 = dy5(Q(3,i,j+3,k-2), Q(3,i,j+3,k-1), Q(3,i,j+3,k+1), Q(3,i,j+3,k+2), dz(k))
            vz3(:) = interpolation6_scalar(vz_1, vz_2, vz_3, vz_4, vz_5, vz_6)
          end block
          call tauxy_4(mu(:), vz3(:), wy3(:), &
                      Q(4,i,j-2,k), Q(4,i,j-1,k), Q(4,i,j,k), Q(4,i,j+1,k), Q(4,i,j+2,k), Q(4,i,j+3,k), tyz, wtyz)
        end block
      end block
    else
      block
        real(8) my, muy, mvy, mwy, mvz, mwz, mux, mvx
        block
          real(8), device :: Tx(3,2), mx(2)
          Tx(:,:) = Q(5,i-1:i+1,j:j+1,k) / (R * Q(1,i-1:i+1,j:j+1,k))
          mx(:)   = mu23(Tx(:,:))
          mux = 0.25d0 * (mx(1) * (-Q(2,i-1,j,k) + Q(2,i,j,k) - Q(2,i-1,j+1,k) + Q(2,i,j+1,k)) &
                        + mx(2) * (-Q(2,i,j,k) + Q(2,i+1,j,k) - Q(2,i,j+1,k) + Q(2,i+1,j+1,k))) * dx(i)
          mvx = 0.25d0 * (mx(1) * (-Q(3,i-1,j,k) + Q(3,i,j,k) - Q(3,i-1,j+1,k) + Q(3,i,j+1,k)) &
                        + mx(2) * (-Q(3,i,j,k) + Q(3,i+1,j,k) - Q(3,i,j+1,k) + Q(3,i+1,j+1,k))) * dx(i)
        end block
        block
          real(8), device :: Ty(2)
          Ty(:) = Q(5,i,j:j+1,k) / (R * Q(1,i,j:j+1,k))
          my    = mu2(Ty(:))
          kTy   = Cp_over_Pr * my * (-Ty(1) + Ty(2)) * dy(j)
        end block
        block
          real(8), device :: Tz(2,3), mz(2)
          Tz(:,:) = Q(5,i,j:j+1,k-1:k+1) / (R * Q(1,i,j:j+1,k-1:k+1))
          mz(:)   = mu32(Tz(:,:))
          mvz = 0.25d0 * (mz(1) * (-Q(3,i,j,k-1) + Q(3,i,j,k) - Q(3,i,j+1,k-1) + Q(3,i,j+1,k)) &
                        + mz(2) * (-Q(3,i,j,k) + Q(3,i,j,k+1) - Q(3,i,j+1,k) + Q(3,i,j+1,k+1))) * dz(k)
          mwz = 0.25d0 * (mz(1) * (-Q(4,i,j,k-1) + Q(4,i,j,k) - Q(4,i,j+1,k-1) + Q(4,i,j+1,k)) &
                        + mz(2) * (-Q(4,i,j,k) + Q(4,i,j,k+1) - Q(4,i,j+1,k) + Q(4,i,j+1,k+1))) * dz(k)
        end block
        muy = my * (-Q(2,i,j,k) + Q(2,i,j+1,k)) * dy(j)
        mvy = my * (-Q(3,i,j,k) + Q(3,i,j+1,k)) * dy(j)
        mwy = my * (-Q(4,i,j,k) + Q(4,i,j+1,k)) * dy(j)
        tyx  = muy + mvx
        tyy  = 2.d0 * (2.d0 * mvy - mwz - mux) * one_third
        tyz  = mvz + mwy
        utyx = 0.5d0 * (Q(2,i,j,k) + Q(2,i,j+1,k)) * tyx
        vtyy = 0.5d0 * (Q(3,i,j,k) + Q(3,i,j+1,k)) * tyy
        wtyz = 0.5d0 * (Q(4,i,j,k) + Q(4,i,j+1,k)) * tyz
      end block
    endif
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy)
  end subroutine calc_Fv4
 

  attributes(global) subroutine calc_Fv_LES4(nx, ny, nz, dy, dx, dz, Q, mut, qc2, F)
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
    if (3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8), device :: mu(3)
        real(8) muty
        muty = 0.0625d0 * (-mut(i,j-1,k) + 9.d0 * (mut(i,j,k) + mut(i,j+1,k)) -mut(i,j+2,k))
        block ! dQdy
          real(8), device :: T6(6)
          T6(:) = Q(5,i,j-2:j+3,k) / (R * Q(1,i,j-2:j+3,k))
          mu(:) = mu6(T6(:))
          kTy   = heat_conduction6(mu(:), T6(:), dy(j))
        end block
        

        block ! tyx
          real(8), device :: uy3(3), vx3(3)
          uy3(:) = 0.125d0 * (9.d0 * (-Q(2,i,j-1:j+1,k) + Q(2,i,j:j+1,k)) &
                                   - (-Q(2,i,j-2:j,k)   + Q(2,i,j+1:j+3,k)) * one_third) * dy(j)
          block
            real(8) vx_1, vx_2, vx_3, vx_4, vx_5, vx_6
            vx_1 = dy5(Q(3,i-2,j-2,k), Q(3,i-1,j-2,k), Q(3,i+1,j-2,k), Q(3,i+2,j-2,k), dx(i))
            vx_2 = dy5(Q(3,i-2,j-1,k), Q(3,i-1,j-1,k), Q(3,i+1,j-1,k), Q(3,i+2,j-1,k), dx(i))
            vx_3 = dy5(Q(3,i-2,j  ,k), Q(3,i-1,j  ,k), Q(3,i+1,j  ,k), Q(3,i+2,j  ,k), dx(i))
            vx_4 = dy5(Q(3,i-2,j+1,k), Q(3,i-1,j+1,k), Q(3,i+1,j+1,k), Q(3,i+2,j+1,k), dx(i))
            vx_5 = dy5(Q(3,i-2,j+2,k), Q(3,i-1,j+2,k), Q(3,i+1,j+2,k), Q(3,i+2,j+2,k), dx(i))
            vx_6 = dy5(Q(3,i-2,j+3,k), Q(3,i-1,j+3,k), Q(3,i+1,j+3,k), Q(3,i+2,j+3,k), dx(i))
            vx3(:) = interpolation6_scalar(vx_1, vx_2, vx_3, vx_4, vx_5, vx_6)
          end block
          call tauxy_4(mu(:), uy3(:), vx3(:), &
                      Q(2,i,j-2,k), Q(2,i,j-1,k), Q(2,i,j,k), Q(2,i,j+1,k), Q(2,i,j+2,k), Q(2,i,j+3,k), tyx, utyx)
          tyx = tyx + muty * (uy3(2) + vx3(2))
        end block


        block ! tyy
          real(8), device :: vy3(3), wz3(3), ux3(3)
          vy3(:) = 0.125d0 * (9.d0 * (-Q(3,i,j-1:j+1,k) + Q(3,i,j:j+1,k)) &
                                   - (-Q(3,i,j-2:j,k)   + Q(3,i,j+1:j+3,k)) * one_third) * dy(j)
          block
            real(8) wz_1, wz_2, wz_3, wz_4, wz_5, wz_6
            wz_1 = dy5(Q(4,i,j-2,k-2), Q(4,i,j-2,k-1), Q(4,i,j-2,k+1), Q(4,i,j-2,k+2), dz(k))
            wz_2 = dy5(Q(4,i,j-1,k-2), Q(4,i,j-1,k-1), Q(4,i,j-1,k+1), Q(4,i,j-1,k+2), dz(k))
            wz_3 = dy5(Q(4,i,j  ,k-2), Q(4,i,j  ,k-1), Q(4,i,j  ,k+1), Q(4,i,j  ,k+2), dz(k))
            wz_4 = dy5(Q(4,i,j+1,k-2), Q(4,i,j+1,k-1), Q(4,i,j+1,k+1), Q(4,i,j+1,k+2), dz(k))
            wz_5 = dy5(Q(4,i,j+2,k-2), Q(4,i,j+2,k-1), Q(4,i,j+2,k+1), Q(4,i,j+2,k+2), dz(k))
            wz_6 = dy5(Q(4,i,j+3,k-2), Q(4,i,j+3,k-1), Q(4,i,j+3,k+1), Q(4,i,j+3,k+2), dz(k))
            wz3(:) = interpolation6_scalar(wz_1, wz_2, wz_3, wz_4, wz_5, wz_6)
          end block
          block
            real(8) ux_1, ux_2, ux_3, ux_4, ux_5, ux_6
            ux_1 = dy5(Q(2,i-2,j-2,k), Q(2,i-1,j-2,k), Q(2,i+1,j-2,k), Q(2,i+2,j-2,k), dx(i))
            ux_2 = dy5(Q(2,i-2,j-1,k), Q(2,i-1,j-1,k), Q(2,i+1,j-1,k), Q(2,i+2,j-1,k), dx(i))
            ux_3 = dy5(Q(2,i-2,j  ,k), Q(2,i-1,j  ,k), Q(2,i+1,j  ,k), Q(2,i+2,j  ,k), dx(i))
            ux_4 = dy5(Q(2,i-2,j+1,k), Q(2,i-1,j+1,k), Q(2,i+1,j+1,k), Q(2,i+2,j+1,k), dx(i))
            ux_5 = dy5(Q(2,i-2,j+2,k), Q(2,i-1,j+2,k), Q(2,i+1,j+2,k), Q(2,i+2,j+2,k), dx(i))
            ux_6 = dy5(Q(2,i-2,j+3,k), Q(2,i-1,j+3,k), Q(2,i+1,j+3,k), Q(2,i+2,j+3,k), dx(i))
            ux3(:) = interpolation6_scalar(ux_1, ux_2, ux_3, ux_4, ux_5, ux_6)
          end block
          call tauxx_4(mu(:), vy3(:), wz3(:), ux3(:), &
                      Q(3,i,j-2,k), Q(3,i,j-1,k), Q(3,i,j,k), Q(3,i,j+1,k), Q(3,i,j+2,k), Q(3,i,j+3,k), tyy, vtyy)
          tyy = tyy + 2.d0 * muty * (2.d0 * vy3(2) - ux3(2) - wz3(2)) * one_third
        end block


        block ! tyz
          real(8), device :: wy3(3), vz3(3)
          wy3(:) = 0.125d0 * (9.d0 * (-Q(4,i,j-1:j+1,k) + Q(4,i,j:j+1,k)) &
                                   - (-Q(4,i,j-2:j,k)   + Q(4,i,j+1:j+3,k)) * one_third) * dy(j)
          block
            real(8) vz_1, vz_2, vz_3, vz_4, vz_5, vz_6
            vz_1 = dy5(Q(3,i,j-2,k-2), Q(3,i,j-2,k-1), Q(3,i,j-2,k+1), Q(3,i,j-2,k+2), dz(k))
            vz_2 = dy5(Q(3,i,j-1,k-2), Q(3,i,j-1,k-1), Q(3,i,j-1,k+1), Q(3,i,j-1,k+2), dz(k))
            vz_3 = dy5(Q(3,i,j,  k-2), Q(3,i,j  ,k-1), Q(3,i,j  ,k+1), Q(3,i,j  ,k+2), dz(k))
            vz_4 = dy5(Q(3,i,j+1,k-2), Q(3,i,j+1,k-1), Q(3,i,j+1,k+1), Q(3,i,j+1,k+2), dz(k))
            vz_5 = dy5(Q(3,i,j+2,k-2), Q(3,i,j+2,k-1), Q(3,i,j+2,k+1), Q(3,i,j+2,k+2), dz(k))
            vz_6 = dy5(Q(3,i,j+3,k-2), Q(3,i,j+3,k-1), Q(3,i,j+3,k+1), Q(3,i,j+3,k+2), dz(k))
            vz3(:) = interpolation6_scalar(vz_1, vz_2, vz_3, vz_4, vz_5, vz_6)
          end block
          call tauxy_4(mu(:), vz3(:), wy3(:), &
                      Q(4,i,j-2,k), Q(4,i,j-1,k), Q(4,i,j,k), Q(4,i,j+1,k), Q(4,i,j+2,k), Q(4,i,j+3,k), tyz, wtyz)
          tyz = tyz + muty * (vz3(2) + wy3(2))
        end block


        block
          real(8), device :: H(4)
          H(:)   = (gamma * Q(5,i,j-1:j+2,k) / (Q(1,i,j-1:j+2,k) * gamma_1)) &
                   + 0.5d0 * (Q(2,i,j-1:j+2,k)**2 + Q(3,i,j-1:j+2,k)**2 + Q(4,i,j-1:j+2,k)**2) + qc2(i,j-1:j+2,k)
          Hsgs   = -muty * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) * one_third) * dy(j) / Prt
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
          real(8), device :: Tx(3,2)
          Tx(:,:) = Q(5,i-1:i+1,j:j+1,k) / (R * Q(1,i-1:i+1,j:j+1,k))
          mx(:)   = mu23(Tx(:,:))
          mux    = 0.25d0 * (mx(1)    * (-Q(2,i-1,j,k) + Q(2,i,j,k) - Q(2,i-1,j+1,k) + Q(2,i,j+1,k)) &
                           + mx(2)    * (-Q(2,i,j,k) + Q(2,i+1,j,k) - Q(2,i,j+1,k) + Q(2,i+1,j+1,k))) * dx(i)
          muxsgs = 0.25d0 * (mxsgs(1) * (-Q(2,i-1,j,k) + Q(2,i,j,k) - Q(2,i-1,j+1,k) + Q(2,i,j+1,k)) &
                           + mxsgs(2) * (-Q(2,i,j,k) + Q(2,i+1,j,k) - Q(2,i,j+1,k) + Q(2,i+1,j+1,k))) * dx(i)
          mvx    = 0.25d0 * (mx(1)    * (-Q(3,i-1,j,k) + Q(3,i,j,k) - Q(3,i-1,j+1,k) + Q(3,i,j+1,k)) &
                           + mx(2)    * (-Q(3,i,j,k) + Q(3,i+1,j,k) - Q(3,i,j+1,k) + Q(3,i+1,j+1,k))) * dx(i)
          mvxsgs = 0.25d0 * (mxsgs(1) * (-Q(3,i-1,j,k) + Q(3,i,j,k) - Q(3,i-1,j+1,k) + Q(3,i,j+1,k)) &
                           + mxsgs(2) * (-Q(3,i,j,k) + Q(3,i+1,j,k) - Q(3,i,j+1,k) + Q(3,i+1,j+1,k))) * dx(i)
        end block
        block
          real(8), device :: Ty(2)
          Ty(:) = Q(5,i,j:j+1,k) / (R * Q(1,i,j:j+1,k))
          my    = mu2(Ty(:))
          kTy   = Cp_over_Pr * my * (-Ty(1) + Ty(2)) * dy(j)
        end block
        block
          real(8), device :: Tz(2,3)
          Tz(:,:) = Q(5,i,j:j+1,k-1:k+1) / (R * Q(1,i,j:j+1,k-1:k+1))
          mz(:)   = mu32(Tz(:,:))
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
        tyy    = 2.d0 * (2.d0 * mvy - mwz - mux) * one_third
        tyz    = mvz + mwy
        utyx   = 0.5d0 * (Q(2,i,j,k) + Q(2,i,j+1,k)) * tyx
        vtyy   = 0.5d0 * (Q(3,i,j,k) + Q(3,i,j+1,k)) * tyy
        wtyz   = 0.5d0 * (Q(4,i,j,k) + Q(4,i,j+1,k)) * tyz
        tyx    = tyx + muysgs + mvxsgs
        tyy    = tyy + 2.d0 * (2.d0 * mvysgs - mwzsgs - muxsgs) * one_third
        tyz    = tyz + mvzsgs + mwysgs
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
  end subroutine calc_Fv_LES4
 

  attributes(global) subroutine calc_Gv4(nx, ny, nz, dx, dy, dz, Q, G, seed)
    use calc_sutherland, only : mu6, mu2, mu32
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(inout), device :: G(5,nx-2,ny-2,nz-1)
    integer(8), intent(inout), device, optional :: seed(nx,ny,nz)
    integer i, j, k
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      block
        real(8), device :: mu(3)
        block ! dQdz
          real(8), device :: T6(6)
          T6(:) = Q(5,i,j,k-2:k+3) / (R * Q(1,i,j,k-2:k+3))
          mu(:) = mu6(T6(:))
          kTz   = heat_conduction6(mu(:), T6(:), dz(k))
        end block
        

        block ! tzx
          real(8), device :: uz3(3), wx3(3)
          uz3(:) = 0.125d0 * (9.d0 * (-Q(2,i,j,k-1:k+1) + Q(2,i,j,k:k+1)) &
                                   - (-Q(2,i,j,k-2:k)   + Q(2,i,j,k+1:k+3)) * one_third) * dz(k)
          block
            real(8) wx_1, wx_2, wx_3, wx_4, wx_5, wx_6
            wx_1 = dy5(Q(4,i-2,j,k-2), Q(4,i-1,j,k-2), Q(4,i+1,j,k-2), Q(4,i+2,j,k-2), dx(i))
            wx_2 = dy5(Q(4,i-2,j,k-1), Q(4,i-1,j,k-1), Q(4,i+1,j,k-1), Q(4,i+2,j,k-1), dx(i))
            wx_3 = dy5(Q(4,i-2,j,k  ), Q(4,i-1,j,k  ), Q(4,i+1,j,k  ), Q(4,i+2,j,k  ), dx(i))
            wx_4 = dy5(Q(4,i-2,j,k+1), Q(4,i-1,j,k+1), Q(4,i+1,j,k+1), Q(4,i+2,j,k+1), dx(i))
            wx_5 = dy5(Q(4,i-2,j,k+2), Q(4,i-1,j,k+2), Q(4,i+1,j,k+2), Q(4,i+2,j,k+2), dx(i))
            wx_6 = dy5(Q(4,i-2,j,k+3), Q(4,i-1,j,k+3), Q(4,i+1,j,k+3), Q(4,i+2,j,k+3), dx(i))
            wx3(:) = interpolation6_scalar(wx_1, wx_2, wx_3, wx_4, wx_5, wx_6)
          end block
          call tauxy_4(mu(:), wx3(:), uz3(:), &
                      Q(2,i,j,k-2), Q(2,i,j,k-1), Q(2,i,j,k), Q(2,i,j,k+1), Q(2,i,j,k+2), Q(2,i,j,k+3), tzx, utzx)
        end block


        block ! tzy
          real(8), device :: vz3(3), wy3(3)
          vz3(:) = 0.125d0 * (9.d0 * (-Q(3,i,j,k-1:k+1) + Q(3,i,j,k:k+1)) &
                                   - (-Q(3,i,j,k-2:k)   + Q(3,i,j,k+1:k+3)) * one_third) * dz(k)
          block
            real(8) wy_1, wy_2, wy_3, wy_4, wy_5, wy_6
            wy_1 = dy5(Q(4,i,j-2,k-2), Q(4,i,j-1,k-2), Q(4,i,j+1,k-2), Q(4,i,j+2,k-2), dy(j))
            wy_2 = dy5(Q(4,i,j-2,k-1), Q(4,i,j-1,k-1), Q(4,i,j+1,k-1), Q(4,i,j+2,k-1), dy(j))
            wy_3 = dy5(Q(4,i,j-2,k  ), Q(4,i,j-1,k  ), Q(4,i,j+1,k  ), Q(4,i,j+2,k  ), dy(j))
            wy_4 = dy5(Q(4,i,j-2,k+1), Q(4,i,j-1,k+1), Q(4,i,j+1,k+1), Q(4,i,j+2,k+1), dy(j))
            wy_5 = dy5(Q(4,i,j-2,k+2), Q(4,i,j-1,k+2), Q(4,i,j+1,k+2), Q(4,i,j+2,k+2), dy(j))
            wy_6 = dy5(Q(4,i,j-2,k+3), Q(4,i,j-1,k+3), Q(4,i,j+1,k+3), Q(4,i,j+2,k+3), dy(j))
            wy3(:) = interpolation6_scalar(wy_1, wy_2, wy_3, wy_4, wy_5, wy_6)
          end block
          call tauxy_4(mu(:), vz3(:), wy3(:), &
                      Q(3,i,j,k-2), Q(3,i,j,k-1), Q(3,i,j,k), Q(3,i,j,k+1), Q(3,i,j,k+2), Q(3,i,j,k+3), tzy, vtzy)
        end block


        block ! tzz
          real(8), device :: wz3(3), ux3(3), vy3(3)
          wz3(:) = 0.125d0 * (9.d0 * (-Q(4,i,j,k-1:k+1) + Q(4,i,j,k:k+1)) &
                                   - (-Q(4,i,j,k-2:k)   + Q(4,i,j,k+1:k+3)) * one_third) * dz(k)
          block
            real(8) ux_1, ux_2, ux_3, ux_4, ux_5, ux_6
            ux_1 = dy5(Q(2,i-2,j,k-2), Q(2,i-1,j,k-2), Q(2,i+1,j,k-2), Q(2,i+2,j,k-2), dx(i))
            ux_2 = dy5(Q(2,i-2,j,k-1), Q(2,i-1,j,k-1), Q(2,i+1,j,k-1), Q(2,i+2,j,k-1), dx(i))
            ux_3 = dy5(Q(2,i-2,j,k  ), Q(2,i-1,j,k  ), Q(2,i+1,j,k  ), Q(2,i+2,j,k  ), dx(i))
            ux_4 = dy5(Q(2,i-2,j,k+1), Q(2,i-1,j,k+1), Q(2,i+1,j,k+1), Q(2,i+2,j,k+1), dx(i))
            ux_5 = dy5(Q(2,i-2,j,k+2), Q(2,i-1,j,k+2), Q(2,i+1,j,k+2), Q(2,i+2,j,k+2), dx(i))
            ux_6 = dy5(Q(2,i-2,j,k+3), Q(2,i-1,j,k+3), Q(2,i+1,j,k+3), Q(2,i+2,j,k+3), dx(i))
            ux3(:) = interpolation6_scalar(ux_1, ux_2, ux_3, ux_4, ux_5, ux_6)
          end block
          block
            real(8) vy_1, vy_2, vy_3, vy_4, vy_5, vy_6
            vy_1 = dy5(Q(3,i,j-2,k-2), Q(3,i,j-1,k-2), Q(3,i,j+1,k-2), Q(3,i,j+2,k-2), dy(j))
            vy_2 = dy5(Q(3,i,j-2,k-1), Q(3,i,j-1,k-1), Q(3,i,j+1,k-1), Q(3,i,j+2,k-1), dy(j))
            vy_3 = dy5(Q(3,i,j-2,k  ), Q(3,i,j-1,k  ), Q(3,i,j+1,k  ), Q(3,i,j+2,k  ), dy(j))
            vy_4 = dy5(Q(3,i,j-2,k+1), Q(3,i,j-1,k+1), Q(3,i,j+1,k+1), Q(3,i,j+2,k+1), dy(j))
            vy_5 = dy5(Q(3,i,j-2,k+2), Q(3,i,j-1,k+2), Q(3,i,j+1,k+2), Q(3,i,j+2,k+2), dy(j))
            vy_6 = dy5(Q(3,i,j-2,k+3), Q(3,i,j-1,k+3), Q(3,i,j+1,k+3), Q(3,i,j+2,k+3), dy(j))
            vy3(:) = interpolation6_scalar(vy_1, vy_2, vy_3, vy_4, vy_5, vy_6)
          end block
          call tauxx_4(mu(:), wz3(:), ux3(:), vy3(:), &
                      Q(4,i,j,k-2), Q(4,i,j,k-1), Q(4,i,j,k), Q(4,i,j,k+1), Q(4,i,j,k+2), Q(4,i,j,k+3), tzz, wtzz)
        end block
      end block
    else
      block
        real(8) mz, muz, mvz, mwz, mwx, mux, mvy, mwy
        block
          real(8), device :: Tx(3,2), mx(2)
          Tx(:,:) = Q(5,i-1:i+1,j,k:k+1) / (R * Q(1,i-1:i+1,j,k:k+1))
          mx(:)   = mu32(Tx(:,:))
          mux = 0.25d0 * (mx(1) * (-Q(2,i-1,j,k) + Q(2,i,j,k) - Q(2,i-1,j,k+1) + Q(2,i,j,k+1)) &
                        + mx(2) * (-Q(2,i,j,k) + Q(2,i+1,j,k) - Q(2,i,j,k+1) + Q(2,i+1,j,k+1))) * dx(i)
          mwx = 0.25d0 * (mx(1) * (-Q(4,i-1,j,k) + Q(4,i,j,k) - Q(4,i-1,j,k+1) + Q(4,i,j,k+1)) &
                        + mx(2) * (-Q(4,i,j,k) + Q(4,i+1,j,k) - Q(4,i,j,k+1) + Q(4,i+1,j,k+1))) * dx(i)
        end block
        block
          real(8), device :: Ty(3,2), my(2)
          Ty(:,:) = Q(5,i,j-1:j+1,k:k+1) / (R * Q(1,i,j-1:j+1,k:k+1))
          my(:)   = mu32(Ty(:,:))
          mvy = 0.25d0 * (my(1) * (-Q(3,i,j-1,k) + Q(3,i,j,k) - Q(3,i,j-1,k+1) + Q(3,i,j,k+1)) &
                        + my(2) * (-Q(3,i,j,k) + Q(3,i,j+1,k) - Q(3,i,j,k+1) + Q(3,i,j+1,k+1))) * dy(j)
          mwy = 0.25d0 * (my(1) * (-Q(4,i,j-1,k) + Q(4,i,j,k) - Q(4,i,j-1,k+1) + Q(4,i,j,k+1)) &
                        + my(2) * (-Q(4,i,j,k) + Q(4,i,j+1,k) - Q(4,i,j,k+1) + Q(4,i,j+1,k+1))) * dy(j)
        end block
        block
          real(8), device :: Tz(2)
          Tz(:)   = Q(5,i,j,k:k+1) / (R * Q(1,i,j,k:k+1))
          mz      = mu2(Tz(:))
          kTz     = Cp_over_Pr * mz * (-Tz(1) + Tz(2)) * dz(k)
        end block
        muz = mz * (-Q(2,i,j,k) + Q(2,i,j,k+1)) * dz(k)
        mvz = mz * (-Q(3,i,j,k) + Q(3,i,j,k+1)) * dz(k)
        mwz = mz * (-Q(4,i,j,k) + Q(4,i,j,k+1)) * dz(k)
        tzx  = mwx + muz
        tzy  = mvz + mwy
        tzz  = 2.d0 * (2.d0 * mwz - mux - mvy) * one_third
        utzx = 0.5d0 * (Q(2,i,j,k) + Q(2,i,j,k+1)) * tzx
        vtzy = 0.5d0 * (Q(3,i,j,k) + Q(3,i,j,k+1)) * tzy
        wtzz = 0.5d0 * (Q(4,i,j,k) + Q(4,i,j,k+1)) * tzz
      end block
    endif
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz)
  end subroutine calc_Gv4


  attributes(global) subroutine calc_Gv_LES4(nx, ny, nz, dx, dy, dz, Q, mut, qc2, G)
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
    if (3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      block
        real(8), device :: mu(3)
        real(8) mutz
        mutz = 0.0625d0 * (-mut(i,j,k-1) + 9.d0 * (mut(i,j,k) + mut(i,j,k+1)) -mut(i,j,k+2))
        block ! dQdz
          real(8), device :: T6(6)
          T6(:) = Q(5,i,j,k-2:k+3) / (R * Q(1,i,j,k-2:k+3))
          mu(:) = mu6(T6(:))
          kTz   = heat_conduction6(mu(:), T6(:), dz(k))
        end block
        

        block ! tzx
          real(8), device :: uz3(3), wx3(3)
          uz3(:) = 0.125d0 * (9.d0 * (-Q(2,i,j,k-1:k+1) + Q(2,i,j,k:k+1)) &
                                   - (-Q(2,i,j,k-2:k)   + Q(2,i,j,k+1:k+3)) * one_third) * dz(k)
          block
            real(8) wx_1, wx_2, wx_3, wx_4, wx_5, wx_6
            wx_1 = dy5(Q(4,i-2,j,k-2), Q(4,i-1,j,k-2), Q(4,i+1,j,k-2), Q(4,i+2,j,k-2), dx(i))
            wx_2 = dy5(Q(4,i-2,j,k-1), Q(4,i-1,j,k-1), Q(4,i+1,j,k-1), Q(4,i+2,j,k-1), dx(i))
            wx_3 = dy5(Q(4,i-2,j,k  ), Q(4,i-1,j,k  ), Q(4,i+1,j,k  ), Q(4,i+2,j,k  ), dx(i))
            wx_4 = dy5(Q(4,i-2,j,k+1), Q(4,i-1,j,k+1), Q(4,i+1,j,k+1), Q(4,i+2,j,k+1), dx(i))
            wx_5 = dy5(Q(4,i-2,j,k+2), Q(4,i-1,j,k+2), Q(4,i+1,j,k+2), Q(4,i+2,j,k+2), dx(i))
            wx_6 = dy5(Q(4,i-2,j,k+3), Q(4,i-1,j,k+3), Q(4,i+1,j,k+3), Q(4,i+2,j,k+3), dx(i))
            wx3(:) = interpolation6_scalar(wx_1, wx_2, wx_3, wx_4, wx_5, wx_6)
          end block
          call tauxy_4(mu(:), wx3(:), uz3(:), &
                      Q(2,i,j,k-2), Q(2,i,j,k-1), Q(2,i,j,k), Q(2,i,j,k+1), Q(2,i,j,k+2), Q(2,i,j,k+3), tzx, utzx)
          tzx  = tzx + mutz * (wx3(2) + uz3(2))
        end block


        block ! tzy
          real(8), device :: vz3(3), wy3(3)
          vz3(:) = 0.125d0 * (9.d0 * (-Q(3,i,j,k-1:k+1) + Q(3,i,j,k:k+1)) &
                                   - (-Q(3,i,j,k-2:k)   + Q(3,i,j,k+1:k+3)) * one_third) * dz(k)
          block
            real(8) wy_1, wy_2, wy_3, wy_4, wy_5, wy_6
            wy_1 = dy5(Q(4,i,j-2,k-2), Q(4,i,j-1,k-2), Q(4,i,j+1,k-2), Q(4,i,j+2,k-2), dy(j))
            wy_2 = dy5(Q(4,i,j-2,k-1), Q(4,i,j-1,k-1), Q(4,i,j+1,k-1), Q(4,i,j+2,k-1), dy(j))
            wy_3 = dy5(Q(4,i,j-2,k  ), Q(4,i,j-1,k  ), Q(4,i,j+1,k  ), Q(4,i,j+2,k  ), dy(j))
            wy_4 = dy5(Q(4,i,j-2,k+1), Q(4,i,j-1,k+1), Q(4,i,j+1,k+1), Q(4,i,j+2,k+1), dy(j))
            wy_5 = dy5(Q(4,i,j-2,k+2), Q(4,i,j-1,k+2), Q(4,i,j+1,k+2), Q(4,i,j+2,k+2), dy(j))
            wy_6 = dy5(Q(4,i,j-2,k+3), Q(4,i,j-1,k+3), Q(4,i,j+1,k+3), Q(4,i,j+2,k+3), dy(j))
            wy3(:) = interpolation6_scalar(wy_1, wy_2, wy_3, wy_4, wy_5, wy_6)
          end block
          call tauxy_4(mu(:), vz3(:), wy3(:), &
                      Q(3,i,j,k-2), Q(3,i,j,k-1), Q(3,i,j,k), Q(3,i,j,k+1), Q(3,i,j,k+2), Q(3,i,j,k+3), tzy, vtzy)
          tzy  = tzy + mutz * (vz3(2) + wy3(2))
        end block


        block ! tzz
          real(8), device :: wz3(3), ux3(3), vy3(3)
          wz3(:) = 0.125d0 * (9.d0 * (-Q(4,i,j,k-1:k+1) + Q(4,i,j,k:k+1)) &
                                   - (-Q(4,i,j,k-2:k)   + Q(4,i,j,k+1:k+3)) * one_third) * dz(k)
          block
            real(8) ux_1, ux_2, ux_3, ux_4, ux_5, ux_6
            ux_1 = dy5(Q(2,i-2,j,k-2), Q(2,i-1,j,k-2), Q(2,i+1,j,k-2), Q(2,i+2,j,k-2), dx(i))
            ux_2 = dy5(Q(2,i-2,j,k-1), Q(2,i-1,j,k-1), Q(2,i+1,j,k-1), Q(2,i+2,j,k-1), dx(i))
            ux_3 = dy5(Q(2,i-2,j,k  ), Q(2,i-1,j,k  ), Q(2,i+1,j,k  ), Q(2,i+2,j,k  ), dx(i))
            ux_4 = dy5(Q(2,i-2,j,k+1), Q(2,i-1,j,k+1), Q(2,i+1,j,k+1), Q(2,i+2,j,k+1), dx(i))
            ux_5 = dy5(Q(2,i-2,j,k+2), Q(2,i-1,j,k+2), Q(2,i+1,j,k+2), Q(2,i+2,j,k+2), dx(i))
            ux_6 = dy5(Q(2,i-2,j,k+3), Q(2,i-1,j,k+3), Q(2,i+1,j,k+3), Q(2,i+2,j,k+3), dx(i))
            ux3(:) = interpolation6_scalar(ux_1, ux_2, ux_3, ux_4, ux_5, ux_6)
          end block
          block
            real(8) vy_1, vy_2, vy_3, vy_4, vy_5, vy_6
            vy_1 = dy5(Q(3,i,j-2,k-2), Q(3,i,j-1,k-2), Q(3,i,j+1,k-2), Q(3,i,j+2,k-2), dy(j))
            vy_2 = dy5(Q(3,i,j-2,k-1), Q(3,i,j-1,k-1), Q(3,i,j+1,k-1), Q(3,i,j+2,k-1), dy(j))
            vy_3 = dy5(Q(3,i,j-2,k  ), Q(3,i,j-1,k  ), Q(3,i,j+1,k  ), Q(3,i,j+2,k  ), dy(j))
            vy_4 = dy5(Q(3,i,j-2,k+1), Q(3,i,j-1,k+1), Q(3,i,j+1,k+1), Q(3,i,j+2,k+1), dy(j))
            vy_5 = dy5(Q(3,i,j-2,k+2), Q(3,i,j-1,k+2), Q(3,i,j+1,k+2), Q(3,i,j+2,k+2), dy(j))
            vy_6 = dy5(Q(3,i,j-2,k+3), Q(3,i,j-1,k+3), Q(3,i,j+1,k+3), Q(3,i,j+2,k+3), dy(j))
            vy3(:) = interpolation6_scalar(vy_1, vy_2, vy_3, vy_4, vy_5, vy_6)
          end block
          call tauxx_4(mu(:), wz3(:), ux3(:), vy3(:), &
                      Q(4,i,j,k-2), Q(4,i,j,k-1), Q(4,i,j,k), Q(4,i,j,k+1), Q(4,i,j,k+2), Q(4,i,j,k+3), tzz, wtzz)
          tzz  = tzz + 2.d0 * mutz * (2.d0 * wz3(2) - ux3(2) - vy3(2)) * one_third
        end block
 

        block
          real(8), device :: H(4)
          H(:) = (gamma * Q(5,i,j,k-1:k+2) / (Q(1,i,j,k-1:k+2) * gamma_1)) &
                 + 0.5d0 * (Q(2,i,j,k-1:k+2)**2 + Q(3,i,j,k-1:k+2)**2 + Q(4,i,j,k-1:k+2)**2) + qc2(i,j,k-1:k+2)
          Hsgs = -mutz * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) * one_third) * dz(k) / Prt
        end block
      end block
    else
      block
        real(8), dimension(2), device :: mx, mxsgs, my, mysgs
        real(8) mz, mzsgs, muz, muzsgs, mvz, mvzsgs, mwz, mwzsgs, mwx, mwxsgs, mux, muxsgs, mvy, mvysgs, mwy, mwysgs
        ! SGS
        mz    = 0.5d0 * (mut(i,j,k) + mut(i,j,k+1))
        mx(:) = (/0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j,k+1) + mut(i,j,k+1)), &
                  0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j,k+1) + mut(i+1,j,k+1))/)
        my(:) = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i,j-1,k+1) + mut(i,j,k+1)), &
                  0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i,j,k+1) + mut(i,j+1,k+1))/)
        block
          real(8), device :: Tx(3,2)
          Tx(:,:) = Q(5,i-1:i+1,j,k:k+1) / (R * Q(1,i-1:i+1,j,k:k+1))
          mx(:)   = mu32(Tx(:,:))
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
          real(8), device :: Ty(3,2)
          Ty(:,:) = Q(5,i,j-1:j+1,k:k+1) / (R * Q(1,i,j-1:j+1,k:k+1))
          my(:)   = mu32(Ty(:,:))
          mvy    = 0.25d0 * (my(1)    * (-Q(3,i,j-1,k) + Q(3,i,j,k) - Q(3,i,j-1,k+1) + Q(3,i,j,k+1)) &
                           + my(2)    * (-Q(3,i,j,k) + Q(3,i,j+1,k) - Q(3,i,j,k+1) + Q(3,i,j+1,k+1))) * dy(j)
          mvysgs = 0.25d0 * (mysgs(1) * (-Q(3,i,j-1,k) + Q(3,i,j,k) - Q(3,i,j-1,k+1) + Q(3,i,j,k+1)) &
                           + mysgs(2) * (-Q(3,i,j,k) + Q(3,i,j+1,k) - Q(3,i,j,k+1) + Q(3,i,j+1,k+1))) * dy(j)
          mwy    = 0.25d0 * (my(1)    * (-Q(4,i,j-1,k) + Q(4,i,j,k) - Q(4,i,j-1,k+1) + Q(4,i,j,k+1)) &
                           + my(2)    * (-Q(4,i,j,k) + Q(4,i,j+1,k) - Q(4,i,j,k+1) + Q(4,i,j+1,k+1))) * dy(j)
          mwysgs = 0.25d0 * (mysgs(1) * (-Q(4,i,j-1,k) + Q(4,i,j,k) - Q(4,i,j-1,k+1) + Q(4,i,j,k+1)) &
                           + mysgs(2) * (-Q(4,i,j,k) + Q(4,i,j+1,k) - Q(4,i,j,k+1) + Q(4,i,j+1,k+1))) * dy(j)
        end block
        block
          real(8), device :: Tz(2)
          Tz(:)   = Q(5,i,j,k:k+1) / (R * Q(1,i,j,k:k+1))
          mz      = mu2(Tz(:))
          kTz     = Cp_over_Pr * mz * (-Tz(1) + Tz(2)) * dz(k)
        end block
        muz    = mz    * (-Q(2,i,j,k) + Q(2,i,j,k+1)) * dz(k)
        muzsgs = mzsgs * (-Q(2,i,j,k) + Q(2,i,j,k+1)) * dz(k)
        mvz    = mz    * (-Q(3,i,j,k) + Q(3,i,j,k+1)) * dz(k)
        mvzsgs = mzsgs * (-Q(3,i,j,k) + Q(3,i,j,k+1)) * dz(k)
        mwz    = mz    * (-Q(4,i,j,k) + Q(4,i,j,k+1)) * dz(k)
        mwzsgs = mzsgs * (-Q(4,i,j,k) + Q(4,i,j,k+1)) * dz(k)
        tzx    = mwx + muz
        tzy    = mvz + mwy
        tzz    = 2.d0 * (2.d0 * mwz - mux - mvy) * one_third
        utzx = 0.5d0 * (Q(2,i,j,k) + Q(2,i,j,k+1)) * tzx
        vtzy = 0.5d0 * (Q(3,i,j,k) + Q(3,i,j,k+1)) * tzy
        wtzz = 0.5d0 * (Q(4,i,j,k) + Q(4,i,j,k+1)) * tzz
        tzx    = tzx + mwx + muz
        tzy    = tzy + mvz + mwy
        tzz    = tzz + 2.d0 * (2.d0 * mwz - mux - mvy) * one_third
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
  end subroutine calc_Gv_LES4
end module calc_visc4

