module calc_visc
  use mod_globals, only : accuracy, offset, id_visc, id_turbulence, id_av, gamma, R, Pr, Prt
  use mod_constant, only : Cp, gamma_1
  implicit none
contains
  !dir$ inline
  attributes(device) function interpolation6(a) result(ans)
    real(8), intent(in), device :: a(6)
    real(8) ans(3)
    ans(:) = 0.0625d0 * (-a(1:3) + 9.d0 * (a(2:4) + a(3:5)) -a(4:6))
  end function interpolation6

  !dir$ inline
  attributes(device) function dx6(a, dx) result(ans)
    real(8), intent(in), device :: a(6)
    real(8), intent(in), value  :: dx
    real(8) ans(3)
    ans(:) = 0.125d0 * (9.d0 * (-a(2:4) + a(3:5)) - (-a(1:3) + a(4:6)) / 3.d0) * dx
  end function dx6

  !dir$ inline
  attributes(device) function dy5(a, dy) result(ans)
    real(8), intent(in), device :: a(5)
    real(8), intent(in), value  :: dy
    real(8) ans
    ans = (2.d0 * (-a(2) + a(4)) - 0.25d0 * (-a(1) + a(5))) * dy / 3.d0
  end function dy5

  !dir$ inline
  attributes(device) function dy23(mu, a, dy) result(ans)
    real(8), intent(in), device :: mu(2), a(2,3)
    real(8), intent(in), value  :: dy
    real(8) ans
    ans = 0.25d0 * (mu(1) * (-a(1,1) + a(1,2) -a(2,1) + a(2,2)) &
                  + mu(2) * (-a(1,2) + a(1,3) -a(2,2) + a(2,3)))
  end function dy23

  !dir$ inline
  attributes(device) function dy65(a, dy) result(ans)
    real(8), intent(in), device :: a(6,5)
    real(8), intent(in), value  :: dy
    real(8) ans(3), a1(5), a2(5), a3(5), a4(5), a5(5), a6(5), ay(6)
    a1 = a(1,:); a2 = a(2,:); a3 = a(3,:); a4 = a(4,:); a5 = a(5,:); a6 = a(6,:)
    ay(1) = dy5(a1(:), dy); ay(2) = dy5(a2(:), dy); ay(3) = dy5(a3(:), dy)
    ay(4) = dy5(a4(:), dy); ay(5) = dy5(a5(:), dy); ay(6) = dy5(a6(:), dy)
    ans = interpolation6(ay(:))
  end function dy65
  
  !dir$ inline
  attributes(device) function dy32(mu, a, dy) result(ans)
    real(8), intent(in), device :: mu(2), a(3,2)
    real(8), intent(in), value  :: dy
    real(8) ans
    ans = 0.25d0 * (mu(1) * (-a(1,1) + a(2,1) - a(1,2) + a(2,2)) &
                  + mu(2) * (-a(2,1) + a(3,1) - a(2,2) + a(3,2)))
  end function dy32

  !dir$ inline
  attributes(device) function dy56(a, dy) result(ans)
    real(8), intent(in), device :: a(5,6)
    real(8), intent(in), value  :: dy
    real(8) ans(3), a1(5), a2(5), a3(5), a4(5), a5(5), a6(5), ay(6)
    a1 = a(:,1); a2 = a(:,2); a3 = a(:,3); a4 = a(:,4); a5 = a(:,5); a6 = a(:,6)
    ay(1) = dy5(a1(:), dy); ay(2) = dy5(a2(:), dy); ay(3) = dy5(a3(:), dy)
    ay(4) = dy5(a4(:), dy); ay(5) = dy5(a5(:), dy); ay(6) = dy5(a6(:), dy)
    ans = interpolation6(ay(:))
  end function dy56

  !dir$ inline
  attributes(device) function flux4(a) result(ans)
    real(8), intent(in), device :: a(3)
    real(8) ans
    ans = 0.125d0 * ((9.d0 - 1.d0 / 3.d0) * a(2) - (a(1) + a(3)) / 3.d0)
  end function flux4

  !dir$ inline
  attributes(device) subroutine tauxx4(mu, ux, vy, wz, u6, txx, utxx)
    real(8), intent(in), dimension(3), device :: mu, ux, vy, wz
    real(8), intent(in), dimension(6), device :: u6
    real(8), intent(out) :: txx, utxx
    real(8), dimension(3) :: tau, utau
    tau(:)  = 2.d0 * mu(:) * (2.d0 * ux(:) - vy(:) - wz(:)) / 3.d0
    utau(:) = interpolation6(u6(:)) * tau(:)
    txx     = flux4(tau(:))
    utxx    = flux4(utau(:))
  end subroutine tauxx4

  !dir$ inline
  attributes(device) subroutine tauxy4(mu, uy, vx, v6, txy, vtxy)
    real(8), intent(in), dimension(3), device :: mu, uy, vx
    real(8), intent(in), dimension(6), device :: v6
    real(8), intent(out) :: txy, vtxy
    real(8), dimension(3) :: tau, vtau
    tau(:)  = mu(:) * (uy(:) + vx(:))
    vtau(:) = interpolation6(v6(:)) * tau(:)
    txy     = flux4(tau(:))
    vtxy    = flux4(vtau(:))
  end subroutine tauxy4

  !dir$ inline
  attributes(device) function heat_conduction6(mu, T, dx) result(ans)
    real(8), intent(in), device :: mu(3), T(6)
    real(8), intent(in), value  :: dx
    real(8), dimension(3) ::  kTx
    real(8) :: ans
    kTx(:) = Cp * mu(:) * dx6(T(:), dx) / Pr
    ans = flux4(kTx(:))
  end function heat_conduction6


  !artificial viscosity!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  attributes(device) function av_vonNeumann(sensor, rho, u, p) result(q)
    use mod_globals, only : threshold
    real(8), intent(in), value                :: sensor
    real(8), intent(in), dimension(2), device :: rho, u, p
    real(8) du, c, q
    du = (-u(1) + u(2))
    if (sensor > threshold .and. du < 0.d0) then
      c = 0.5d0 * (sqrt(gamma * p(1) / rho(1)) + sqrt(gamma * p(2) / rho(2)))
      q = 0.6d0 * 0.5d0 * (-(rho(1) + rho(2)) * c * du + (gamma + 1.d0) * du**2)
    else
      q = 0.d0
    endif
  end function av_vonNeumann

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(global) subroutine calc_Ev(nx, ny, nz, dx, dy, dz, Q, E)
    use calc_sutherland, only : mu6, mu2, mu23
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(inout), device :: E(5,nx-accuracy+1,ny-accuracy,nz-accuracy)
    integer i, j, k
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx
    ! 4th-order accuracy
    real(8), dimension(6,5), device :: u651, v651, u615, w615
    real(8), dimension(6), device   :: T6, u6, v6, w6
    real(8), dimension(3), device   :: ux3, vx3, wx3, uy3, vy3, uz3, wz3, mu
    ! 2nd-order accuracy
    real(8), dimension(2,3,3) , device :: T233
    real(8), dimension(2,3), device    :: u231, v231, u213, w213, Ty, Tz
    real(8), dimension(2), device      :: Tx, u2, v2, w2, my, mz
    real(8) mx, mux, mvx, mwx, muy, mvy, muz, mwz
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset - 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    if (id_visc ==2 .and. 3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      u651(:,:) = Q(2,i-2:i+3,j-2:j+2,k)
      u615(:,:) = Q(2,i-2:i+3,j,k-2:k+2)
      v651(:,:) = Q(3,i-2:i+3,j-2:j+2,k)
      w615(:,:) = Q(4,i-2:i+3,j,k-2:k+2)
      T6(:)     = Q(5,i-2:i+3,j,k) / (R * Q(1,i-2:i+3,j,k))
      u6(:)     = u651(:,3) 
      v6(:)     = v651(:,3) 
      w6(:)     = w615(:,3) 
      mu(:)     = mu6(T6(:))
      ux3(:)    = dx6(u6(:), dx(i))
      vx3(:)    = dx6(v6(:), dx(i))
      wx3(:)    = dx6(w6(:), dx(i))
      uy3(:)    = dy65(u651(:,:), dy(j))
      vy3(:)    = dy65(v651(:,:), dy(j))
      uz3(:)    = dy65(u615(:,:), dz(k))
      wz3(:)    = dy65(w615(:,:), dz(k))
      call tauxx4(mu(:), ux3(:), vy3(:), wz3(:), u6(:), txx, utxx)
      call tauxy4(mu(:), uy3(:), vx3(:), v6(:), txy, vtxy)
      call tauxy4(mu(:), wx3(:), uz3(:), w6(:), txz, wtxz)
      kTx = heat_conduction6(mu(:), T6(:), dx(i))
    else
      T233(:,:,:) = Q(5,i:i+1,j-1:j+1,k-1:k+1) / (R * Q(1,i:i+1,j-1:j+1,k-1:k+1))
      u231(:,:)   = Q(2,i:i+1,j-1:j+1,k)
      u213(:,:)   = Q(2,i:i+1,j,k-1:k+1)
      v231(:,:)   = Q(3,i:i+1,j-1:j+1,k)
      w213(:,:)   = Q(4,i:i+1,j,k-1:k+1)
      Tx(:)       = T233(:,2,2)
      Ty(:,:)     = T233(:,:,2)
      Tz(:,:)     = T233(:,2,:)
      u2(:)       = u231(:,2)
      v2(:)       = v231(:,2)
      w2(:)       = w213(:,2)
      mx          = mu2(Tx(:))
      my(:)       = mu23(Ty(:,:))
      mz(:)       = mu23(Tz(:,:))
      mux         = mx * (-u2(1) + u2(2)) * dx(i)
      mvx         = mx * (-v2(1) + v2(2)) * dx(i)
      mwx         = mx * (-w2(1) + w2(2)) * dx(i)
      muy         = dy23(my(:), u231(:,:), dy(j))
      mvy         = dy23(my(:), v231(:,:), dy(j))
      muz         = dy23(mz(:), u213(:,:), dz(k))
      mwz         = dy23(mz(:), w213(:,:), dz(k))
      txx         = 2.d0 * (2.d0 * mux - mvy - mwz) / 3.d0
      txy         = muy + mvx
      txz         = mwx + muz
      utxx        = 0.5d0 * (u2(1) + u2(2)) * txx
      vtxy        = 0.5d0 * (v2(1) + v2(2)) * txy
      wtxz        = 0.5d0 * (w2(1) + w2(2)) * txz
      kTx         = Cp * mx * (-Tx(1) + Tx(2)) * dx(i) / Pr
    endif
    E(2,i-offset+1,j-offset,k-offset) = E(2,i-offset+1,j-offset,k-offset) - txx
    E(3,i-offset+1,j-offset,k-offset) = E(3,i-offset+1,j-offset,k-offset) - txy
    E(4,i-offset+1,j-offset,k-offset) = E(4,i-offset+1,j-offset,k-offset) - txz
    E(5,i-offset+1,j-offset,k-offset) = E(5,i-offset+1,j-offset,k-offset) - (utxx + vtxy + wtxz + kTx)
  end subroutine calc_Ev
  
  attributes(global) subroutine calc_Ev_LES(nx, ny, nz, dx, dy, dz, Q, mut, qc2, E)
    use calc_sutherland, only : mu6, mu2, mu23
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(in), device    :: mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device :: E(5,nx-accuracy+1,ny-accuracy,nz-accuracy)
    integer i, j, k
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx, mutx, H(4), txxsgs = 0.d0, txysgs = 0.d0, txzsgs = 0.d0, Hsgs = 0.d0
    ! 4th-order accuracy
    real(8), dimension(6,5), device :: u651, v651, u615, w615
    real(8), dimension(6), device   :: T6, u6, v6, w6
    real(8), dimension(3), device   :: ux3, vx3, wx3, uy3, vy3, uz3, wz3, mu
    ! 2nd-order accuracy
    real(8), dimension(2,3,3) , device :: T233
    real(8), dimension(2,3), device    :: u231, v231, u213, w213, Ty, Tz
    real(8), dimension(2), device      :: Tx, u2, v2, w2, my, mz
    real(8) mx, mux, mvx, mwx, muy, mvy, muz, mwz
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset - 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    if (id_visc ==2 .and. 3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      u651(:,:) = Q(2,i-2:i+3,j-2:j+2,k)
      u615(:,:) = Q(2,i-2:i+3,j,k-2:k+2)
      v651(:,:) = Q(3,i-2:i+3,j-2:j+2,k)
      w615(:,:) = Q(4,i-2:i+3,j,k-2:k+2)
      T6(:)     = Q(5,i-2:i+3,j,k) / (R * Q(1,i-2:i+3,j,k))
      u6(:)     = u651(:,3) 
      v6(:)     = v651(:,3) 
      w6(:)     = w615(:,3) 
      mu(:)     = mu6(T6(:))
      ux3(:)    = dx6(u6(:), dx(i))
      vx3(:)    = dx6(v6(:), dx(i))
      wx3(:)    = dx6(w6(:), dx(i))
      uy3(:)    = dy65(u651(:,:), dy(j))
      vy3(:)    = dy65(v651(:,:), dy(j))
      uz3(:)    = dy65(u615(:,:), dz(k))
      wz3(:)    = dy65(w615(:,:), dz(k))
      call tauxx4(mu(:), ux3(:), vy3(:), wz3(:), u6(:), txx, utxx)
      call tauxy4(mu(:), uy3(:), vx3(:), v6(:), txy, vtxy)
      call tauxy4(mu(:), wx3(:), uz3(:), w6(:), txz, wtxz)
      kTx = heat_conduction6(mu(:), T6(:), dx(i))
      mutx   = 0.0625d0 * (-mut(i-1,j,k) + 9.d0 * (mut(i,j,k) + mut(i+1,j,k)) -mut(i+2,j,k))
      txxsgs = 2.d0 * mutx * (2.d0 * ux3(2) - vy3(2) - wz3(2)) / 3.d0
      txysgs = mutx * (uy3(2) + vx3(2))
      txzsgs = mutx * (wx3(2) + uz3(2))
      H(:)   = (gamma * Q(5,i-1:i+2,j,k) / (Q(1,i-1:i+2,j,k) * gamma_1)) &
               + 0.5d0 * (Q(2,i-1:i+2,j,k)**2 + Q(3,i-1:i+2,j,k)**2 + Q(4,i-1:i+2,j,k)**2) + qc2(i-1:i+2,j,k)
      Hsgs   = -mutx * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) / 3.d0) * dx(i) / Prt
    else
      T233(:,:,:) = Q(5,i:i+1,j-1:j+1,k-1:k+1) / (R * Q(1,i:i+1,j-1:j+1,k-1:k+1))
      u231(:,:)   = Q(2,i:i+1,j-1:j+1,k)
      u213(:,:)   = Q(2,i:i+1,j,k-1:k+1)
      v231(:,:)   = Q(3,i:i+1,j-1:j+1,k)
      w213(:,:)   = Q(4,i:i+1,j,k-1:k+1)
      Tx(:)       = T233(:,2,2)
      Ty(:,:)     = T233(:,:,2)
      Tz(:,:)     = T233(:,2,:)
      u2(:)       = u231(:,2)
      v2(:)       = v231(:,2)
      w2(:)       = w213(:,2)
      mx          = mu2(Tx(:))
      my(:)       = mu23(Ty(:,:))
      mz(:)       = mu23(Tz(:,:))
      mux         = mx * (-u2(1) + u2(2)) * dx(i)
      mvx         = mx * (-v2(1) + v2(2)) * dx(i)
      mwx         = mx * (-w2(1) + w2(2)) * dx(i)
      muy         = dy23(my(:), u231(:,:), dy(j))
      mvy         = dy23(my(:), v231(:,:), dy(j))
      muz         = dy23(mz(:), u213(:,:), dz(k))
      mwz         = dy23(mz(:), w213(:,:), dz(k))
      txx         = 2.d0 * (2.d0 * mux - mvy - mwz) / 3.d0
      txy         = muy + mvx
      txz         = mwx + muz
      utxx        = 0.5d0 * (u2(1) + u2(2)) * txx
      vtxy        = 0.5d0 * (v2(1) + v2(2)) * txy
      wtxz        = 0.5d0 * (w2(1) + w2(2)) * txz
      kTx         = Cp * mx * (-Tx(1) + Tx(2)) * dx(i) / Pr
      mx     = 0.5d0 * (mut(i,j,k) + mut(i+1,j,k))
      my(:)  = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i+1,j-1,k) + mut(i+1,j,k)), &
                 0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i+1,j,k) + mut(i+1,j+1,k))/)
      mz(:)  = (/0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i+1,j,k-1) + mut(i+1,j,k)), &
                 0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i+1,j,k) + mut(i+1,j,k+1))/)
      mux    = mx * (-u2(1) + u2(2)) * dx(i)
      mvx    = mx * (-v2(1) + v2(2)) * dx(i)
      mwx    = mx * (-w2(1) + w2(2)) * dx(i)
      muy    = dy23(my(:), u231(:,:), dy(j))
      mvy    = dy23(my(:), v231(:,:), dy(j))
      muz    = dy23(mz(:), u213(:,:), dz(k))
      mwz    = dy23(mz(:), w213(:,:), dz(k))
      txxsgs = 2.d0 * (2.d0 * mux - mvy - mwz) / 3.d0
      txysgs = muy + mvx
      txzsgs = mwx + muz
      H(2:3) = (gamma * Q(5,i:i+1,j,k) / (Q(1,i:i+1,j,k) * gamma_1)) &
               + 0.5d0 * (u2(:)**2 + v2(:)**2 + w2(:)**2) + qc2(i:i+1,j,k)
      Hsgs   = -mx * (-H(2) + H(3)) * dx(i) / Prt
    endif
    E(2,i-offset+1,j-offset,k-offset) = E(2,i-offset+1,j-offset,k-offset) - (txx+txxsgs)
    E(3,i-offset+1,j-offset,k-offset) = E(3,i-offset+1,j-offset,k-offset) - (txy+txysgs)
    E(4,i-offset+1,j-offset,k-offset) = E(4,i-offset+1,j-offset,k-offset) - (txz+txzsgs)
    E(5,i-offset+1,j-offset,k-offset) = E(5,i-offset+1,j-offset,k-offset) - (utxx + vtxy + wtxz + kTx + Hsgs)
  end subroutine calc_Ev_LES
  
  attributes(global) subroutine calc_Fv(nx, ny, nz, dy, dx, dz, Q, F)
    use calc_sutherland, only : mu6, mu2, mu23, mu32
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(inout), device :: F(5,nx-accuracy,ny-accuracy+1,nz-accuracy)
    integer i, j, k
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, muty
    ! 4th-order accuracy
    real(8), dimension(5,6), device :: u561, v561
    real(8), dimension(6,5), device :: v165, w165
    real(8), dimension(6), device   :: T6, u6, v6, w6
    real(8), dimension(3), device   :: uy3, vy3, wy3, vz3, wz3, ux3, vx3, mu
    ! 2nd-order accuracy
    real(8), dimension(3,2,3), device :: T323
    real(8), dimension(3,2), device   :: u321, v321, Tx
    real(8), dimension(2,3), device   :: v123, w123, Tz
    real(8), dimension(2), device     :: Ty, u2, v2, w2, mz, mx
    real(8) my, muy, mvy, mwy, mvz, mwz, mux, mvx
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset - 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      u561(:,:) = Q(2,i-2:i+2,j-2:j+3,k)
      v561(:,:) = Q(3,i-2:i+2,j-2:j+3,k)
      v165(:,:) = Q(3,i,j-2:j+3,k-2:k+2)
      w165(:,:) = Q(4,i,j-2:j+3,k-2:k+2)
      T6(:)     = Q(5,i,j-2:j+3,k) / (R * Q(1,i,j-2:j+3,k))
      u6(:)     = u561(3,:)
      v6(:)     = v561(3,:)
      w6(:)     = w165(:,3)
      mu(:)     = mu6(T6(:))
      uy3(:)    = dx6(u6(:), dy(j))
      vy3(:)    = dx6(v6(:), dy(j))
      wy3(:)    = dx6(w6(:), dy(j))
      vz3(:)    = dy65(v165(:,:), dz(k))
      wz3(:)    = dy65(w165(:,:), dz(k))
      ux3(:)    = dy56(u561(:,:), dx(i))
      vx3(:)    = dy56(v561(:,:), dx(i))
      call tauxy4(mu(:), uy3(:), vx3(:), u6(:), tyx, utyx)
      call tauxx4(mu(:), vy3(:), wz3(:), ux3(:), v6(:), tyy, vtyy)
      call tauxy4(mu(:), vz3(:), wy3(:), w6(:), tyz, wtyz)
      kTy = heat_conduction6(mu(:), T6(:), dy(j))
    else
      T323(:,:,:) = Q(5,i-1:i+1,j:j+1,k-1:k+1) / (R * Q(1,i-1:i+1,j:j+1,k-1:k+1))
      u321(:,:)   = Q(2,i-1:i+1,j:j+1,k)
      v321(:,:)   = Q(3,i-1:i+1,j:j+1,k)
      v123(:,:)   = Q(3,i,j:j+1,k-1:k+1)
      w123(:,:)   = Q(4,i,j:j+1,k-1:k+1)
      Tx(:,:)     = T323(:,:,2)
      Ty(:)       = T323(2,:,2)
      Tz(:,:)     = T323(2,:,:)
      u2(:)       = u321(2,:)
      v2(:)       = v321(2,:)
      w2(:)       = w123(:,2)
      mx(:)       = mu23(Tx(:,:))
      my          = mu2(Ty(:))
      mz(:)       = mu32(Tz(:,:))
      muy         = my * (-u2(1) + u2(2)) * dy(j)
      mvy         = my * (-v2(1) + v2(2)) * dy(j)
      mwy         = my * (-w2(1) + w2(2)) * dy(j)
      mvz         = dy23(mz(:), v123(:,:), dz(k))
      mwz         = dy23(mz(:), w123(:,:), dz(k))
      mux         = dy32(mx(:), u321(:,:), dx(i))
      mvx         = dy32(mx(:), v321(:,:), dx(i))
      tyx         = muy + mvx
      tyy         = 2.d0 * (2.d0 * mvy - mwz - mux) / 3.d0
      tyz         = mvz + mwy
      utyx        = 0.5d0 * (u2(1) + u2(2)) * tyx
      vtyy        = 0.5d0 * (v2(1) + v2(2)) * tyy
      wtyz        = 0.5d0 * (w2(1) + w2(2)) * tyz
      kTy         = Cp * my * (-Ty(1) + Ty(2)) * dy(j) / Pr
    endif
    F(2,i-offset,j-offset+1,k-offset) = F(2,i-offset,j-offset+1,k-offset) - tyx
    F(3,i-offset,j-offset+1,k-offset) = F(3,i-offset,j-offset+1,k-offset) - tyy
    F(4,i-offset,j-offset+1,k-offset) = F(4,i-offset,j-offset+1,k-offset) - tyz
    F(5,i-offset,j-offset+1,k-offset) = F(5,i-offset,j-offset+1,k-offset) - (utyx + vtyy + wtyz + kTy)
  end subroutine calc_Fv
  
  attributes(global) subroutine calc_Fv_LES(nx, ny, nz, dy, dx, dz, Q, mut, qc2, F)
    use calc_sutherland, only : mu6, mu2, mu23, mu32
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(in), device    :: mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device :: F(5,nx-accuracy,ny-accuracy+1,nz-accuracy)
    integer i, j, k
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, muty, H(4), tyxsgs = 0.d0, tyysgs = 0.d0, tyzsgs = 0.d0, Hsgs = 0.d0
    ! 4th-order accuracy
    real(8), dimension(5,6), device :: u561, v561
    real(8), dimension(6,5), device :: v165, w165
    real(8), dimension(6), device   :: T6, u6, v6, w6
    real(8), dimension(3), device   :: uy3, vy3, wy3, vz3, wz3, ux3, vx3, mu
    ! 2nd-order accuracy
    real(8), dimension(3,2,3), device :: T323
    real(8), dimension(3,2), device   :: u321, v321, Tx
    real(8), dimension(2,3), device   :: v123, w123, Tz
    real(8), dimension(2), device     :: Ty, u2, v2, w2, mz, mx
    real(8) my, muy, mvy, mwy, mvz, mwz, mux, mvx
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset - 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      u561(:,:) = Q(2,i-2:i+2,j-2:j+3,k)
      v561(:,:) = Q(3,i-2:i+2,j-2:j+3,k)
      v165(:,:) = Q(3,i,j-2:j+3,k-2:k+2)
      w165(:,:) = Q(4,i,j-2:j+3,k-2:k+2)
      T6(:)     = Q(5,i,j-2:j+3,k) / (R * Q(1,i,j-2:j+3,k))
      u6(:)     = u561(3,:)
      v6(:)     = v561(3,:)
      w6(:)     = w165(:,3)
      mu(:)     = mu6(T6(:))
      uy3(:)    = dx6(u6(:), dy(j))
      vy3(:)    = dx6(v6(:), dy(j))
      wy3(:)    = dx6(w6(:), dy(j))
      vz3(:)    = dy65(v165(:,:), dz(k))
      wz3(:)    = dy65(w165(:,:), dz(k))
      ux3(:)    = dy56(u561(:,:), dx(i))
      vx3(:)    = dy56(v561(:,:), dx(i))
      call tauxy4(mu(:), uy3(:), vx3(:), u6(:), tyx, utyx)
      call tauxx4(mu(:), vy3(:), wz3(:), ux3(:), v6(:), tyy, vtyy)
      call tauxy4(mu(:), vz3(:), wy3(:), w6(:), tyz, wtyz)
      kTy = heat_conduction6(mu(:), T6(:), dy(j))
      muty   = 0.0625d0 * (-mut(i,j-1,k) + 9.d0 * (mut(i,j,k) + mut(i,j+1,k)) -mut(i,j+2,k))
      tyxsgs = muty * (uy3(2) + vx3(2))
      tyysgs = 2.d0 * muty * (2.d0 * vy3(2) - ux3(2) - wz3(2)) / 3.d0
      tyzsgs = muty * (vz3(2) + wy3(2))
      H(:)   = (gamma * Q(5,i,j-1:j+2,k) / (Q(1,i,j-1:j+2,k) * gamma_1)) &
               + 0.5d0 * (Q(2,i,j-1:j+2,k)**2 + Q(3,i,j-1:j+2,k)**2 + Q(4,i,j-1:j+2,k)**2) + qc2(i,j-1:j+2,k)
      Hsgs   = -muty * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) / 3.d0) * dy(j) / Prt
    else
      T323(:,:,:) = Q(5,i-1:i+1,j:j+1,k-1:k+1) / (R * Q(1,i-1:i+1,j:j+1,k-1:k+1))
      u321(:,:)   = Q(2,i-1:i+1,j:j+1,k)
      v321(:,:)   = Q(3,i-1:i+1,j:j+1,k)
      v123(:,:)   = Q(3,i,j:j+1,k-1:k+1)
      w123(:,:)   = Q(4,i,j:j+1,k-1:k+1)
      Tx(:,:)     = T323(:,:,2)
      Ty(:)       = T323(2,:,2)
      Tz(:,:)     = T323(2,:,:)
      u2(:)       = u321(2,:)
      v2(:)       = v321(2,:)
      w2(:)       = w123(:,2)
      mx(:)       = mu23(Tx(:,:))
      my          = mu2(Ty(:))
      mz(:)       = mu32(Tz(:,:))
      muy         = my * (-u2(1) + u2(2)) * dy(j)
      mvy         = my * (-v2(1) + v2(2)) * dy(j)
      mwy         = my * (-w2(1) + w2(2)) * dy(j)
      mvz         = dy23(mz(:), v123(:,:), dz(k))
      mwz         = dy23(mz(:), w123(:,:), dz(k))
      mux         = dy32(mx(:), u321(:,:), dx(i))
      mvx         = dy32(mx(:), v321(:,:), dx(i))
      tyx         = muy + mvx
      tyy         = 2.d0 * (2.d0 * mvy - mwz - mux) / 3.d0
      tyz         = mvz + mwy
      utyx        = 0.5d0 * (u2(1) + u2(2)) * tyx
      vtyy        = 0.5d0 * (v2(1) + v2(2)) * tyy
      wtyz        = 0.5d0 * (w2(1) + w2(2)) * tyz
      kTy         = Cp * my * (-Ty(1) + Ty(2)) * dy(j) / Pr
      my     = 0.5d0 * (mut(i,j,k) + mut(i,j+1,k))
      mz(:)  = (/0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i,j+1,k-1) + mut(i,j+1,k)), &
                 0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i,j+1,k) + mut(i,j+1,k+1))/)
      mx(:)  = (/0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j+1,k) + mut(i,j+1,k)), &
                 0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j+1,k) + mut(i+1,j+1,k))/)
      muy    = my * (-u2(1) + u2(2)) * dy(j)
      mvy    = my * (-v2(1) + v2(2)) * dy(j)
      mwy    = my * (-w2(1) + w2(2)) * dy(j)
      mvz    = dy23(mz(:), v123(:,:), dz(k))
      mwz    = dy23(mz(:), w123(:,:), dz(k))
      mux    = dy32(mx(:), u321(:,:), dx(i))
      mvx    = dy32(mx(:), v321(:,:), dx(i))
      tyxsgs = muy + mvx
      tyysgs = 2.d0 * (2.d0 * mvy - mwz - mux) / 3.d0
      tyzsgs = mvz + mwy
      H(2:3) = (gamma * Q(5,i,j:j+1,k) / (Q(1,i,j:j+1,k) * gamma_1)) &
               + 0.5d0 * (u2(:)**2 + v2(:)**2 + w2(:)**2) + qc2(i,j:j+1,k)
      Hsgs   = -my * (-H(2) + H(3)) * dy(j) / Prt
    endif
    F(2,i-offset,j-offset+1,k-offset) = F(2,i-offset,j-offset+1,k-offset) - (tyx+tyxsgs)
    F(3,i-offset,j-offset+1,k-offset) = F(3,i-offset,j-offset+1,k-offset) - (tyy+tyysgs)
    F(4,i-offset,j-offset+1,k-offset) = F(4,i-offset,j-offset+1,k-offset) - (tyz+tyzsgs)
    F(5,i-offset,j-offset+1,k-offset) = F(5,i-offset,j-offset+1,k-offset) - (utyx + vtyy + wtyz + kTy + Hsgs)
  end subroutine calc_Fv_LES
  
  attributes(global) subroutine calc_Gv(nx, ny, nz, dx, dy, dz, Q, G)
    use calc_sutherland, only : mu6, mu2, mu32
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(inout), device :: G(5,nx-accuracy,ny-accuracy,nz-accuracy+1)
    integer i, j, k
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz
    ! 4th-order accuracy
    real(8), dimension(5,6), device :: u516, w516, v156, w156
    real(8), dimension(6), device   :: T6, u6, v6, w6
    real(8), dimension(3), device   :: uz3, vz3, wz3, wx3, ux3, vy3, wy3, mu
    ! 2nd-order accuracy
    real(8), dimension(3,3,2), device :: T332
    real(8), dimension(3,2), device   :: u312, w312, v132, w132, Tx, Ty
    real(8), dimension(2), device     :: Tz, u2, v2, w2, mx, my
    real(8) mz, muz, mvz, mwz, mwx, mux, mvy, mwy
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset - 1
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      u516(:,:) = Q(2,i-2:i+2,j,k-2:k+3)
      v156(:,:) = Q(3,i,j-2:j+2,k-2:k+3)
      w516(:,:) = Q(4,i-2:i+2,j,k-2:k+3)
      w156(:,:) = Q(4,i,j-2:j+2,k-2:k+3)
      T6(:)     = Q(5,i,j,k-2:k+3) / (R * Q(1,i,j,k-2:k+3))
      u6(:)     = u516(3,:)
      v6(:)     = v156(3,:)
      w6(:)     = w156(3,:)
      mu(:)     = mu6(T6(:))
      uz3(:)    = dx6(u6(:), dz(k))
      vz3(:)    = dx6(v6(:), dz(k))
      wz3(:)    = dx6(w6(:), dz(k))
      wx3(:)    = dy56(w516(:,:), dx(i))
      ux3(:)    = dy56(u516(:,:), dx(i))
      vy3(:)    = dy56(v156(:,:), dy(j))
      wy3(:)    = dy56(w156(:,:), dy(j))
      call tauxy4(mu(:), wx3(:), uz3(:), u6(:), tzx, utzx)
      call tauxy4(mu(:), vz3(:), wy3(:), v6(:), tzy, vtzy)
      call tauxx4(mu(:), wz3(:), ux3(:), vy3(:), w6(:), tzz, wtzz)
      kTz = heat_conduction6(mu(:), T6(:), dz(k))
    else
      T332(:,:,:) = Q(5,i-1:i+1,j-1:j+1,k:k+1) / (R * Q(1,i-1:i+1,j-1:j+1,k:k+1))
      u312(:,:)   = Q(2,i-1:i+1,j,k:k+1)
      v132(:,:)   = Q(3,i,j-1:j+1,k:k+1)
      w312(:,:)   = Q(4,i-1:i+1,j,k:k+1)
      w132(:,:)   = Q(4,i,j-1:j+1,k:k+1)
      Tx(:,:)     = T332(:,2,:)
      Ty(:,:)     = T332(2,:,:)
      Tz(:)       = T332(2,2,:)
      u2(:)       = u312(2,:)
      v2(:)       = v132(2,:)
      w2(:)       = w132(2,:)
      mx(:)       = mu32(Tx(:,:))
      my(:)       = mu32(Ty(:,:))
      mz          = mu2(Tz(:))
      muz         = mz * (-u2(1) + u2(2)) * dz(k)
      mvz         = mz * (-v2(1) + v2(2)) * dz(k)
      mwz         = mz * (-w2(1) + w2(2)) * dz(k)
      mwx         = dy32(mx(:), w312(:,:), dx(i))
      mux         = dy32(mx(:), u312(:,:), dx(i))
      mvy         = dy32(my(:), v132(:,:), dy(j))
      mwy         = dy32(my(:), w132(:,:), dy(j))
      tzx         = mwx + muz
      tzy         = mvz + mwy
      tzz         = 2.d0 * (2.d0 * mwz - mux - mvy) / 3.d0
      utzx        = 0.5d0 * (u2(1) + u2(2)) * tzx
      vtzy        = 0.5d0 * (v2(1) + v2(2)) * tzy
      wtzz        = 0.5d0 * (w2(1) + w2(2)) * tzz
      kTz         = Cp * mz * (-Tz(1) + Tz(2)) * dz(k) / Pr
    endif
    G(2,i-offset,j-offset,k-offset+1) = G(2,i-offset,j-offset,k-offset+1) - tzx
    G(3,i-offset,j-offset,k-offset+1) = G(3,i-offset,j-offset,k-offset+1) - tzy
    G(4,i-offset,j-offset,k-offset+1) = G(4,i-offset,j-offset,k-offset+1) - tzz
    G(5,i-offset,j-offset,k-offset+1) = G(5,i-offset,j-offset,k-offset+1) - (utzx + vtzy + wtzz + kTz)
  end subroutine calc_Gv

  attributes(global) subroutine calc_Gv_LES(nx, ny, nz, dx, dy, dz, Q, mut, qc2, G)
    use calc_sutherland, only : mu6, mu2, mu32
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(in), device    :: mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device :: G(5,nx-accuracy,ny-accuracy,nz-accuracy+1)
    integer i, j, k
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, mutz, H(4), tzxsgs = 0.d0, tzysgs = 0.d0, tzzsgs = 0.d0, Hsgs = 0.d0
    ! 4th-order accuracy
    real(8), dimension(5,6), device :: u516, w516, v156, w156
    real(8), dimension(6), device   :: T6, u6, v6, w6
    real(8), dimension(3), device   :: uz3, vz3, wz3, wx3, ux3, vy3, wy3, mu
    ! 2nd-order accuracy
    real(8), dimension(3,3,2), device :: T332
    real(8), dimension(3,2), device   :: u312, w312, v132, w132, Tx, Ty
    real(8), dimension(2), device     :: Tz, u2, v2, w2, mx, my
    real(8) mz, muz, mvz, mwz, mwx, mux, mvy, mwy
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset - 1
    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      u516(:,:) = Q(2,i-2:i+2,j,k-2:k+3)
      v156(:,:) = Q(3,i,j-2:j+2,k-2:k+3)
      w516(:,:) = Q(4,i-2:i+2,j,k-2:k+3)
      w156(:,:) = Q(4,i,j-2:j+2,k-2:k+3)
      T6(:)     = Q(5,i,j,k-2:k+3) / (R * Q(1,i,j,k-2:k+3))
      u6(:)     = u516(3,:)
      v6(:)     = v156(3,:)
      w6(:)     = w156(3,:)
      mu(:)     = mu6(T6(:))
      uz3(:)    = dx6(u6(:), dz(k))
      vz3(:)    = dx6(v6(:), dz(k))
      wz3(:)    = dx6(w6(:), dz(k))
      wx3(:)    = dy56(w516(:,:), dx(i))
      ux3(:)    = dy56(u516(:,:), dx(i))
      vy3(:)    = dy56(v156(:,:), dy(j))
      wy3(:)    = dy56(w156(:,:), dy(j))
      call tauxy4(mu(:), wx3(:), uz3(:), u6(:), tzx, utzx)
      call tauxy4(mu(:), vz3(:), wy3(:), v6(:), tzy, vtzy)
      call tauxx4(mu(:), wz3(:), ux3(:), vy3(:), w6(:), tzz, wtzz)
      kTz = heat_conduction6(mu(:), T6(:), dz(k))
      mutz   = 0.0625d0 * (-mut(i,j,k-1) + 9.d0 * (mut(i,j,k) + mut(i,j,k+1)) -mut(i,j,k+2))
      tzxsgs = mutz * (wx3(2) + uz3(2))
      tzysgs = mutz * (vz3(2) + wy3(2))
      tzzsgs = 2.d0 * mutz * (2.d0 * wz3(2) - ux3(2) - vy3(2)) / 3.d0
      H(:)   = (gamma * Q(5,i,j,k-1:k+2) / (Q(1,i,j,k-1:k+2) * gamma_1)) &
               + 0.5d0 * (Q(2,i,j,k-1:k+2)**2 + Q(3,i,j,k-1:k+2)**2 + Q(4,i,j,k-1:k+2)**2) + qc2(i,j,k-1:k+2)
      Hsgs   = -mutz * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) / 3.d0) * dz(k) / Prt
    else
      T332(:,:,:) = Q(5,i-1:i+1,j-1:j+1,k:k+1) / (R * Q(1,i-1:i+1,j-1:j+1,k:k+1))
      u312(:,:)   = Q(2,i-1:i+1,j,k:k+1)
      v132(:,:)   = Q(3,i,j-1:j+1,k:k+1)
      w312(:,:)   = Q(4,i-1:i+1,j,k:k+1)
      w132(:,:)   = Q(4,i,j-1:j+1,k:k+1)
      Tx(:,:)     = T332(:,2,:)
      Ty(:,:)     = T332(2,:,:)
      Tz(:)       = T332(2,2,:)
      u2(:)       = u312(2,:)
      v2(:)       = v132(2,:)
      w2(:)       = w132(2,:)
      mx(:)       = mu32(Tx(:,:))
      my(:)       = mu32(Ty(:,:))
      mz          = mu2(Tz(:))
      muz         = mz * (-u2(1) + u2(2)) * dz(k)
      mvz         = mz * (-v2(1) + v2(2)) * dz(k)
      mwz         = mz * (-w2(1) + w2(2)) * dz(k)
      mwx         = dy32(mx(:), w312(:,:), dx(i))
      mux         = dy32(mx(:), u312(:,:), dx(i))
      mvy         = dy32(my(:), v132(:,:), dy(j))
      mwy         = dy32(my(:), w132(:,:), dy(j))
      tzx         = mwx + muz
      tzy         = mvz + mwy
      tzz         = 2.d0 * (2.d0 * mwz - mux - mvy) / 3.d0
      utzx        = 0.5d0 * (u2(1) + u2(2)) * tzx
      vtzy        = 0.5d0 * (v2(1) + v2(2)) * tzy
      wtzz        = 0.5d0 * (w2(1) + w2(2)) * tzz
      kTz         = Cp * mz * (-Tz(1) + Tz(2)) * dz(k) / Pr
      mz     = 0.5d0 * (mut(i,j,k) + mut(i,j,k+1))
      mx(:)  = (/0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j,k+1) + mut(i,j,k+1)), &
                 0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j,k+1) + mut(i+1,j,k+1))/)
      my(:)  = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i,j-1,k+1) + mut(i,j,k+1)), &
                 0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i,j,k+1) + mut(i,j+1,k+1))/)
      muz    = mz * (-u2(1) + u2(2)) * dz(k)
      mvz    = mz * (-v2(1) + v2(2)) * dz(k)
      mwz    = mz * (-w2(1) + w2(2)) * dz(k)
      mwx    = dy32(mx(:), w312(:,:), dx(i))
      mux    = dy32(mx(:), u312(:,:), dx(i))
      mvy    = dy32(my(:), v132(:,:), dy(j))
      mwy    = dy32(my(:), w132(:,:), dy(j))
      tzxsgs = mwx + muz
      tzysgs = mvz + mwy
      tzzsgs = 2.d0 * (2.d0 * mwz - mux - mvy) / 3.d0
      H(2:3) = (gamma * Q(5,i,j,k:k+1) / (Q(1,i,j,k:k+1) * gamma_1)) &
               + 0.5d0 * (u2(:)**2 + v2(:)**2 + w2(:)**2) + qc2(i,j,k:k+1)
      Hsgs   = -mz * (-H(2) + H(3)) * dz(k) / Prt
    endif
    G(2,i-offset,j-offset,k-offset+1) = G(2,i-offset,j-offset,k-offset+1) - (tzx+tzxsgs)
    G(3,i-offset,j-offset,k-offset+1) = G(3,i-offset,j-offset,k-offset+1) - (tzy+tzysgs)
    G(4,i-offset,j-offset,k-offset+1) = G(4,i-offset,j-offset,k-offset+1) - (tzz+tzzsgs)
    G(5,i-offset,j-offset,k-offset+1) = G(5,i-offset,j-offset,k-offset+1) - (utzx + vtzy + wtzz + kTz + Hsgs)
  end subroutine calc_Gv_LES
end module calc_visc

