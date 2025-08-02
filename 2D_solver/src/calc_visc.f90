module calc_visc
  use mod_globals, only : accuracy, offset, id_visc, id_turbulence, id_av, gamma, R, Pr, Prt
  implicit none
contains
  attributes(device) function interpolation6(a) result(ans)
    real(8), intent(in), device :: a(6)
    real(8) ans(3)
    ans(:) = 0.0625d0 * (-a(1:3) + 9.d0 * (a(2:4) + a(3:5)) -a(4:6))
  end function interpolation6

  attributes(device) function dx6(a, dx) result(ans)
    real(8), intent(in), device :: a(6)
    real(8), intent(in), value  :: dx
    real(8) ans(3)
    ans(:) = 0.125d0 * (9.d0 * (-a(2:4) + a(3:5)) - (-a(1:3) + a(4:6)) / 3.d0) * dx
  end function dx6

  attributes(device) function dy5(a, dy) result(ans)
    real(8), intent(in), device :: a(5)
    real(8), intent(in), value  :: dy
    real(8) ans
    ans = (2.d0 * (-a(2) + a(4)) - 0.25d0 * (-a(1) + a(5))) * dy / 3.d0
  end function dy5

  attributes(device) function dy23(mu, a, dy) result(ans)
    real(8), intent(in), device :: mu(2), a(2,3)
    real(8), intent(in), value  :: dy
    real(8) ans
    ans = 0.25d0 * (mu(1) * (-a(1,1) + a(1,2) -a(2,1) + a(2,2)) &
                  + mu(2) * (-a(1,2) + a(1,3) -a(2,2) + a(2,3)))
  end function dy23

  attributes(device) function dy65(a, dy) result(ans)
    real(8), intent(in), device :: a(6,5)
    real(8), intent(in), value  :: dy
    real(8) ans(3), a1(5), a2(5), a3(5), a4(5), a5(5), a6(5), ay(6)
    a1 = a(1,:); a2 = a(2,:); a3 = a(3,:); a4 = a(4,:); a5 = a(5,:); a6 = a(6,:)
    ay(1) = dy5(a1(:), dy); ay(2) = dy5(a2(:), dy); ay(3) = dy5(a3(:), dy)
    ay(4) = dy5(a4(:), dy); ay(5) = dy5(a5(:), dy); ay(6) = dy5(a6(:), dy)
    ans = interpolation6(ay(:))
  end function dy65
  
  attributes(device) function dy32(mu, a, dy) result(ans)
    real(8), intent(in), device :: mu(2), a(3,2)
    real(8), intent(in), value  :: dy
    real(8) ans
    ans = 0.25d0 * (mu(1) * (-a(1,1) + a(2,1) - a(1,2) + a(2,2)) &
                  + mu(2) * (-a(2,1) + a(3,1) - a(2,2) + a(3,2)))
  end function dy32

  attributes(device) function dy56(a, dy) result(ans)
    real(8), intent(in), device :: a(5,6)
    real(8), intent(in), value  :: dy
    real(8) ans(3), a1(5), a2(5), a3(5), a4(5), a5(5), a6(5), ay(6)
    a1 = a(:,1); a2 = a(:,2); a3 = a(:,3); a4 = a(:,4); a5 = a(:,5); a6 = a(:,6)
    ay(1) = dy5(a1(:), dy); ay(2) = dy5(a2(:), dy); ay(3) = dy5(a3(:), dy)
    ay(4) = dy5(a4(:), dy); ay(5) = dy5(a5(:), dy); ay(6) = dy5(a6(:), dy)
    ans = interpolation6(ay(:))
  end function dy56

  attributes(device) function flux4(a) result(ans)
    real(8), intent(in), device :: a(3)
    real(8), intent(out)        :: ans
    ans = 0.125d0 * ((9.d0 - 1.d0 / 3.d0) * a(2) - (a(1) + a(3)) / 3.d0)
  end function flux4

  attributes(device) subroutine tauxx4(mu, ux, vy, u6, txx, utxx)
    real(8), intent(in), dimension(3), device :: mu, ux, vy
    real(8), intent(in), dimension(6), device :: u6
    real(8), intent(out) :: txx, utxx
    real(8), dimension(3) :: tau, utau
    tau(:)  = 2.d0 * mu(:) * (2.d0 * ux(:) - vy(:)) / 3.d0
    utau(:) = interpolation6(u6(:)) * tau(:)
    txx     = flux4(tau(:))
    utxx    = flux4(utau(:))
  end subroutine tauxx4

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

  attributes(device) function heat_conduction6(mu, T, dx) result(ans)
    real(8), intent(in), device :: mu(3), T(6)
    real(8), intent(in), value  :: dx
    real(8), dimension(3) ::  kTx
    real(8) :: ans, Cp = gamma * R / (gamma - 1.d0)
    kTx(:) = Cp * mu(:) * dx6(T(:), dx) / Pr
    ans = flux4(kTx(:))
  end function heat_conduction6

  attributes(global) subroutine calc_Ev(nx, ny, dx, dy, rho, u, v, p, E)
    use calc_sutherland, only : mu6, mu2, mu23
    integer, intent(in), value                    :: nx, ny
    real(8), intent(in), dimension(nx-1), device  :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1), device  :: dy ! 1 / dy
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(inout), device                :: E(4,nx-accuracy+1,ny-accuracy)
    integer i, j
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    real(8) :: txx, txy, utxx, vtxy, kTx
    ! 4th-order accuracy
    real(8), dimension(6,5), device :: u65, v65
    real(8), dimension(6), device   :: T6, u6, v6
    real(8), dimension(3), device   :: ux3, vx3, uy3, vy3, mu
    ! 2nd-order accuracy
    real(8), dimension(2,3), device :: T23, u23, v23, Ty
    real(8), dimension(2), device   :: Tx, u2, v2, my
    real(8) mx, mux, mvx, muy, mvy
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset - 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset

    if (id_visc ==2 .and. 3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2) then
      u65(:,:) = u(i-2:i+3,j-2:j+2)
      v65(:,:) = v(i-2:i+3,j-2:j+2)
      T6(:)    = p(i-2:i+3,j) / (R * rho(i-2:i+3,j))
      u6(:)    = u65(:,3)
      v6(:)    = v65(:,3)
      mu(:)    = mu6(T6(:))
      ux3(:)   = dx6(u6(:), dx(i))
      vx3(:)   = dx6(v6(:), dx(i))
      uy3(:)   = dy65(u65(:,:), dy(j))
      vy3(:)   = dy65(v65(:,:), dy(j))
      call tauxx4(mu(:), ux3(:), vy3(:), u6(:), txx, utxx)
      call tauxy4(mu(:), uy3(:), vx3(:), v6(:), txy, vtxy)
      kTx = heat_conduction6(mu(:), T6(:), dx(i))
    else
      T23(:,:) = p(i:i+1,j-1:j+1) / (R * rho(i:i+1,j-1:j+1))
      u23(:,:) = u(i:i+1,j-1:j+1)
      v23(:,:) = v(i:i+1,j-1:j+1)
      Tx(:)    = T23(:,2)
      Ty(:,:)  = T23(:,:)
      u2(:)    = u23(:,2)
      v2(:)    = v23(:,2)
      mx       = mu2(Tx(:))
      my(:)    = mu23(Ty(:,:))
      mux      = mx * (-u2(1) + u2(2)) * dx(i)
      mvx      = mx * (-v2(1) + v2(2)) * dx(i)
      muy      = dy23(my(:), u23(:,:), dy(j))
      mvy      = dy23(my(:), v23(:,:), dy(j))
      txx      = 2.d0 * (2.d0 * mux - mvy) / 3.d0
      txy      = muy + mvx
      utxx     = 0.5d0 * (u2(1) + u2(2)) * txx
      vtxy     = 0.5d0 * (v2(1) + v2(2)) * txy
      kTx      = Cp * mx * (-Tx(1) + Tx(2)) * dx(i) / Pr
    endif

    E(2,i-offset+1,j-offset) = E(2,i-offset+1,j-offset) - txx
    E(3,i-offset+1,j-offset) = E(3,i-offset+1,j-offset) - txy
    E(4,i-offset+1,j-offset) = E(4,i-offset+1,j-offset) - (utxx + vtxy + kTx)
  end subroutine calc_Ev
  
  attributes(global) subroutine calc_Fv(nx, ny, dy, dx, rho, u, v, p, F)
    use calc_sutherland, only : mu6, mu2, mu23, mu32
    integer, intent(in), value                    :: nx, ny
    real(8), intent(in), dimension(ny-1), device  :: dy ! 1 / dy
    real(8), intent(in), dimension(nx-1), device  :: dx ! 1 / dx
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(inout), device                :: F(4,nx-accuracy,ny-accuracy+1)
    integer i, j
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    real(8) :: tyx, tyy, utyx, vtyy, kTy
    ! 4th-order accuracy
    real(8), dimension(5,6), device :: u56, v56
    real(8), dimension(6), device   :: T6, u6, v6
    real(8), dimension(3), device   :: uy3, vy3, ux3, vx3, mu
    ! 2nd-order accuracy
    real(8), dimension(3,2), device :: T32, u32, v32, Tx
    real(8), dimension(2), device   :: Ty, u2, v2, mx
    real(8) my, muy, mvy, mux, mvx
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset - 1

    if (id_visc == 2 .and. 3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3) then
      u56(:,:) = u(i-2:i+2,j-2:j+3)
      v56(:,:) = v(i-2:i+2,j-2:j+3)
      T6(:)    = p(i,j-2:j+3) / (R * rho(i,j-2:j+3))
      u6(:)    = u56(3,:)
      v6(:)    = v56(3,:)
      mu(:)    = mu6(T6(:))
      uy3(:)   = dx6(u6(:), dy(j))
      vy3(:)   = dx6(v6(:), dy(j))
      ux3(:)   = dy56(u56(:,:), dx(i))
      vx3(:)   = dy56(v56(:,:), dx(i))
      call tauxy4(mu(:), uy3(:), vx3(:), u6(:), tyx, utyx)
      call tauxx4(mu(:), vy3(:), ux3(:), v6(:), tyy, vtyy)
      kTy = heat_conduction6(mu(:), T6(:), dy(j))
    else
      T32(:,:) = p(i-1:i+1,j:j+1) / (R * rho(i-1:i+1,j:j+1))
      u32(:,:) = u(i-1:i+1,j:j+1)
      v32(:,:) = v(i-1:i+1,j:j+1)
      Tx(:,:)  = T32(:,:)
      Ty(:)    = T32(2,:)
      u2(:)    = u32(2,:)
      v2(:)    = v32(2,:)
      mx(:)    = mu23(Tx(:,:))
      my       = mu2(Ty(:))
      muy      = my * (-u2(1) + u2(2)) * dy(j)
      mvy      = my * (-v2(1) + v2(2)) * dy(j)
      mux      = dy32(mx(:), u32(:,:), dx(i))
      mvx      = dy32(mx(:), v32(:,:), dx(i))
      tyx      = muy + mvx
      tyy      = 2.d0 * (2.d0 * mvy - mux) / 3.d0
      utyx     = 0.5d0 * (u2(1) + u2(2)) * tyx
      vtyy     = 0.5d0 * (v2(1) + v2(2)) * tyy
      kTy      = Cp * my * (-Ty(1) + Ty(2)) * dy(j) / Pr
    endif

    F(2,i-offset,j-offset+1) = F(2,i-offset,j-offset+1) - tyx
    F(3,i-offset,j-offset+1) = F(3,i-offset,j-offset+1) - tyy
    F(4,i-offset,j-offset+1) = F(4,i-offset,j-offset+1) - (utyx + vtyy + kTy)
  end subroutine calc_Fv
end module calc_visc

