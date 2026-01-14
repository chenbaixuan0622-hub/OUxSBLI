module calc_visc4
  use curand
  use curand_device
  use mod_globals, only : id_visc, gamma, R, Pr, Prt, dt, threadsEv, threadsFv
  use mod_constant, only : Cp, gamma_1, Cp_over_Pr, one_third, two_third, one_twelfth
  use calc_rand
  implicit none
  real(8), parameter :: one_24 = 1.d0 / 24.d0
contains
  !dir$ inline
  attributes(device) function flux4(a) result(ans)
    real(8), intent(in), device :: a(3)
    real(8) ans
    ans = (-a(1) + 26.d0 * a(2) - a(3)) * one_24
  end function flux4


  !$dir inline
  attributes(device) subroutine calc_tau_straight(mu, u, vy, d, t11, ut11)
    real(8), intent(in)        :: mu(3), u(6), vy(6)
    real(8), intent(in), value :: d
    real(8), intent(out)       :: t11, ut11
    real(8) tmp(3)
    tmp(:) = (2.25d0 * (-u(2:4) + u(3:5)) - (-u(1:3) + u(4:6)) * one_twelfth) * d
    tmp(:) = tmp(:) - 0.0625d0 * (-vy(1:3) + 9.d0 * (vy(2:4) + vy(3:5)) - vy(4:6))
    tmp(:) = two_third * mu(:) * tmp(:)
    t11    = (-tmp(1) + 26.d0 * tmp(2) - tmp(3)) * one_24
    tmp(1) = 0.0625d0 * (-u(1) + 9.d0 * (u(2) + u(3)) - u(4)) * tmp(1)
    tmp(2) = 0.0625d0 * (-u(2) + 9.d0 * (u(3) + u(4)) - u(5)) * tmp(2)
    tmp(3) = 0.0625d0 * (-u(3) + 9.d0 * (u(4) + u(5)) - u(6)) * tmp(3)
    ut11   = (-tmp(1) + 26.d0 * tmp(2) - tmp(3)) * one_24
  end subroutine calc_tau_straight


  !$dir inline
  attributes(device) subroutine calc_tau_cross(mu, v, uy, d, t12, vt12)
    real(8), intent(in)        :: mu(3), v(6), uy(6)
    real(8), intent(in), value :: d
    real(8), intent(out)       :: t12, vt12
    real(8) tmp(3)
    tmp(:) = (1.125d0 * (-v(2:4) + v(3:5)) - (-v(1:3) + v(4:6)) * one_24) * d
    tmp(:) = tmp(:) + 0.0625d0 * (-uy(1:3) + 9.d0 * (uy(2:4) + uy(3:5)) - uy(4:6))
    tmp(:) = mu(:) * tmp(:)
    t12    = (-tmp(1) + 26.d0 * tmp(2) - tmp(3)) * one_24
    tmp(1) = 0.0625d0 * (-v(1) + 9.d0 * (v(2) + v(3)) - v(4)) * tmp(1)
    tmp(2) = 0.0625d0 * (-v(2) + 9.d0 * (v(3) + v(4)) - v(5)) * tmp(2)
    tmp(3) = 0.0625d0 * (-v(3) + 9.d0 * (v(4) + v(5)) - v(6)) * tmp(3)
    vt12   = (-tmp(1) + 26.d0 * tmp(2) - tmp(3)) * one_24
  end subroutine calc_tau_cross


  attributes(global) subroutine calc_Ev4(nx, ny, dx, dy, Q, T, mu, E, state)
    integer, intent(in), value     :: nx, ny
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: Q(4,nx,ny), T(nx,ny), mu(nx,ny)
    real(8), intent(inout), device :: E(4,nx-1,ny-2)
    type(curandStateXORWOW), intent(inout), device, optional :: state(nx*ny)
    real(8), shared ::  u(-2:threadsEv%x+3,threadsEv%y)
    real(8), shared ::  v(-2:threadsEv%x+3,threadsEv%y)
    real(8), shared :: uy(-2:threadsEv%x+3,threadsEv%y)
    real(8), shared :: vy(-2:threadsEv%x+3,threadsEv%y)
    integer i, j, it, jt, ii, i_base
    real(8) :: txx, txy, utxx, vtxy, kTx
    it = threadIdx%x
    jt = threadIdx%y
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    do ii = it-2, threadsEv%x+3, blockDim%x
      i = i_base + ii
      if (1 <= i .and. i <= nx .and. j <= ny) then
        u(ii,jt) = Q(2,i,j)
        v(ii,jt) = Q(3,i,j)
      endif
      if (1 <= i .and. i <= nx .and. 3 <= j .and. j <= ny-2) then
        uy(ii,jt) = (two_third * (-Q(2,i,j-1) + Q(2,i,j+1)) - one_twelfth * (-Q(2,i,j-2) + Q(2,i,j+2))) * dy(j)
        vy(ii,jt) = (two_third * (-Q(3,i,j-1) + Q(3,i,j+1)) - one_twelfth * (-Q(3,i,j-2) + Q(3,i,j+2))) * dy(j)
      endif
    enddo
    call syncthreads()
    i  = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j) return
    if (3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2) then
      block
        real(8), device :: mu3(3)
        mu3(:) = 0.0625d0 * (9.d0 * (mu(i-1:i+1,j) + mu(i:i+2,j)) - (mu(i-2:i,j) + mu(i+1:i+3,j)))
        block ! dQdx
          real(8), device :: kTx3(3)
          kTx3(:) = Cp_over_Pr * mu3(:) * &
                    (1.125d0 * (-T(i-1:i+1,j) + T(i:i+2,j)) - (-T(i-2:i,j) + T(i+1:i+3,j)) * one_24) * dx(i)
          kTx     = flux4(kTx3(:))
        end block
        call calc_tau_straight(mu3, u(it-2:it+3,jt), vy(it-2:it+3,jt), dx(i), txx, utxx)
        call calc_tau_cross(mu3, v(it-2:it+3,jt), uy(it-2:it+3,jt), dx(i), txy, vtxy)
      end block
    else
      block
        real(8) mx, mux, mvx, muy, mvy 
        mx  = 0.5d0 * (mu(i,j) + mu(i+1,j))
        kTx = Cp_over_Pr * mx * (-T(i,j) + T(i+1,j)) * dx(i)
        block
          real(8), device :: my(2)
          my(1) = 0.25d0 * (mu(i,j-1) + mu(i,j  ) + mu(i+1,j-1) + mu(i+1,j  ))
          my(2) = 0.25d0 * (mu(i,j  ) + mu(i,j+1) + mu(i+1,j  ) + mu(i+1,j+1))
          muy = 0.25d0 * (my(1) * (-Q(2,i,j-1) + Q(2,i,j) - Q(2,i+1,j-1) + Q(2,i+1,j)) &
                        + my(2) * (-Q(2,i,j) + Q(2,i,j+1) - Q(2,i+1,j) + Q(2,i+1,j+1))) * dy(j)
          mvy = 0.25d0 * (my(1) * (-Q(3,i,j-1) + Q(3,i,j) - Q(3,i+1,j-1) + Q(3,i+1,j)) &
                        + my(2) * (-Q(3,i,j) + Q(3,i,j+1) - Q(3,i+1,j) + Q(3,i+1,j+1))) * dy(j)
        end block
        mux  = mx * (-u(it,jt) + u(it+1,jt)) * dx(i)
        mvx  = mx * (-v(it,jt) + v(it+1,jt)) * dx(i)
        txx  = two_third * (2.d0 * mux - mvy)
        txy  = muy + mvx
        utxx = 0.5d0 * (u(it,jt) + u(it+1,jt)) * txx
        vtxy = 0.5d0 * (v(it,jt) + v(it+1,jt)) * txy
      end block
    endif
    E(2,i,j-1) = E(2,i,j-1) - txx
    E(3,i,j-1) = E(3,i,j-1) - txy
    E(4,i,j-1) = E(4,i,j-1) - (utxx + vtxy + kTx)
  end subroutine calc_Ev4
 

  attributes(global) subroutine calc_Fv4(nx, ny, dy, dx, Q, T, mu, F, state)
    integer, intent(in), value     :: nx, ny
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: Q(4,nx,ny), T(nx,ny), mu(nx,ny)
    real(8), intent(inout), device :: F(4,nx-2,ny-1)
    type(curandStateXORWOW), intent(inout), device, optional :: state(nx*ny)
    real(8), shared ::  u(-2:threadsFv%y+3,threadsFv%x)
    real(8), shared ::  v(-2:threadsFv%y+3,threadsFv%x)
    real(8), shared :: ux(-2:threadsFv%y+3,threadsFv%x)
    real(8), shared :: vx(-2:threadsFv%y+3,threadsFv%x)
    integer i, j, it, jt, jj, j_base
    real(8) :: tyx, tyy, utyx, vtyy, kTy
    it = threadIdx%x
    jt = threadIdx%y
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j_base = (blockIdx%y-1)*blockDim%y
    do jj = jt-2, threadsFv%y+3, blockDim%y
      j = j_base + jj
      if (i <= nx .and. 1 <= j .and. j <= ny) then
        u(jj,it) = Q(2,i,j)
        v(jj,it) = Q(3,i,j)
      endif
      if (3 <= i .and. i <= nx-2 .and. 1 <= j .and. j <= ny) then
        ux(jj,it) = (two_third * (-Q(2,i-1,j) + Q(2,i+1,j)) - one_twelfth * (-Q(2,i-2,j) + Q(2,i+2,j))) * dx(i)
        vx(jj,it) = (two_third * (-Q(3,i-1,j) + Q(3,i+1,j)) - one_twelfth * (-Q(3,i-2,j) + Q(3,i+2,j))) * dx(i)
      endif
    enddo
    call syncthreads()
    j  = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j) return
    if (3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3) then
      block
        real(8), device :: mu3(3)
        mu3(:) = 0.0625d0 * (9.d0 * (mu(i,j-1:j+1) + mu(i,j:j+2)) - (mu(i,j-2:j) + mu(i,j+1:j+3)))
        block ! dQdy
          real(8), device :: kTy3(3)
          kTy3(:) = Cp_over_Pr * mu3(:) * &
                    (1.125d0 * (-T(i,j-1:j+1) + T(i,j:j+2)) - (-T(i,j-2:j) + T(i,j+1:j+3)) * one_24) * dy(j)
          kTy     = flux4(kTy3(:))
        end block
        call calc_tau_straight(mu3, v(jt-2:jt+3,it), ux(jt-2:jt+3,it), dy(j), tyy, vtyy)
        call calc_tau_cross(mu3, u(jt-2:jt+3,it), vx(jt-2:jt+3,it), dy(j), tyx, utyx)
      end block
    else
      block
        real(8) my, muy, mvy, mux, mvx
        my  = 0.5d0 * (mu(i,j) + mu(i,j+1))
        kTy = Cp_over_Pr * my * (-T(i,j) + T(i,j+1)) * dy(j)
        block
          real(8), device :: mx(2)
          mx(1) = 0.25d0 * (mu(i-1,j) + mu(i,  j) + mu(i-1,j+1) + mu(i,  j+1))
          mx(2) = 0.25d0 * (mu(i,  j) + mu(i+1,j) + mu(i,  j+1) + mu(i+1,j+1))
          mux = 0.25d0 * (mx(1) * (-Q(2,i-1,j) + Q(2,i,j) - Q(2,i-1,j+1) + Q(2,i,j+1)) &
                        + mx(2) * (-Q(2,i,j) + Q(2,i+1,j) - Q(2,i,j+1) + Q(2,i+1,j+1))) * dx(i)
          mvx = 0.25d0 * (mx(1) * (-Q(3,i-1,j) + Q(3,i,j) - Q(3,i-1,j+1) + Q(3,i,j+1)) &
                        + mx(2) * (-Q(3,i,j) + Q(3,i+1,j) - Q(3,i,j+1) + Q(3,i+1,j+1))) * dx(i)
        end block
        my  = 0.5d0 * (mu(i,j) + mu(i,j+1))
        kTy = Cp_over_Pr * my * (-T(i,j) + T(i,j+1)) * dy(j)
        muy  = my * (-u(jt,it) + u(jt+1,it)) * dy(j)
        mvy  = my * (-v(jt,it) + v(jt+1,it)) * dy(j)
        tyx  = muy + mvx
        tyy  = two_third * (2.d0 * mvy - mux)
        utyx = 0.5d0 * (u(jt,it) + u(jt+1,it)) * tyx
        vtyy = 0.5d0 * (v(jt,it) + v(jt+1,it)) * tyy
      end block
    endif
    F(2,i-1,j) = F(2,i-1,j) - tyx
    F(3,i-1,j) = F(3,i-1,j) - tyy
    F(4,i-1,j) = F(4,i-1,j) - (utyx + vtyy + kTy)
  end subroutine calc_Fv4
end module calc_visc4

