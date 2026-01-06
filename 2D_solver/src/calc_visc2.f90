module calc_visc2
  use curand
  use curand_device
  use mod_globals, only : id_LL, gamma, R, Pr, Prt, dt, threadsEv, threadsFv
  use mod_constant, only : Cp, gamma_1, Cp_over_Pr, one_third, two_third
  use calc_rand
  implicit none
contains
  attributes(global) subroutine calc_Ev2(nx, ny, dx, dy, Q, T, mu, E, state)
    integer, intent(in), value     :: nx, ny
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: Q(4,nx,ny), T(nx,ny), mu(nx,ny)
    real(8), intent(inout), device :: E(4,nx-1,ny-2)
    type(curandStateXORWOW), intent(inout), device, optional :: state(nx*ny)
    real(8), shared :: u(threadsEv%x+1,0:threadsEv%y+1)
    real(8), shared :: v(threadsEv%x+1,0:threadsEv%y+1)
    integer i, j, it, jt
    real(8) :: txx, txy, utxx, vtxy, kTx
    real(8) m1, m2, mx, mux, mvx, muy, mvy 
    it = threadIdx%x
    jt = threadIdx%y
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    if (nx-1 < i .or. ny-1 < j) return
    u(it,jt-1:jt+1) = Q(2,i,j-1:j+1)
    v(it,jt-1:jt+1) = Q(3,i,j-1:j+1)
    if (it == blockDim%x) then
      u(it+1,jt-1:jt+1) = Q(2,i+1,j-1:j+1)
      v(it+1,jt-1:jt+1) = Q(3,i+1,j-1:j+1)
    endif 
    call syncthreads()

    mx  = 0.5d0 * (mu(i,j) + mu(i+1,j))
    kTx = Cp_over_Pr * mx * (-T(i,j) + T(i+1,j)) * dx(i)
    block
      real(8) my1, my2
      my1 = 0.25d0 * (mu(i,j-1) + mu(i,j  ) + mu(i+1,j-1) + mu(i+1,j  ))
      my2 = 0.25d0 * (mu(i,j  ) + mu(i,j+1) + mu(i+1,j  ) + mu(i+1,j+1))
      muy = 0.25d0 * (my1 * (-u(it,jt-1) - u(it+1,jt-1)) + (my1 - my2) * (u(it,jt) + u(it+1,jt)) &
                    + my2 * ( u(it,jt+1) + u(it+1,jt+1))) * dy(j)
      mvy = 0.25d0 * (my1 * (-v(it,jt-1) - v(it+1,jt-1)) + (my1 - my2) * (v(it,jt) + v(it+1,jt)) &
                    + my2 * ( v(it,jt+1) + v(it+1,jt+1))) * dy(j)
    end block
    mux = mx * (-u(it,jt) + u(it+1,jt)) * dx(i)
    mvx = mx * (-v(it,jt) + v(it+1,jt)) * dx(i)
    txx = two_third * (2.d0 * mux - mvy)
    txy = muy + mvx
    if (kind(id_LL) == 4) then
      block
        real(8) std_t, std_q, over_V, Zq
        real(8), device    :: Z(3), Zx(3)
        real(8), parameter :: kb_over_dt = 1.380649d-23 / dt
        over_V = dx(i) * dy(j)
        std_t  = sqrt(kb_over_dt * over_V * mx * (T(i,j) + T(i+1,j)))
        std_q  = sqrt(kb_over_dt * over_V * mx * Cp_over_Pr * (T(i,j)*T(i,j) + T(i+1,j)*T(i+1,j)))
        !Z   = Z_tilde(state(i+nx*(j-1)))
        !Zx  = Z_tilde(state(i+1+nx*(j-1)))
        Z   = 0.5d0 * (Z + Zx)
        !Zq  = Zq_x(state(i+nx*(j-1)))
        txx = txx + std_t * (2.d0 * Z(1) - Z(3)) * one_third
        txy = txy + std_t * Z(2)
        kTx = kTx + std_q * Zq
      end block
    endif
    utxx = 0.5d0 * (u(it,jt) + u(it+1,jt)) * txx
    vtxy = 0.5d0 * (v(it,jt) + v(it+1,jt)) * txy
    E(2,i,j-1) = E(2,i,j-1) - txx
    E(3,i,j-1) = E(3,i,j-1) - txy
    E(4,i,j-1) = E(4,i,j-1) - (utxx + vtxy + kTx)
  end subroutine calc_Ev2
 

  attributes(global) subroutine calc_Fv2(nx, ny, dy, dx, Q, T, mu, F, state)
    integer, intent(in), value     :: nx, ny
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: Q(4,nx,ny), T(nx,ny), mu(nx,ny)
    real(8), intent(inout), device :: F(4,nx-2,ny-1)
    type(curandStateXORWOW), intent(inout), device, optional :: state(nx*ny)
    real(8), shared :: u(threadsFv%y+1,0:threadsFv%x+1)
    real(8), shared :: v(threadsFv%y+1,0:threadsFv%x+1)
    integer i, j, it, jt
    real(8) :: tyx, tyy, utyx, vtyy, kTy
    real(8) m1, m2, my, muy, mvy, mux, mvx
    it = threadIdx%x
    jt = threadIdx%y
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j) return
    u(jt,it-1:it+1) = Q(2,i-1:i+1,j)
    v(jt,it-1:it+1) = Q(3,i-1:i+1,j)
    if (jt == blockDim%y) then
      u(jt+1,it-1:it+1) = Q(2,i-1:i+1,j+1)
      v(jt+1,it-1:it+1) = Q(3,i-1:i+1,j+1)
    endif 
    call syncthreads()

    block
      real(8) mx1, mx2
      mx1 = 0.25d0 * (mu(i-1,j) + mu(i,  j) + mu(i-1,j+1) + mu(i,  j+1))
      mx2 = 0.25d0 * (mu(i,  j) + mu(i+1,j) + mu(i,  j+1) + mu(i+1,j+1))
      mux = 0.25d0 * (mx1 * (-u(jt,it-1) - u(jt+1,it-1)) + (mx1 - mx2) * (u(jt,it) + u(jt+1,it)) &
                    + mx2 * ( u(jt,it+1) + u(jt+1,it+1))) * dx(i)
      mvx = 0.25d0 * (mx1 * (-v(jt,it-1) - v(jt+1,it-1)) + (mx1 - mx2) * (v(jt,it) + v(jt+1,it)) &
                    + mx2 * ( v(jt,it+1) + v(jt+1,it+1))) * dx(i)
    end block
    my  = 0.5d0 * (mu(i,j) + mu(i,j+1))
    kTy = Cp_over_Pr * my * (-T(i,j) + T(i,j+1)) * dy(j)
    muy = my * (-u(jt,it) + u(jt+1,it)) * dy(j)
    mvy = my * (-v(jt,it) + v(jt+1,it)) * dy(j)
    tyx = muy + mvx
    tyy = two_third * (2.d0 * mvy - mux)
    if (kind(id_LL) == 4) then
      block
        real(8) std_t, std_q, over_V, Zq
        real(8), device    :: Z(3), Zy(3)
        real(8), parameter :: kb_over_dt = 1.380649d-23 / dt
        over_V = dx(i) * dy(j)
        std_t  = sqrt(kb_over_dt * over_V * my * (T(i,j) + T(i,j+1)))
        std_q  = sqrt(kb_over_dt * over_V * my * Cp_over_Pr * (T(i,j)*T(i,j) + T(i,j+1)*T(i,j+1)))
        !Z   = Z_tilde(state(i+nx*(j-1)))
        !Zy  = Z_tilde(state(i+nx*j))
        Z   = 0.5d0 * (Z + Zy)
        !Zq  = Zq_y(state(i+nx*(j-1)))
        tyx = tyx + std_t * Z(2)
        tyy = tyy + std_t * (2.d0 * Z(3) - Z(1)) * one_third
        kTy = kTy + std_q * Zq
      end block
    endif
    utyx = 0.5d0 * (u(jt,it) + u(jt+1,it)) * tyx
    vtyy = 0.5d0 * (v(jt,it) + v(jt+1,it)) * tyy
    F(2,i-1,j) = F(2,i-1,j) - tyx
    F(3,i-1,j) = F(3,i-1,j) - tyy
    F(4,i-1,j) = F(4,i-1,j) - (utyx + vtyy + kTy)
  end subroutine calc_Fv2
end module calc_visc2

