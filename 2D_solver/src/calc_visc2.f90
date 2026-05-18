module calc_visc2
  use mod_globals, only : gamma, R, Pr, Prt, dt, threadsEv, threadsFv
  use mod_constant, only : Cp, gamma_1, Cp_over_Pr, one_third, two_third
  use load_smem_visc2
  implicit none
  private
  public calc_Ev2, calc_Fv2
contains
  attributes(global) subroutine calc_Ev2(nx, ny, dx, dy, Q, T, mu, E)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    real(8), intent(in), device, contiguous    :: dx(nx-1)
    real(8), intent(in), device, contiguous    :: dy(ny-1)
    real(8), intent(in), device, contiguous    :: Q(nx,4,ny)
    real(8), intent(in), device, contiguous    :: T(nx,ny)
    real(8), intent(in), device, contiguous    :: mu(nx,ny)
    real(8), intent(inout), device, contiguous :: E(4,nx-1,ny-2)
    integer, parameter :: sx = threadsEv%x + 1
    integer, parameter :: sy = threadsEv%y
    real(8), dimension(0:sx*sy-1), shared :: u, v
    integer i, j, it, jt, idx
    real(8) viscous_work, txx, txy
    real(8) mux, mvx, muy, mvy 
    it = threadIdx%x
    jt = threadIdx%y
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    idx = (it-1) + (jt-1)*sx
    call load_smem_visc2_x(it, jt, j, idx, nx, ny, Q, u, v)
    if (nx-1 < i .or. ny-1 < j) return
    block
      real(8) mudx
      mudx         = 0.5d0 * (mu(i,j) + mu(i+1,j)) * dx(i)
      viscous_work = Cp_over_Pr * mudx * (-T(i,j) + T(i+1,j))
      mux          = mudx * (-u(idx) + u(idx+1))
      mvx          = mudx * (-v(idx) + v(idx+1))
    end block
    block
      real(8) my1, my2
      my1 = 0.0625d0 * (mu(i,j-1) + mu(i,j  ) + mu(i+1,j-1) + mu(i+1,j  )) ! my / 4
      my2 = 0.0625d0 * (mu(i,j  ) + mu(i,j+1) + mu(i+1,j  ) + mu(i+1,j+1)) ! my / 4
      muy = (my1 * (-Q(i,2,j-1) - Q(i+1,2,j-1)) + (my1 - my2) * (u(idx) + u(idx+1)) &
           + my2 * ( Q(i,2,j+1) + Q(i+1,2,j+1))) * dy(j)
      mvy = (my1 * (-Q(i,3,j-1) - Q(i+1,3,j-1)) + (my1 - my2) * (v(idx) + v(idx+1)) &
           + my2 * ( Q(i,3,j+1) + Q(i+1,3,j+1))) * dy(j)
    end block
    txx = two_third * (2.d0 * mux - mvy)
    txy = muy + mvx
    viscous_work = viscous_work + 0.5d0 * ((u(idx) + u(idx+1)) * txx + (v(idx) + v(idx+1)) * txy)
    E(2,i,j-1) = E(2,i,j-1) - txx
    E(3,i,j-1) = E(3,i,j-1) - txy
    E(4,i,j-1) = E(4,i,j-1) - viscous_work
  end subroutine calc_Ev2


  attributes(global) subroutine calc_Fv2(nx, ny, dy, dx, Q, T, mu, F)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    real(8), intent(in), device, contiguous    :: dy(ny-1)
    real(8), intent(in), device, contiguous    :: dx(nx-1)
    real(8), intent(in), device, contiguous    :: Q(nx,4,ny)
    real(8), intent(in), device, contiguous    :: T(nx,ny)
    real(8), intent(in), device, contiguous    :: mu(nx,ny)
    real(8), intent(inout), device, contiguous :: F(4,nx-2,ny-1)
    integer, parameter :: sx = threadsFv%x
    integer, parameter :: sy = threadsFv%y + 1
    real(8), dimension(0:sx*sy-1), shared :: u, v
    integer i, j, it, jt, idx
    real(8) viscous_work, tyx, tyy
    real(8) muy, mvy, mux, mvx
    it = threadIdx%x
    jt = threadIdx%y
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    idx = (jt-1) + (it-1)*sy
    call load_smem_visc2_y(it, jt, i, idx, nx, ny, Q, u, v)
    if (nx-1 < i .or. ny-1 < j) return
    block
      real(8) mudy
      mudy         = 0.5d0 * (mu(i,j) + mu(i,j+1)) * dy(j)
      viscous_work = Cp_over_Pr * mudy * (-T(i,j) + T(i,j+1))
      muy          = mudy * (-u(idx) + u(idx+1))
      mvy          = mudy * (-v(idx) + v(idx+1))
    end block
    block
      real(8) mx1, mx2
      mx1 = 0.0625d0 * (mu(i-1,j) + mu(i,  j) + mu(i-1,j+1) + mu(i,  j+1)) ! mx / 4
      mx2 = 0.0625d0 * (mu(i,  j) + mu(i+1,j) + mu(i,  j+1) + mu(i+1,j+1)) ! mx / 4
      mux = (mx1 * (-Q(i-1,2,j) - Q(i-1,2,j+1)) + (mx1 - mx2) * (u(idx) + u(idx+1)) &
           + mx2 * ( Q(i+1,2,j) + Q(i+1,2,j+1))) * dx(i)
      mvx = (mx1 * (-Q(i-1,3,j) - Q(i-1,3,j+1)) + (mx1 - mx2) * (v(idx) + v(idx+1)) &
           + mx2 * ( Q(i+1,3,j) + Q(i+1,3,j+1))) * dx(i)
    end block
    tyx = muy + mvx
    tyy = two_third * (2.d0 * mvy - mux)
    viscous_work = viscous_work + 0.5d0 * ((u(idx) + u(idx+1)) * tyx + (v(idx) + v(idx+1)) * tyy)
    F(2,i-1,j) = F(2,i-1,j) - tyx
    F(3,i-1,j) = F(3,i-1,j) - tyy
    F(4,i-1,j) = F(4,i-1,j) - viscous_work
  end subroutine calc_Fv2
end module calc_visc2

