module calc_visc2
  use mod_globals, only : gamma, R, Pr, Prt, dt, threadsEv, threadsFv, threadsGv
  use mod_constant, only : Cp, gamma_1, Cp_over_Pr, one_third, two_third
  use load_smem_visc2
  implicit none
  private
  public calc_Ev2, calc_Ev_LES2, calc_Fv2, calc_Fv_LES2, calc_Gv2, calc_Gv_LES2, &
         calc_Ev2_koff, calc_Fv2_koff, calc_Gv2_koff, &
         calc_Ev_LES2_koff, calc_Fv_LES2_koff, calc_Gv_LES2_koff
contains
  attributes(global) subroutine calc_Ev2(nx, ny, nz, dx, dy, dz, Q, T, mu, E)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), device, contiguous    :: dx(nx-1)
    real(8), intent(in), device, contiguous    :: dy(ny-1)
    real(8), intent(in), device, contiguous    :: dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2)
    integer, parameter :: sx = threadsEv%x + 1
    integer, parameter :: sy = threadsEv%y
    integer, parameter :: sz = threadsEv%z
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, txx, txy, txz
    real(8) mux, mvx, mwx, muy, mvy, muz, mwz
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    idx = (it-1) + (jt-1)*sx + (kt-1)*sx*sy
    call load_smem_visc2_x(it, jt, kt, j, k, idx, nx, ny, nz, Q, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    block
      real(8) mudx
      mudx         = 0.5d0 * (mu(i,j,k) + mu(i+1,j,k)) * dx(i)
      viscous_work = Cp_over_Pr * mudx * (-T(i,j,k) + T(i+1,j,k))
      mux          = mudx * (-u(idx) + u(idx+1))
      mvx          = mudx * (-v(idx) + v(idx+1))
      mwx          = mudx * (-w(idx) + w(idx+1))
    end block
    block
      real(8) my1, my2
      my1 = 0.0625d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i+1,j-1,k) + mu(i+1,j,  k)) ! my / 4
      my2 = 0.0625d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i+1,j,  k) + mu(i+1,j+1,k)) ! my / 4
      muy = (my1 * (-Q(i,2,j-1,k) - Q(i+1,2,j-1,k)) + (my1 - my2) * (u(idx) + u(idx+1)) &
           + my2 * ( Q(i,2,j+1,k) + Q(i+1,2,j+1,k))) * dy(j)
      mvy = (my1 * (-Q(i,3,j-1,k) - Q(i+1,3,j-1,k)) + (my1 - my2) * (v(idx) + v(idx+1)) &
           + my2 * ( Q(i,3,j+1,k) + Q(i+1,3,j+1,k))) * dy(j)
    end block
    block
      real(8) mz1, mz2
      mz1 = 0.0625d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i+1,j,k-1) + mu(i+1,j,k  )) ! mz / 4
      mz2 = 0.0625d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i+1,j,k  ) + mu(i+1,j,k+1)) ! mz / 4
      muz = (mz1 * (-Q(i,2,j,k-1) - Q(i+1,2,j,k-1)) + (mz1 - mz2) * (u(idx) + u(idx+1)) &
           + mz2 * ( Q(i,2,j,k+1) + Q(i+1,2,j,k+1))) * dz(k)
      mwz = (mz1 * (-Q(i,4,j,k-1) - Q(i+1,4,j,k-1)) + (mz1 - mz2) * (w(idx) + w(idx+1)) &
           + mz2 * ( Q(i,4,j,k+1) + Q(i+1,4,j,k+1))) * dz(k)
    end block
    txx = two_third * (2.d0 * mux - mvy - mwz)
    txy = muy + mvx
    txz = mwx + muz
    viscous_work = viscous_work + 0.5d0 * ((u(idx) + u(idx+1)) * txx + (v(idx) + v(idx+1)) * txy + (w(idx) + w(idx+1)) * txz)
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - viscous_work
  end subroutine calc_Ev2


  attributes(global) subroutine calc_Ev_LES2(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, E)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), device, contiguous    :: dx(nx-1)
    real(8), intent(in), device, contiguous    :: dy(ny-1)
    real(8), intent(in), device, contiguous    :: dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2)
    integer, parameter :: sx = threadsEv%x + 1
    integer, parameter :: sy = threadsEv%y
    integer, parameter :: sz = threadsEv%z
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, Hsgs, txx, txy, txz
    real(8) mux, muxsgs, mvx, mvxsgs, mwx, mwxsgs
    real(8) muy, muysgs, mvy, mvysgs, muz, muzsgs, mwz, mwzsgs
    real(8) mysgs1, mysgs2, mzsgs1, mzsgs2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    idx = (it-1) + (jt-1)*sx + (kt-1)*sx*sy
    call load_smem_visc2_x(it, jt, kt, j, k, idx, nx, ny, nz, Q, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    block
      real(8) mudx
      mudx         = 0.5d0 * (mu(i,j,k) + mu(i+1,j,k)) * dx(i)
      viscous_work = Cp_over_Pr * mudx * (-T(i,j,k) + T(i+1,j,k))
      mux          = mudx * (-u(idx) + u(idx+1))
      mvx          = mudx * (-v(idx) + v(idx+1))
      mwx          = mudx * (-w(idx) + w(idx+1))
    end block
    mysgs1 = 0.0625d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i+1,j-1,k) + mut(i+1,j,k))   ! mysgs / 4
    mysgs2 = 0.0625d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i+1,j,k) + mut(i+1,j+1,k))   ! mysgs / 4
    mzsgs1 = 0.0625d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i+1,j,k-1) + mut(i+1,j,k))   ! mzsgs / 4
    mzsgs2 = 0.0625d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i+1,j,k) + mut(i+1,j,k+1))   ! mzsgs / 4
    block
      real(8) mxsgsdx, H1, H2
      mxsgsdx = 0.5d0 * (mut(i,j,k) + mut(i+1,j,k)) * dx(i)
      H1 = Cp * T(i,j,k)   + 0.5d0 * (u(idx)**2   + v(idx)**2   + w(idx)**2)   + qc2(i,j,k)
      H2 = Cp * T(i+1,j,k) + 0.5d0 * (u(idx+1)**2 + v(idx+1)**2 + w(idx+1)**2) + qc2(i+1,j,k)
      Hsgs = -mxsgsdx * (-H1 + H2) / Prt
      muxsgs = mxsgsdx * (-u(idx) + u(idx+1))
      mvxsgs = mxsgsdx * (-v(idx) + v(idx+1))
      mwxsgs = mxsgsdx * (-w(idx) + w(idx+1))
    end block
    block
      real(8) my1, my2
      my1    = 0.0625d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i+1,j-1,k) + mu(i+1,j,  k)) ! my / 4
      my2    = 0.0625d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i+1,j,  k) + mu(i+1,j+1,k)) ! my / 4
      muy    = (my1    * (-Q(i,2,j-1,k) - Q(i+1,2,j-1,k)) + (my1 - my2) * (u(idx) + u(idx+1)) &
              + my2    * ( Q(i,2,j+1,k) + Q(i+1,2,j+1,k))) * dy(j)
      muysgs = (mysgs1 * (-Q(i,2,j-1,k) - Q(i+1,2,j-1,k)) + (mysgs1 - mysgs2) * (u(idx) + u(idx+1)) &
              + mysgs2 * ( Q(i,2,j+1,k) + Q(i+1,2,j+1,k))) * dy(j)
      mvy    = (my1    * (-Q(i,3,j-1,k) - Q(i+1,3,j-1,k)) + (my1 - my2) * (v(idx) + v(idx+1)) &
              + my2    * ( Q(i,3,j+1,k) + Q(i+1,3,j+1,k))) * dy(j)
      mvysgs = (mysgs1 * (-Q(i,3,j-1,k) - Q(i+1,3,j-1,k)) + (mysgs1 - mysgs2) * (v(idx) + v(idx+1)) &
              + mysgs2 * ( Q(i,3,j+1,k) + Q(i+1,3,j+1,k))) * dy(j)
    end block
    block
      real(8) mz1, mz2
      mz1    = 0.0625d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i+1,j,k-1) + mu(i+1,j,k  )) ! mz / 4
      mz2    = 0.0625d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i+1,j,k  ) + mu(i+1,j,k+1)) ! mz / 4
      muz    = (mz1    * (-Q(i,2,j,k-1) - Q(i+1,2,j,k-1)) + (mz1 - mz2) * (u(idx) + u(idx+1)) &
              + mz2    * ( Q(i,2,j,k+1) + Q(i+1,2,j,k+1))) * dz(k)
      muzsgs = (mzsgs1 * (-Q(i,2,j,k-1) - Q(i+1,2,j,k-1)) + (mzsgs1 - mzsgs2) * (u(idx) + u(idx+1)) &
              + mzsgs2 * ( Q(i,2,j,k+1) + Q(i+1,2,j,k+1))) * dz(k)
      mwz    = (mz1    * (-Q(i,4,j,k-1) - Q(i+1,4,j,k-1)) + (mz1 - mz2) * (w(idx) + w(idx+1)) &
              + mz2    * ( Q(i,4,j,k+1) + Q(i+1,4,j,k+1))) * dz(k)
      mwzsgs = (mzsgs1 * (-Q(i,4,j,k-1) - Q(i+1,4,j,k-1)) + (mzsgs1 - mzsgs2) * (w(idx) + w(idx+1)) &
              + mzsgs2 * ( Q(i,4,j,k+1) + Q(i+1,4,j,k+1))) * dz(k)
    end block
    txx = two_third * (2.d0 * mux - mvy - mwz)
    txy = muy + mvx
    txz = mwx + muz
    viscous_work = viscous_work + 0.5d0 * ((u(idx) + u(idx+1)) * txx + (v(idx) + v(idx+1)) * txy + (w(idx) + w(idx+1)) * txz)
    txx = txx + two_third * (2.d0 * muxsgs - mvysgs - mwzsgs)
    txy = txy + muysgs + mvxsgs
    txz = txz + mwxsgs + muzsgs
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (viscous_work + Hsgs)
  end subroutine calc_Ev_LES2


  attributes(global) subroutine calc_Fv2(nx, ny, nz, dy, dx, dz, Q, T, mu, F)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), device, contiguous    :: dy(ny-1)
    real(8), intent(in), device, contiguous    :: dx(nx-1)
    real(8), intent(in), device, contiguous    :: dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: F(5,nx-2,ny-1,nz-2)
    integer, parameter :: sx = threadsFv%x
    integer, parameter :: sy = threadsFv%y + 1
    integer, parameter :: sz = threadsFv%z
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, tyx, tyy, tyz
    real(8) muy, mvy, mwy, mvz, mwz, mux, mvx
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    idx = (jt-1) + (it-1)*sy + (kt-1)*sy*sx
    call load_smem_visc2_y(it, jt, kt, i, k, idx, nx, ny, nz, Q, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    block
      real(8) mudy
      mudy         = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k)) * dy(j)
      viscous_work = Cp_over_Pr * mudy * (-T(i,j,k) + T(i,j+1,k))
      muy          = mudy * (-u(idx) + u(idx+1))
      mvy          = mudy * (-v(idx) + v(idx+1))
      mwy          = mudy * (-w(idx) + w(idx+1))
    end block
    block
      real(8) mx1, mx2
      mx1 = 0.0625d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j+1,k) + mu(i,  j+1,k)) ! mx / 4
      mx2 = 0.0625d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j+1,k) + mu(i+1,j+1,k)) ! mx / 4
      mux = (mx1 * (-Q(i-1,2,j,k) - Q(i-1,2,j+1,k)) + (mx1 - mx2) * (u(idx) + u(idx+1)) &
           + mx2 * ( Q(i+1,2,j,k) + Q(i+1,2,j+1,k))) * dx(i)
      mvx = (mx1 * (-Q(i-1,3,j,k) - Q(i-1,3,j+1,k)) + (mx1 - mx2) * (v(idx) + v(idx+1)) &
           + mx2 * ( Q(i+1,3,j,k) + Q(i+1,3,j+1,k))) * dx(i)
    end block
    block
      real(8) mz1, mz2
      mz1 = 0.0625d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i,j+1,k-1) + mu(i,j+1,k  )) ! mz / 4
      mz2 = 0.0625d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i,j+1,k  ) + mu(i,j+1,k+1)) ! mz / 4
      mvz = (mz1 * (-Q(i,3,j,k-1) - Q(i,3,j+1,k-1)) + (mz1 - mz2) * (v(idx) + v(idx+1)) &
           + mz2 * ( Q(i,3,j,k+1) + Q(i,3,j+1,k+1))) * dz(k)
      mwz = (mz1 * (-Q(i,4,j,k-1) - Q(i,4,j+1,k-1)) + (mz1 - mz2) * (w(idx) + w(idx+1)) &
           + mz2 * ( Q(i,4,j,k+1) + Q(i,4,j+1,k+1))) * dz(k)
    end block
    tyx = muy + mvx
    tyy = two_third * (2.d0 * mvy - mwz - mux)
    tyz = mvz + mwy
    viscous_work = viscous_work + 0.5d0 * ((u(idx) + u(idx+1)) * tyx + (v(idx) + v(idx+1)) * tyy + (w(idx) + w(idx+1)) * tyz)
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - viscous_work
  end subroutine calc_Fv2


  attributes(global) subroutine calc_Fv_LES2(nx, ny, nz, dy, dx, dz, Q, T, mu, mut, qc2, F)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), device, contiguous    :: dy(ny-1)
    real(8), intent(in), device, contiguous    :: dx(nx-1)
    real(8), intent(in), device, contiguous    :: dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: F(5,nx-2,ny-1,nz-2)
    integer, parameter :: sx = threadsFv%x
    integer, parameter :: sy = threadsFv%y + 1
    integer, parameter :: sz = threadsFv%z
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, Hsgs, tyx, tyy, tyz
    real(8) muy, muysgs, mvy, mvysgs, mwy, mwysgs
    real(8) mvz, mvzsgs, mwz, mwzsgs, mux, muxsgs, mvx, mvxsgs
    real(8) mzsgs1, mzsgs2, mxsgs1, mxsgs2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    idx = (jt-1) + (it-1)*sy + (kt-1)*sy*sx
    call load_smem_visc2_y(it, jt, kt, i, k, idx, nx, ny, nz, Q, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    block
      real(8) mudy
      mudy         = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k)) * dy(j)
      viscous_work = Cp_over_Pr * mudy * (-T(i,j,k) + T(i,j+1,k))
      muy          = mudy * (-u(idx) + u(idx+1))
      mvy          = mudy * (-v(idx) + v(idx+1))
      mwy          = mudy * (-w(idx) + w(idx+1))
    end block
    mzsgs1 = 0.0625d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i,j+1,k-1) + mut(i,j+1,k))     ! mzsgs / 4
    mzsgs2 = 0.0625d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i,j+1,k) + mut(i,j+1,k+1))     ! mzsgs / 4
    mxsgs1 = 0.0625d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j+1,k) + mut(i,j+1,k))     ! mxsgs / 4
    mxsgs2 = 0.0625d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j+1,k) + mut(i+1,j+1,k))     ! mxsgs / 4
    block
      real(8) mysgsdy, H1, H2
      mysgsdy  = 0.5d0 * (mut(i,j,k) + mut(i,j+1,k)) * dy(j)
      H1 = Cp * T(i,j,k)   + 0.5d0 * (u(idx)**2   + v(idx)**2   + w(idx)**2)   + qc2(i,j,k)
      H2 = Cp * T(i,j+1,k) + 0.5d0 * (u(idx+1)**2 + v(idx+1)**2 + w(idx+1)**2) + qc2(i,j+1,k)
      Hsgs = -mysgsdy * (-H1 + H2) / Prt
      muysgs = mysgsdy * (-u(idx) + u(idx+1))
      mvysgs = mysgsdy * (-v(idx) + v(idx+1))
      mwysgs = mysgsdy * (-w(idx) + w(idx+1))
    end block
    block
      real(8) mx1, mx2
      mx1    = 0.0625d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j+1,k) + mu(i,  j+1,k)) ! mx / 4
      mx2    = 0.0625d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j+1,k) + mu(i+1,j+1,k)) ! mx / 4
      mux    = (mx1    * (-Q(i-1,2,j,k) - Q(i-1,2,j+1,k)) + (mx1 - mx2) * (u(idx) + u(idx+1)) &
              + mx2    * ( Q(i+1,2,j,k) + Q(i+1,2,j+1,k))) * dx(i)
      muxsgs = (mxsgs1 * (-Q(i-1,2,j,k) - Q(i-1,2,j+1,k)) + (mxsgs1 - mxsgs2) * (u(idx) + u(idx+1)) &
              + mxsgs2 * ( Q(i+1,2,j,k) + Q(i+1,2,j+1,k))) * dx(i)
      mvx    = (mx1    * (-Q(i-1,3,j,k) - Q(i-1,3,j+1,k)) + (mx1 - mx2) * (v(idx) + v(idx+1)) &
              + mx2    * ( Q(i+1,3,j,k) + Q(i+1,3,j+1,k))) * dx(i)
      mvxsgs = (mxsgs1 * (-Q(i-1,3,j,k) - Q(i-1,3,j+1,k)) + (mxsgs1 - mxsgs2) * (v(idx) + v(idx+1)) &
              + mxsgs2 * ( Q(i+1,3,j,k) + Q(i+1,3,j+1,k))) * dx(i)
    end block
    block
      real(8) mz1, mz2
      mz1    = 0.0625d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i,j+1,k-1) + mu(i,j+1,k  )) ! mz / 4
      mz2    = 0.0625d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i,j+1,k  ) + mu(i,j+1,k+1)) ! mz / 4
      mvz    = (mz1    * (-Q(i,3,j,k-1) - Q(i,3,j+1,k-1)) + (mz1 - mz2) * (v(idx) + v(idx+1)) &
              + mz2    * ( Q(i,3,j,k+1) + Q(i,3,j+1,k+1))) * dz(k)
      mvzsgs = (mzsgs1 * (-Q(i,3,j,k-1) - Q(i,3,j+1,k-1)) + (mzsgs1 - mzsgs2) * (v(idx) + v(idx+1)) &
              + mzsgs2 * ( Q(i,3,j,k+1) + Q(i,3,j+1,k+1))) * dz(k)
      mwz    = (mz1    * (-Q(i,4,j,k-1) - Q(i,4,j+1,k-1)) + (mz1 - mz2) * (w(idx) + w(idx+1)) &
              + mz2    * ( Q(i,4,j,k+1) + Q(i,4,j+1,k+1))) * dz(k)
      mwzsgs = (mzsgs1 * (-Q(i,4,j,k-1) - Q(i,4,j+1,k-1)) + (mzsgs1 - mzsgs2) * (w(idx) + w(idx+1)) &
              + mzsgs2 * ( Q(i,4,j,k+1) + Q(i,4,j+1,k+1))) * dz(k)
    end block
    tyx = muy + mvx
    tyy = two_third * (2.d0 * mvy - mwz - mux)
    tyz = mvz + mwy
    viscous_work = viscous_work + 0.5d0 * ((u(idx) + u(idx+1)) * tyx + (v(idx) + v(idx+1)) * tyy + (w(idx) + w(idx+1)) * tyz)
    tyx = tyx + muysgs + mvxsgs
    tyy = tyy + two_third * (2.d0 * mvysgs - mwzsgs - muxsgs)
    tyz = tyz + mvzsgs + mwysgs
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (viscous_work + Hsgs)
  end subroutine calc_Fv_LES2


  attributes(global) subroutine calc_Gv2(nx, ny, nz, dx, dy, dz, Q, T, mu, G)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), device, contiguous    :: dx(nx-1)
    real(8), intent(in), device, contiguous    :: dy(ny-1)
    real(8), intent(in), device, contiguous    :: dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    integer, parameter :: sx = threadsGv%x
    integer, parameter :: sy = threadsGv%y
    integer, parameter :: sz = threadsGv%z + 1
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, tzx, tzy, tzz
    real(8) muz, mvz, mwz, mwx, mux, mvy, mwy
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt
    idx = (kt-1) + (jt-1)*sz + (it-1)*sz*sy
    call load_smem_visc2_z(it, jt, kt, i, j, idx, nx, ny, nz, Q, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    block
      real(8) mudz
      mudz         = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1)) * dz(k)
      viscous_work = Cp_over_Pr * mudz * (-T(i,j,k) + T(i,j,k+1))
      muz          = mudz * (-u(idx) + u(idx+1))
      mvz          = mudz * (-v(idx) + v(idx+1))
      mwz          = mudz * (-w(idx) + w(idx+1))
    end block
    block
      real(8) mx1, mx2
      mx1 = 0.0625d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j,k+1) + mu(i,  j,k+1)) ! mx / 4
      mx2 = 0.0625d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j,k+1) + mu(i+1,j,k+1)) ! mx  4
      mux = (mx1 * (-Q(i-1,2,j,k) - Q(i-1,2,j,k+1)) + (mx1 - mx2) * (u(idx) + u(idx+1)) &
           + mx2 * ( Q(i+1,2,j,k) + Q(i+1,2,j,k+1))) * dx(i)
      mwx = (mx1 * (-Q(i-1,4,j,k) - Q(i-1,4,j,k+1)) + (mx1 - mx2) * (w(idx) + w(idx+1)) &
           + mx2 * ( Q(i+1,4,j,k) + Q(i+1,4,j,k+1))) * dx(i)
    end block
    block
      real(8) my1, my2
      my1 = 0.0625d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i,j-1,k+1) + mu(i,j,  k+1)) ! my / 4
      my2 = 0.0625d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i,j,  k+1) + mu(i,j+1,k+1)) ! my / 4
      mvy = (my1 * (-Q(i,3,j-1,k) - Q(i,3,j-1,k+1)) + (my1 - my2) * (v(idx) + v(idx+1)) &
           + my2 * ( Q(i,3,j+1,k) + Q(i,3,j+1,k+1))) * dy(j)
      mwy = (my1 * (-Q(i,4,j-1,k) - Q(i,4,j-1,k+1)) + (my1 - my2) * (w(idx) + w(idx+1)) &
           + my2 * ( Q(i,4,j+1,k) + Q(i,4,j+1,k+1))) * dy(j)
    end block
    tzx  = mwx + muz
    tzy  = mvz + mwy
    tzz  = two_third * (2.d0 * mwz - mux - mvy)
    viscous_work = viscous_work + 0.5d0 * ((u(idx) + u(idx+1)) * tzx + (v(idx) + v(idx+1)) * tzy + (w(idx) + w(idx+1)) * tzz)
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - viscous_work
  end subroutine calc_Gv2


  attributes(global) subroutine calc_Gv_LES2(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, G)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), device, contiguous    :: dx(nx-1)
    real(8), intent(in), device, contiguous    :: dy(ny-1)
    real(8), intent(in), device, contiguous    :: dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    integer, parameter :: sx = threadsGv%x
    integer, parameter :: sy = threadsGv%y
    integer, parameter :: sz = threadsGv%z + 1
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, Hsgs, tzx, tzy, tzz
    real(8) muz, muzsgs, mvz, mvzsgs, mwz, mwzsgs
    real(8) mwx, mwxsgs, mux, muxsgs, mvy, mvysgs, mwy, mwysgs
    real(8) mxsgs1, mxsgs2, mysgs1, mysgs2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt
    idx = (kt-1) + (jt-1)*sz + (it-1)*sz*sy
    call load_smem_visc2_z(it, jt, kt, i, j, idx, nx, ny, nz, Q, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    block
      real(8) mudz
      mudz         = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1)) * dz(k)
      viscous_work = Cp_over_Pr * mudz * (-T(i,j,k) + T(i,j,k+1))
      muz          = mudz * (-u(idx) + u(idx+1))
      mvz          = mudz * (-v(idx) + v(idx+1))
      mwz          = mudz * (-w(idx) + w(idx+1))
    end block
    mxsgs1 = 0.0625d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j,k+1) + mut(i,j,k+1))     ! mxsgs / 4
    mxsgs2 = 0.0625d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j,k+1) + mut(i+1,j,k+1))     ! mxsgs / 4
    mysgs1 = 0.0625d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i,j-1,k+1) + mut(i,j,k+1))     ! mysgs / 4
    mysgs2 = 0.0625d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i,j,k+1) + mut(i,j+1,k+1))     ! mysgs / 4
    block
      real(8) mzsgsdz, H1, H2
      mzsgsdz  = 0.5d0 * (mut(i,j,k) + mut(i,j,k+1)) * dz(k)
      H1 = Cp * T(i,j,k)   + 0.5d0 * (u(idx)**2   + v(idx)**2   + w(idx)**2)   + qc2(i,j,k)
      H2 = Cp * T(i,j,k+1) + 0.5d0 * (u(idx+1)**2 + v(idx+1)**2 + w(idx+1)**2) + qc2(i,j,k+1)
      Hsgs = -mzsgsdz * (-H1 + H2) / Prt
      muzsgs = mzsgsdz * (-u(idx) + u(idx+1))
      mvzsgs = mzsgsdz * (-v(idx) + v(idx+1))
      mwzsgs = mzsgsdz * (-w(idx) + w(idx+1))
    end block
    block
      real(8) mx1, mx2
      mx1    = 0.0625d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j,k+1) + mu(i,  j,k+1)) ! mx / 4
      mx2    = 0.0625d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j,k+1) + mu(i+1,j,k+1)) ! mx / 4
      mux    = (mx1    * (-Q(i-1,2,j,k) - Q(i-1,2,j,k+1)) + (mx1 - mx2) * (u(idx) + u(idx+1)) &
              + mx2    * ( Q(i+1,2,j,k) + Q(i+1,2,j,k+1))) * dx(i)
      muxsgs = (mxsgs1 * (-Q(i-1,2,j,k) - Q(i-1,2,j,k+1)) + (mxsgs1 - mxsgs2) * (u(idx) + u(idx+1)) &
              + mxsgs2 * ( Q(i+1,2,j,k) + Q(i+1,2,j,k+1))) * dx(i)
      mwx    = (mx1    * (-Q(i-1,4,j,k) - Q(i-1,4,j,k+1)) + (mx1 - mx2) * (w(idx) + w(idx+1)) &
              + mx2    * ( Q(i+1,4,j,k) + Q(i+1,4,j,k+1))) * dx(i)
      mwxsgs = (mxsgs1 * (-Q(i-1,4,j,k) - Q(i-1,4,j,k+1)) + (mxsgs1 - mxsgs2) * (w(idx) + w(idx+1)) &
              + mxsgs2 * ( Q(i+1,4,j,k) + Q(i+1,4,j,k+1))) * dx(i)
    end block
    block
      real(8) my1, my2
      my1    = 0.0625d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i,j-1,k+1) + mu(i,j,  k+1)) ! my / 4
      my2    = 0.0625d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i,j,  k+1) + mu(i,j+1,k+1)) ! my / 4
      mvy    = (my1    * (-Q(i,3,j-1,k) - Q(i,3,j-1,k+1)) + (my1 - my2) * (v(idx) + v(idx+1)) &
              + my2    * ( Q(i,3,j+1,k) + Q(i,3,j+1,k+1))) * dy(j)
      mvysgs = (mysgs1 * (-Q(i,3,j-1,k) - Q(i,3,j-1,k+1)) + (mysgs1 - mysgs2) * (v(idx) + v(idx+1)) &
              + mysgs2 * ( Q(i,3,j+1,k) + Q(i,3,j+1,k+1))) * dy(j)
      mwy    = (my1    * (-Q(i,4,j-1,k) - Q(i,4,j-1,k+1)) + (my1 - my2) * (w(idx) + w(idx+1)) &
              + my2    * ( Q(i,4,j+1,k) + Q(i,4,j+1,k+1))) * dy(j)
      mwysgs = (mysgs1 * (-Q(i,4,j-1,k) - Q(i,4,j-1,k+1)) + (mysgs1 - mysgs2) * (w(idx) + w(idx+1)) &
              + mysgs2 * ( Q(i,4,j+1,k) + Q(i,4,j+1,k+1))) * dy(j)
    end block
    tzx = mwx + muz
    tzy = mvz + mwy
    tzz = two_third * (2.d0 * mwz - mux - mvy)
    viscous_work = viscous_work + 0.5d0 * ((u(idx) + u(idx+1)) * tzx + (v(idx) + v(idx+1)) * tzy + (w(idx) + w(idx+1)) * tzz)
    tzx = tzx + mwxsgs + muzsgs
    tzy = tzy + mvzsgs + mwysgs
    tzz = tzz + two_third * (2.d0 * mwzsgs - muxsgs - mvysgs)
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (viscous_work + Hsgs)
  end subroutine calc_Gv_LES2


  !> Like calc_Ev2 but restricted to k in [k_lo, k_hi]; k_lo absorbs the +1 offset.
  attributes(global) subroutine calc_Ev2_koff(nx, ny, nz, dx, dy, dz, Q, T, mu, E, k_lo, k_hi)
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: dx(nx-1), dy(ny-1), dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2)
    integer, parameter :: sx = threadsEv%x + 1
    integer, parameter :: sy = threadsEv%y
    integer, parameter :: sz = threadsEv%z
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, txx, txy, txz
    real(8) mux, mvx, mwx, muy, mvy, muz, mwz
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt
    idx = (it-1) + (jt-1)*sx + (kt-1)*sx*sy
    call load_smem_visc2_x(it, jt, kt, j, k, idx, nx, ny, nz, Q, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. k_hi < k) return
    block
      real(8) mudx
      mudx         = 0.5d0 * (mu(i,j,k) + mu(i+1,j,k)) * dx(i)
      viscous_work = Cp_over_Pr * mudx * (-T(i,j,k) + T(i+1,j,k))
      mux          = mudx * (-u(idx) + u(idx+1))
      mvx          = mudx * (-v(idx) + v(idx+1))
      mwx          = mudx * (-w(idx) + w(idx+1))
    end block
    block
      real(8) my1, my2
      my1 = 0.0625d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i+1,j-1,k) + mu(i+1,j,  k))
      my2 = 0.0625d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i+1,j,  k) + mu(i+1,j+1,k))
      muy = (my1 * (-Q(i,2,j-1,k) - Q(i+1,2,j-1,k)) + (my1 - my2) * (u(idx) + u(idx+1)) &
           + my2 * ( Q(i,2,j+1,k) + Q(i+1,2,j+1,k))) * dy(j)
      mvy = (my1 * (-Q(i,3,j-1,k) - Q(i+1,3,j-1,k)) + (my1 - my2) * (v(idx) + v(idx+1)) &
           + my2 * ( Q(i,3,j+1,k) + Q(i+1,3,j+1,k))) * dy(j)
    end block
    block
      real(8) mz1, mz2
      mz1 = 0.0625d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i+1,j,k-1) + mu(i+1,j,k  ))
      mz2 = 0.0625d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i+1,j,k  ) + mu(i+1,j,k+1))
      muz = (mz1 * (-Q(i,2,j,k-1) - Q(i+1,2,j,k-1)) + (mz1 - mz2) * (u(idx) + u(idx+1)) &
           + mz2 * ( Q(i,2,j,k+1) + Q(i+1,2,j,k+1))) * dz(k)
      mwz = (mz1 * (-Q(i,4,j,k-1) - Q(i+1,4,j,k-1)) + (mz1 - mz2) * (w(idx) + w(idx+1)) &
           + mz2 * ( Q(i,4,j,k+1) + Q(i+1,4,j,k+1))) * dz(k)
    end block
    txx = two_third * (2.d0 * mux - mvy - mwz)
    txy = muy + mvx
    txz = mwx + muz
    viscous_work = viscous_work + 0.5d0 * ((u(idx) + u(idx+1)) * txx + (v(idx) + v(idx+1)) * txy + (w(idx) + w(idx+1)) * txz)
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - viscous_work
  end subroutine calc_Ev2_koff


  !> Like calc_Fv2 but restricted to k in [k_lo, k_hi]; k_lo absorbs the +1 offset.
  attributes(global) subroutine calc_Fv2_koff(nx, ny, nz, dy, dx, dz, Q, T, mu, F, k_lo, k_hi)
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: dy(ny-1), dx(nx-1), dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: F(5,nx-2,ny-1,nz-2)
    integer, parameter :: sx = threadsFv%x
    integer, parameter :: sy = threadsFv%y + 1
    integer, parameter :: sz = threadsFv%z
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, tyx, tyy, tyz
    real(8) muy, mvy, mwy, mvz, mwz, mux, mvx
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt
    idx = (jt-1) + (it-1)*sy + (kt-1)*sy*sx
    call load_smem_visc2_y(it, jt, kt, i, k, idx, nx, ny, nz, Q, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. k_hi < k) return
    block
      real(8) mudy
      mudy         = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k)) * dy(j)
      viscous_work = Cp_over_Pr * mudy * (-T(i,j,k) + T(i,j+1,k))
      muy          = mudy * (-u(idx) + u(idx+1))
      mvy          = mudy * (-v(idx) + v(idx+1))
      mwy          = mudy * (-w(idx) + w(idx+1))
    end block
    block
      real(8) mx1, mx2
      mx1 = 0.0625d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j+1,k) + mu(i,  j+1,k))
      mx2 = 0.0625d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j+1,k) + mu(i+1,j+1,k))
      mux = (mx1 * (-Q(i-1,2,j,k) - Q(i-1,2,j+1,k)) + (mx1 - mx2) * (u(idx) + u(idx+1)) &
           + mx2 * ( Q(i+1,2,j,k) + Q(i+1,2,j+1,k))) * dx(i)
      mvx = (mx1 * (-Q(i-1,3,j,k) - Q(i-1,3,j+1,k)) + (mx1 - mx2) * (v(idx) + v(idx+1)) &
           + mx2 * ( Q(i+1,3,j,k) + Q(i+1,3,j+1,k))) * dx(i)
    end block
    block
      real(8) mz1, mz2
      mz1 = 0.0625d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i,j+1,k-1) + mu(i,j+1,k  ))
      mz2 = 0.0625d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i,j+1,k  ) + mu(i,j+1,k+1))
      mvz = (mz1 * (-Q(i,3,j,k-1) - Q(i,3,j+1,k-1)) + (mz1 - mz2) * (v(idx) + v(idx+1)) &
           + mz2 * ( Q(i,3,j,k+1) + Q(i,3,j+1,k+1))) * dz(k)
      mwz = (mz1 * (-Q(i,4,j,k-1) - Q(i,4,j+1,k-1)) + (mz1 - mz2) * (w(idx) + w(idx+1)) &
           + mz2 * ( Q(i,4,j,k+1) + Q(i,4,j+1,k+1))) * dz(k)
    end block
    tyx = muy + mvx
    tyy = two_third * (2.d0 * mvy - mwz - mux)
    tyz = mvz + mwy
    viscous_work = viscous_work + 0.5d0 * ((u(idx) + u(idx+1)) * tyx + (v(idx) + v(idx+1)) * tyy + (w(idx) + w(idx+1)) * tyz)
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - viscous_work
  end subroutine calc_Fv2_koff


  !> Like calc_Gv2 but restricted to k in [k_lo, k_hi]; uses load_smem_visc2_z_koff.
  attributes(global) subroutine calc_Gv2_koff(nx, ny, nz, dx, dy, dz, Q, T, mu, G, k_lo, k_hi)
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: dx(nx-1), dy(ny-1), dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    integer, parameter :: sx = threadsGv%x
    integer, parameter :: sy = threadsGv%y
    integer, parameter :: sz = threadsGv%z + 1
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, tzx, tzy, tzz
    real(8) muz, mvz, mwz, mwx, mux, mvy, mwy
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt
    idx = (kt-1) + (jt-1)*sz + (it-1)*sz*sy
    call load_smem_visc2_z_koff(it, jt, kt, i, j, idx, nx, ny, nz, Q, u, v, w, k_lo)
    if (nx-1 < i .or. ny-1 < j .or. k_hi < k) return
    block
      real(8) mudz
      mudz         = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1)) * dz(k)
      viscous_work = Cp_over_Pr * mudz * (-T(i,j,k) + T(i,j,k+1))
      muz          = mudz * (-u(idx) + u(idx+1))
      mvz          = mudz * (-v(idx) + v(idx+1))
      mwz          = mudz * (-w(idx) + w(idx+1))
    end block
    block
      real(8) mx1, mx2
      mx1 = 0.0625d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j,k+1) + mu(i,  j,k+1))
      mx2 = 0.0625d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j,k+1) + mu(i+1,j,k+1))
      mux = (mx1 * (-Q(i-1,2,j,k) - Q(i-1,2,j,k+1)) + (mx1 - mx2) * (u(idx) + u(idx+1)) &
           + mx2 * ( Q(i+1,2,j,k) + Q(i+1,2,j,k+1))) * dx(i)
      mwx = (mx1 * (-Q(i-1,4,j,k) - Q(i-1,4,j,k+1)) + (mx1 - mx2) * (w(idx) + w(idx+1)) &
           + mx2 * ( Q(i+1,4,j,k) + Q(i+1,4,j,k+1))) * dx(i)
    end block
    block
      real(8) my1, my2
      my1 = 0.0625d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i,j-1,k+1) + mu(i,j,  k+1))
      my2 = 0.0625d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i,j,  k+1) + mu(i,j+1,k+1))
      mvy = (my1 * (-Q(i,3,j-1,k) - Q(i,3,j-1,k+1)) + (my1 - my2) * (v(idx) + v(idx+1)) &
           + my2 * ( Q(i,3,j+1,k) + Q(i,3,j+1,k+1))) * dy(j)
      mwy = (my1 * (-Q(i,4,j-1,k) - Q(i,4,j-1,k+1)) + (my1 - my2) * (w(idx) + w(idx+1)) &
           + my2 * ( Q(i,4,j+1,k) + Q(i,4,j+1,k+1))) * dy(j)
    end block
    tzx  = mwx + muz
    tzy  = mvz + mwy
    tzz  = two_third * (2.d0 * mwz - mux - mvy)
    viscous_work = viscous_work + 0.5d0 * ((u(idx) + u(idx+1)) * tzx + (v(idx) + v(idx+1)) * tzy + (w(idx) + w(idx+1)) * tzz)
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - viscous_work
  end subroutine calc_Gv2_koff


  !> Like calc_Ev_LES2 but restricted to k in [k_lo, k_hi]; k_lo absorbs the +1 offset.
  attributes(global) subroutine calc_Ev_LES2_koff(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, E, k_lo, k_hi)
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: dx(nx-1), dy(ny-1), dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: E(5,nx-1,ny-2,nz-2)
    integer, parameter :: sx = threadsEv%x + 1
    integer, parameter :: sy = threadsEv%y
    integer, parameter :: sz = threadsEv%z
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, Hsgs, txx, txy, txz
    real(8) mux, muxsgs, mvx, mvxsgs, mwx, mwxsgs
    real(8) muy, muysgs, mvy, mvysgs, muz, muzsgs, mwz, mwzsgs
    real(8) mysgs1, mysgs2, mzsgs1, mzsgs2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt
    idx = (it-1) + (jt-1)*sx + (kt-1)*sx*sy
    call load_smem_visc2_x(it, jt, kt, j, k, idx, nx, ny, nz, Q, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. k_hi < k) return
    block
      real(8) mudx
      mudx         = 0.5d0 * (mu(i,j,k) + mu(i+1,j,k)) * dx(i)
      viscous_work = Cp_over_Pr * mudx * (-T(i,j,k) + T(i+1,j,k))
      mux          = mudx * (-u(idx) + u(idx+1))
      mvx          = mudx * (-v(idx) + v(idx+1))
      mwx          = mudx * (-w(idx) + w(idx+1))
    end block
    mysgs1 = 0.0625d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i+1,j-1,k) + mut(i+1,j,k))
    mysgs2 = 0.0625d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i+1,j,k) + mut(i+1,j+1,k))
    mzsgs1 = 0.0625d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i+1,j,k-1) + mut(i+1,j,k))
    mzsgs2 = 0.0625d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i+1,j,k) + mut(i+1,j,k+1))
    block
      real(8) mxsgsdx, H1, H2
      mxsgsdx = 0.5d0 * (mut(i,j,k) + mut(i+1,j,k)) * dx(i)
      H1 = Cp * T(i,j,k)   + 0.5d0 * (u(idx)**2   + v(idx)**2   + w(idx)**2)   + qc2(i,j,k)
      H2 = Cp * T(i+1,j,k) + 0.5d0 * (u(idx+1)**2 + v(idx+1)**2 + w(idx+1)**2) + qc2(i+1,j,k)
      Hsgs = -mxsgsdx * (-H1 + H2) / Prt
      muxsgs = mxsgsdx * (-u(idx) + u(idx+1))
      mvxsgs = mxsgsdx * (-v(idx) + v(idx+1))
      mwxsgs = mxsgsdx * (-w(idx) + w(idx+1))
    end block
    block
      real(8) my1, my2
      my1    = 0.0625d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i+1,j-1,k) + mu(i+1,j,  k))
      my2    = 0.0625d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i+1,j,  k) + mu(i+1,j+1,k))
      muy    = (my1    * (-Q(i,2,j-1,k) - Q(i+1,2,j-1,k)) + (my1 - my2) * (u(idx) + u(idx+1)) &
              + my2    * ( Q(i,2,j+1,k) + Q(i+1,2,j+1,k))) * dy(j)
      muysgs = (mysgs1 * (-Q(i,2,j-1,k) - Q(i+1,2,j-1,k)) + (mysgs1 - mysgs2) * (u(idx) + u(idx+1)) &
              + mysgs2 * ( Q(i,2,j+1,k) + Q(i+1,2,j+1,k))) * dy(j)
      mvy    = (my1    * (-Q(i,3,j-1,k) - Q(i+1,3,j-1,k)) + (my1 - my2) * (v(idx) + v(idx+1)) &
              + my2    * ( Q(i,3,j+1,k) + Q(i+1,3,j+1,k))) * dy(j)
      mvysgs = (mysgs1 * (-Q(i,3,j-1,k) - Q(i+1,3,j-1,k)) + (mysgs1 - mysgs2) * (v(idx) + v(idx+1)) &
              + mysgs2 * ( Q(i,3,j+1,k) + Q(i+1,3,j+1,k))) * dy(j)
    end block
    block
      real(8) mz1, mz2
      mz1    = 0.0625d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i+1,j,k-1) + mu(i+1,j,k  ))
      mz2    = 0.0625d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i+1,j,k  ) + mu(i+1,j,k+1))
      muz    = (mz1    * (-Q(i,2,j,k-1) - Q(i+1,2,j,k-1)) + (mz1 - mz2) * (u(idx) + u(idx+1)) &
              + mz2    * ( Q(i,2,j,k+1) + Q(i+1,2,j,k+1))) * dz(k)
      muzsgs = (mzsgs1 * (-Q(i,2,j,k-1) - Q(i+1,2,j,k-1)) + (mzsgs1 - mzsgs2) * (u(idx) + u(idx+1)) &
              + mzsgs2 * ( Q(i,2,j,k+1) + Q(i+1,2,j,k+1))) * dz(k)
      mwz    = (mz1    * (-Q(i,4,j,k-1) - Q(i+1,4,j,k-1)) + (mz1 - mz2) * (w(idx) + w(idx+1)) &
              + mz2    * ( Q(i,4,j,k+1) + Q(i+1,4,j,k+1))) * dz(k)
      mwzsgs = (mzsgs1 * (-Q(i,4,j,k-1) - Q(i+1,4,j,k-1)) + (mzsgs1 - mzsgs2) * (w(idx) + w(idx+1)) &
              + mzsgs2 * ( Q(i,4,j,k+1) + Q(i+1,4,j,k+1))) * dz(k)
    end block
    txx = two_third * (2.d0 * mux - mvy - mwz)
    txy = muy + mvx
    txz = mwx + muz
    viscous_work = viscous_work + 0.5d0 * ((u(idx) + u(idx+1)) * txx + (v(idx) + v(idx+1)) * txy + (w(idx) + w(idx+1)) * txz)
    txx = txx + two_third * (2.d0 * muxsgs - mvysgs - mwzsgs)
    txy = txy + muysgs + mvxsgs
    txz = txz + mwxsgs + muzsgs
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (viscous_work + Hsgs)
  end subroutine calc_Ev_LES2_koff


  !> Like calc_Fv_LES2 but restricted to k in [k_lo, k_hi]; k_lo absorbs the +1 offset.
  attributes(global) subroutine calc_Fv_LES2_koff(nx, ny, nz, dy, dx, dz, Q, T, mu, mut, qc2, F, k_lo, k_hi)
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: dy(ny-1), dx(nx-1), dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: F(5,nx-2,ny-1,nz-2)
    integer, parameter :: sx = threadsFv%x
    integer, parameter :: sy = threadsFv%y + 1
    integer, parameter :: sz = threadsFv%z
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, Hsgs, tyx, tyy, tyz
    real(8) muy, muysgs, mvy, mvysgs, mwy, mwysgs
    real(8) mvz, mvzsgs, mwz, mwzsgs, mux, muxsgs, mvx, mvxsgs
    real(8) mzsgs1, mzsgs2, mxsgs1, mxsgs2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt
    idx = (jt-1) + (it-1)*sy + (kt-1)*sy*sx
    call load_smem_visc2_y(it, jt, kt, i, k, idx, nx, ny, nz, Q, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. k_hi < k) return
    block
      real(8) mudy
      mudy         = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k)) * dy(j)
      viscous_work = Cp_over_Pr * mudy * (-T(i,j,k) + T(i,j+1,k))
      muy          = mudy * (-u(idx) + u(idx+1))
      mvy          = mudy * (-v(idx) + v(idx+1))
      mwy          = mudy * (-w(idx) + w(idx+1))
    end block
    mzsgs1 = 0.0625d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i,j+1,k-1) + mut(i,j+1,k))
    mzsgs2 = 0.0625d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i,j+1,k) + mut(i,j+1,k+1))
    mxsgs1 = 0.0625d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j+1,k) + mut(i,j+1,k))
    mxsgs2 = 0.0625d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j+1,k) + mut(i+1,j+1,k))
    block
      real(8) mysgsdy, H1, H2
      mysgsdy  = 0.5d0 * (mut(i,j,k) + mut(i,j+1,k)) * dy(j)
      H1 = Cp * T(i,j,k)   + 0.5d0 * (u(idx)**2   + v(idx)**2   + w(idx)**2)   + qc2(i,j,k)
      H2 = Cp * T(i,j+1,k) + 0.5d0 * (u(idx+1)**2 + v(idx+1)**2 + w(idx+1)**2) + qc2(i,j+1,k)
      Hsgs = -mysgsdy * (-H1 + H2) / Prt
      muysgs = mysgsdy * (-u(idx) + u(idx+1))
      mvysgs = mysgsdy * (-v(idx) + v(idx+1))
      mwysgs = mysgsdy * (-w(idx) + w(idx+1))
    end block
    block
      real(8) mx1, mx2
      mx1    = 0.0625d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j+1,k) + mu(i,  j+1,k))
      mx2    = 0.0625d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j+1,k) + mu(i+1,j+1,k))
      mux    = (mx1    * (-Q(i-1,2,j,k) - Q(i-1,2,j+1,k)) + (mx1 - mx2) * (u(idx) + u(idx+1)) &
              + mx2    * ( Q(i+1,2,j,k) + Q(i+1,2,j+1,k))) * dx(i)
      muxsgs = (mxsgs1 * (-Q(i-1,2,j,k) - Q(i-1,2,j+1,k)) + (mxsgs1 - mxsgs2) * (u(idx) + u(idx+1)) &
              + mxsgs2 * ( Q(i+1,2,j,k) + Q(i+1,2,j+1,k))) * dx(i)
      mvx    = (mx1    * (-Q(i-1,3,j,k) - Q(i-1,3,j+1,k)) + (mx1 - mx2) * (v(idx) + v(idx+1)) &
              + mx2    * ( Q(i+1,3,j,k) + Q(i+1,3,j+1,k))) * dx(i)
      mvxsgs = (mxsgs1 * (-Q(i-1,3,j,k) - Q(i-1,3,j+1,k)) + (mxsgs1 - mxsgs2) * (v(idx) + v(idx+1)) &
              + mxsgs2 * ( Q(i+1,3,j,k) + Q(i+1,3,j+1,k))) * dx(i)
    end block
    block
      real(8) mz1, mz2
      mz1    = 0.0625d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i,j+1,k-1) + mu(i,j+1,k  ))
      mz2    = 0.0625d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i,j+1,k  ) + mu(i,j+1,k+1))
      mvz    = (mz1    * (-Q(i,3,j,k-1) - Q(i,3,j+1,k-1)) + (mz1 - mz2) * (v(idx) + v(idx+1)) &
              + mz2    * ( Q(i,3,j,k+1) + Q(i,3,j+1,k+1))) * dz(k)
      mvzsgs = (mzsgs1 * (-Q(i,3,j,k-1) - Q(i,3,j+1,k-1)) + (mzsgs1 - mzsgs2) * (v(idx) + v(idx+1)) &
              + mzsgs2 * ( Q(i,3,j,k+1) + Q(i,3,j+1,k+1))) * dz(k)
      mwz    = (mz1    * (-Q(i,4,j,k-1) - Q(i,4,j+1,k-1)) + (mz1 - mz2) * (w(idx) + w(idx+1)) &
              + mz2    * ( Q(i,4,j,k+1) + Q(i,4,j+1,k+1))) * dz(k)
      mwzsgs = (mzsgs1 * (-Q(i,4,j,k-1) - Q(i,4,j+1,k-1)) + (mzsgs1 - mzsgs2) * (w(idx) + w(idx+1)) &
              + mzsgs2 * ( Q(i,4,j,k+1) + Q(i,4,j+1,k+1))) * dz(k)
    end block
    tyx = muy + mvx
    tyy = two_third * (2.d0 * mvy - mwz - mux)
    tyz = mvz + mwy
    viscous_work = viscous_work + 0.5d0 * ((u(idx) + u(idx+1)) * tyx + (v(idx) + v(idx+1)) * tyy + (w(idx) + w(idx+1)) * tyz)
    tyx = tyx + muysgs + mvxsgs
    tyy = tyy + two_third * (2.d0 * mvysgs - mwzsgs - muxsgs)
    tyz = tyz + mvzsgs + mwysgs
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (viscous_work + Hsgs)
  end subroutine calc_Fv_LES2_koff


  !> Like calc_Gv_LES2 but restricted to k in [k_lo, k_hi]; uses load_smem_visc2_z_koff.
  attributes(global) subroutine calc_Gv_LES2_koff(nx, ny, nz, dx, dy, dz, Q, T, mu, mut, qc2, G, k_lo, k_hi)
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: dx(nx-1), dy(ny-1), dz(nz-1)
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: G(5,nx-2,ny-2,nz-1)
    integer, parameter :: sx = threadsGv%x
    integer, parameter :: sy = threadsGv%y
    integer, parameter :: sz = threadsGv%z + 1
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, Hsgs, tzx, tzy, tzz
    real(8) muz, muzsgs, mvz, mvzsgs, mwz, mwzsgs
    real(8) mwx, mwxsgs, mux, muxsgs, mvy, mvysgs, mwy, mwysgs
    real(8) mxsgs1, mxsgs2, mysgs1, mysgs2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt
    idx = (kt-1) + (jt-1)*sz + (it-1)*sz*sy
    call load_smem_visc2_z_koff(it, jt, kt, i, j, idx, nx, ny, nz, Q, u, v, w, k_lo)
    if (nx-1 < i .or. ny-1 < j .or. k_hi < k) return
    block
      real(8) mudz
      mudz         = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1)) * dz(k)
      viscous_work = Cp_over_Pr * mudz * (-T(i,j,k) + T(i,j,k+1))
      muz          = mudz * (-u(idx) + u(idx+1))
      mvz          = mudz * (-v(idx) + v(idx+1))
      mwz          = mudz * (-w(idx) + w(idx+1))
    end block
    mxsgs1 = 0.0625d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j,k+1) + mut(i,j,k+1))
    mxsgs2 = 0.0625d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j,k+1) + mut(i+1,j,k+1))
    mysgs1 = 0.0625d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i,j-1,k+1) + mut(i,j,k+1))
    mysgs2 = 0.0625d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i,j,k+1) + mut(i,j+1,k+1))
    block
      real(8) mzsgsdz, H1, H2
      mzsgsdz  = 0.5d0 * (mut(i,j,k) + mut(i,j,k+1)) * dz(k)
      H1 = Cp * T(i,j,k)   + 0.5d0 * (u(idx)**2   + v(idx)**2   + w(idx)**2)   + qc2(i,j,k)
      H2 = Cp * T(i,j,k+1) + 0.5d0 * (u(idx+1)**2 + v(idx+1)**2 + w(idx+1)**2) + qc2(i,j,k+1)
      Hsgs = -mzsgsdz * (-H1 + H2) / Prt
      muzsgs = mzsgsdz * (-u(idx) + u(idx+1))
      mvzsgs = mzsgsdz * (-v(idx) + v(idx+1))
      mwzsgs = mzsgsdz * (-w(idx) + w(idx+1))
    end block
    block
      real(8) mx1, mx2
      mx1    = 0.0625d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j,k+1) + mu(i,  j,k+1))
      mx2    = 0.0625d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j,k+1) + mu(i+1,j,k+1))
      mux    = (mx1    * (-Q(i-1,2,j,k) - Q(i-1,2,j,k+1)) + (mx1 - mx2) * (u(idx) + u(idx+1)) &
              + mx2    * ( Q(i+1,2,j,k) + Q(i+1,2,j,k+1))) * dx(i)
      muxsgs = (mxsgs1 * (-Q(i-1,2,j,k) - Q(i-1,2,j,k+1)) + (mxsgs1 - mxsgs2) * (u(idx) + u(idx+1)) &
              + mxsgs2 * ( Q(i+1,2,j,k) + Q(i+1,2,j,k+1))) * dx(i)
      mwx    = (mx1    * (-Q(i-1,4,j,k) - Q(i-1,4,j,k+1)) + (mx1 - mx2) * (w(idx) + w(idx+1)) &
              + mx2    * ( Q(i+1,4,j,k) + Q(i+1,4,j,k+1))) * dx(i)
      mwxsgs = (mxsgs1 * (-Q(i-1,4,j,k) - Q(i-1,4,j,k+1)) + (mxsgs1 - mxsgs2) * (w(idx) + w(idx+1)) &
              + mxsgs2 * ( Q(i+1,4,j,k) + Q(i+1,4,j,k+1))) * dx(i)
    end block
    block
      real(8) my1, my2
      my1    = 0.0625d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i,j-1,k+1) + mu(i,j,  k+1))
      my2    = 0.0625d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i,j,  k+1) + mu(i,j+1,k+1))
      mvy    = (my1    * (-Q(i,3,j-1,k) - Q(i,3,j-1,k+1)) + (my1 - my2) * (v(idx) + v(idx+1)) &
              + my2    * ( Q(i,3,j+1,k) + Q(i,3,j+1,k+1))) * dy(j)
      mvysgs = (mysgs1 * (-Q(i,3,j-1,k) - Q(i,3,j-1,k+1)) + (mysgs1 - mysgs2) * (v(idx) + v(idx+1)) &
              + mysgs2 * ( Q(i,3,j+1,k) + Q(i,3,j+1,k+1))) * dy(j)
      mwy    = (my1    * (-Q(i,4,j-1,k) - Q(i,4,j-1,k+1)) + (my1 - my2) * (w(idx) + w(idx+1)) &
              + my2    * ( Q(i,4,j+1,k) + Q(i,4,j+1,k+1))) * dy(j)
      mwysgs = (mysgs1 * (-Q(i,4,j-1,k) - Q(i,4,j-1,k+1)) + (mysgs1 - mysgs2) * (w(idx) + w(idx+1)) &
              + mysgs2 * ( Q(i,4,j+1,k) + Q(i,4,j+1,k+1))) * dy(j)
    end block
    tzx = mwx + muz
    tzy = mvz + mwy
    tzz = two_third * (2.d0 * mwz - mux - mvy)
    viscous_work = viscous_work + 0.5d0 * ((u(idx) + u(idx+1)) * tzx + (v(idx) + v(idx+1)) * tzy + (w(idx) + w(idx+1)) * tzz)
    tzx = tzx + mwxsgs + muzsgs
    tzy = tzy + mvzsgs + mwysgs
    tzz = tzz + two_third * (2.d0 * mwzsgs - muxsgs - mvysgs)
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (viscous_work + Hsgs)
  end subroutine calc_Gv_LES2_koff
end module calc_visc2
