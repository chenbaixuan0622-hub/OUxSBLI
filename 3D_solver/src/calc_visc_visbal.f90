module calc_visc_visbal
  use mod_globals, only : threadsEv, threadsFv, threadsGv
  implicit none
contains
  attributes(device) subroutine calc_derivative_x2(nx, ny, nz, i, j, k, dx, dy, dz, mx, my, mz, u, v, w, txx, txy, txz, utxx, vtxy, wtxz)
    integer, intent(in), value :: nx, ny, nz, i, j, k
    real(8), intent(in), value :: dx, dy, dz, mx
    real(8), intent(in)        :: my(2), mz(2)
    real(8), intent(in)        :: u(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), intent(in)        :: v(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), intent(in)        :: w(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), intent(out)       :: txx, txy, txz, utxx, vtxy, wtxz
    real(8) mux, mvx, mwx, muy, mvy, muz, mwz
    mux = mx * (-u(i,j,k) + u(i+1,j,k)) * dx
    mvx = mx * (-v(i,j,k) + v(i+1,j,k)) * dx
    mwx = mx * (-w(i,j,k) + w(i+1,j,k)) * dx
    muy = 0.25d0 * (my(1) * (-u(i,j-1,k) + u(i,j,k) - u(i+1,j-1,k) + u(i+1,j,k)) &
                  + my(2) * (-u(i,j,k) + u(i,j+1,k) - u(i+1,j,k) + u(i+1,j+1,k))) * dy
    mvy = 0.25d0 * (my(1) * (-v(i,j-1,k) + v(i,j,k) - v(i+1,j-1,k) + v(i+1,j,k)) &
                  + my(2) * (-v(i,j,k) + v(i,j+1,k) - v(i+1,j,k) + v(i+1,j+1,k))) * dy
    muz = 0.25d0 * (mz(1) * (-u(i,j,k-1) + u(i,j,k) - u(i+1,j,k-1) + u(i+1,j,k)) &
                  + mz(2) * (-u(i,j,k) + u(i,j,k+1) - u(i+1,j,k) + u(i+1,j,k+1))) * dz
    mwz = 0.25d0 * (mz(1) * (-w(i,j,k-1) + w(i,j,k) - w(i+1,j,k-1) + w(i+1,j,k)) &
                  + mz(2) * (-w(i,j,k) + w(i,j,k+1) - w(i+1,j,k) + w(i+1,j,k+1))) * dz
    txx  = 2.d0 * (2.d0 * mux - mvy - mwz) / 3.d0
    txy  = muy + mvx
    txz  = mwx + muz
    utxx = 0.5d0 * (u(i,j,k) + u(i+1,j,k)) * txx
    vtxy = 0.5d0 * (v(i,j,k) + v(i+1,j,k)) * txy
    wtxz = 0.5d0 * (w(i,j,k) + w(i+1,j,k)) * txz
  end subroutine calc_derivative_x2

  attributes(device) subroutine calc_derivative_sgs_x2(nx, ny, nz, i, j, k, dx, dy, dz, mx, mtx, my, mty, mz, mtz, u, v, w, txx_txxsgs, txy_txysgs, txz_txzsgs, utxx, vtxy, wtxz)
    integer, intent(in), value :: nx, ny, nz, i, j, k
    real(8), intent(in), value :: dx, dy, dz, mx, mtx
    real(8), intent(in)        :: my(2), mty(2), mz(2), mtz(2)
    real(8), intent(in)        :: u(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), intent(in)        :: v(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), intent(in)        :: w(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), intent(out)       :: txx_txxsgs, txy_txysgs, txz_txzsgs, utxx, vtxy, wtxz
    real(8) mux, muxsgs, mvx, mvxsgs, mwx, mwxsgs, muy, muysgs, mvy, mvysgs, muz, muzsgs, mwz, mwzsgs
    real(8) txx, txxsgs, txy, txysgs, txz, txzsgs, tmp1, tmp2
    tmp1 = (-u(i,j,k) + u(i+1,j,k)) * dx
    mux    = mx  * tmp1
    muxsgs = mtx * tmp1
    tmp1 = (-v(i,j,k) + v(i+1,j,k)) * dx
    mvx    = mx  * tmp1
    mvxsgs = mtx * tmp1
    tmp1 = (-w(i,j,k) + w(i+1,j,k)) * dx
    mwx    = mx  * tmp1
    mwxsgs = mtx * tmp1
    tmp1 = -u(i,j-1,k) + u(i,j,k) - u(i+1,j-1,k) + u(i+1,j,k)
    tmp2 = -u(i,j,k) + u(i,j+1,k) - u(i+1,j,k) + u(i+1,j+1,k)
    muy    = 0.25d0 * (my(1)  * tmp1 + my(2)  * tmp2) * dy
    muysgs = 0.25d0 * (mty(1) * tmp1 + mty(2) * tmp2) * dy
    
    tmp1 = -v(i,j-1,k) + v(i,j,k) - v(i+1,j-1,k) + v(i+1,j,k)
    tmp2 = -v(i,j,k) + v(i,j+1,k) - v(i+1,j,k) + v(i+1,j+1,k)
    mvy    = 0.25d0 * (my(1)  * tmp1 + my(2)  * tmp2) * dy
    mvysgs = 0.25d0 * (mty(1) * tmp1 + mty(2) * tmp2) * dy

    tmp1 = -u(i,j,k-1) + u(i,j,k) - u(i+1,j,k-1) + u(i+1,j,k)
    tmp2 = -u(i,j,k) + u(i,j,k+1) - u(i+1,j,k) + u(i+1,j,k+1)
    muz    = 0.25d0 * (mz(1)  * tmp1 + mz(2)  * tmp2) * dz
    muzsgs = 0.25d0 * (mtz(1) * tmp1 + mtz(2) * tmp2) * dz

    tmp1 = -w(i,j,k-1) + w(i,j,k) - w(i+1,j,k-1) + w(i+1,j,k)
    tmp2 = -w(i,j,k) + w(i,j,k+1) - w(i+1,j,k) + w(i+1,j,k+1)
    mwz    = 0.25d0 * (mz(1)  * tmp1 + mz(2)  * tmp2) * dz
    mwzsgs = 0.25d0 * (mtz(1) * tmp1 + mtz(2) * tmp2) * dz
    txx  = 2.d0 * (2.d0 * mux - mvy - mwz) / 3.d0
    txy  = muy + mvx
    txz  = mwx + muz
    txxsgs = 2.d0 * (2.d0 * muxsgs - mvysgs - mwzsgs) / 3.d0
    txysgs = muysgs + mvxsgs
    txzsgs = mwxsgs + muzsgs
    txx_txxsgs = txx + txxsgs
    txy_txysgs = txy + txysgs
    txz_txzsgs = txz + txzsgs
    utxx = 0.5d0 * (u(i,j,k) + u(i+1,j,k)) * txx
    vtxy = 0.5d0 * (v(i,j,k) + v(i+1,j,k)) * txy
    wtxz = 0.5d0 * (w(i,j,k) + w(i+1,j,k)) * txz
  end subroutine calc_derivative_sgs_x2

  attributes(device) subroutine calc_derivative_y2(nx, ny, nz, i, j, k, dx, dy, dz, my, mx, mz, u, v, w, tyy, tyx, tyz, vtyy, utyx, wtyz)
    integer, intent(in), value :: nx, ny, nz, i, j, k
    real(8), intent(in), value :: dx, dy, dz, my
    real(8), intent(in)        :: mx(2), mz(2)
    real(8), intent(in)        :: u(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), intent(in)        :: v(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), intent(in)        :: w(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), intent(out)       :: tyy, tyx, tyz, vtyy, utyx, wtyz
    real(8) muy, mvy, mwy, mux, mvx, mvz, mwz
    muy = my * (-u(i,j,k) + u(i,j+1,k)) * dy
    mvy = my * (-v(i,j,k) + v(i,j+1,k)) * dy
    mwy = my * (-w(i,j,k) + w(i,j+1,k)) * dy
    mux = 0.25d0 * (mx(1) * (-u(i-1,j,k) + u(i,j,k) - u(i-1,j+1,k) + u(i,j+1,k)) &
                  + mx(2) * (-u(i,j,k) + u(i+1,j,k) - u(i,j+1,k) + u(i+1,j+1,k))) * dx
    mvx = 0.25d0 * (mx(1) * (-v(i-1,j,k) + v(i,j,k) - v(i-1,j+1,k) + v(i,j+1,k)) &
                  + mx(2) * (-v(i,j,k) + v(i+1,j,k) - v(i,j+1,k) + v(i+1,j+1,k))) * dx
    mvz = 0.25d0 * (mz(1) * (-v(i,j,k-1) + v(i,j,k) - v(i,j+1,k-1) + v(i,j+1,k)) &
                  + mz(2) * (-v(i,j,k) + v(i,j,k+1) - v(i,j+1,k) + v(i,j+1,k+1))) * dz
    mwz = 0.25d0 * (mz(1) * (-w(i,j,k-1) + w(i,j,k) - w(i,j+1,k-1) + w(i,j+1,k)) &
                  + mz(2) * (-w(i,j,k) + w(i,j,k+1) - w(i,j+1,k) + w(i,j+1,k+1))) * dz
    tyy  = 2.d0 * (2.d0 * mvy - mwz - mux) / 3.d0
    tyx  = muy + mvx
    tyz  = mvz + mwy
    utyx = 0.5d0 * (u(i,j,k) + u(i,j+1,k)) * tyx
    vtyy = 0.5d0 * (v(i,j,k) + v(i,j+1,k)) * tyy
    wtyz = 0.5d0 * (w(i,j,k) + w(i,j+1,k)) * tyz
  end subroutine calc_derivative_y2

  attributes(device) subroutine calc_derivative_sgs_y2(nx, ny, nz, i, j, k, dx, dy, dz, my, mty, mx, mtx, mz, mtz, u, v, w, tyy_tyysgs, tyx_tyxsgs, tyz_tyzsgs, vtyy, utyx, wtyz)
    integer, intent(in), value :: nx, ny, nz, i, j, k
    real(8), intent(in), value :: dx, dy, dz, my, mty
    real(8), intent(in)        :: mx(2), mtx(2), mz(2), mtz(2)
    real(8), intent(in)        :: u(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), intent(in)        :: v(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), intent(in)        :: w(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), intent(out)       :: tyy_tyysgs, tyx_tyxsgs, tyz_tyzsgs, vtyy, utyx, wtyz
    real(8) muy, muysgs, mvy, mvysgs, mwy, mwysgs, mux, muxsgs, mvx, mvxsgs, mvz, mvzsgs, mwz, mwzsgs
    real(8) tyy, tyysgs, tyx, tyxsgs, tyz, tyzsgs, tmp1, tmp2
    tmp1 = (-u(i,j,k) + u(i,j+1,k)) * dy
    muy    = my  * tmp1
    muysgs = mty * tmp1
    tmp1 = (-v(i,j,k) + v(i,j+1,k)) * dy
    mvy    = my  * tmp1
    mvysgs = mty * tmp1
    tmp1 = (-w(i,j,k) + w(i,j+1,k)) * dy
    mwy    = my  * tmp1
    mwysgs = mty * tmp1
    tmp1 = -u(i-1,j,k) + u(i,j,k) - u(i-1,j+1,k) + u(i,j+1,k)
    tmp2 = -u(i,j,k) + u(i+1,j,k) - u(i,j+1,k) + u(i+1,j+1,k)
    mux    = 0.25d0 * (mx(1)  * tmp1 + mx(2)  * tmp2) * dx
    muxsgs = 0.25d0 * (mtx(1) * tmp1 + mtx(2) * tmp2) * dx
    
    tmp1 = -v(i-1,j,k) + v(i,j,k) - v(i-1,j+1,k) + v(i,j+1,k)
    tmp2 = -v(i,j,k) + v(i+1,j,k) - v(i,j+1,k) + v(i+1,j+1,k)
    mvx    = 0.25d0 * (mx(1)  * tmp1 + mx(2)  * tmp2) * dx
    mvxsgs = 0.25d0 * (mtx(1) * tmp1 + mtx(2) * tmp2) * dx
    
    tmp1 = -v(i,j,k-1) + v(i,j,k) - v(i,j+1,k-1) + v(i,j+1,k)
    tmp2 = -v(i,j,k) + v(i,j,k+1) - v(i,j+1,k) + v(i,j+1,k+1)
    mvz    = 0.25d0 * (mz(1)  * tmp1 + mz(2)  * tmp2) * dz
    mvzsgs = 0.25d0 * (mtz(1) * tmp1 + mtz(2) * tmp2) * dz
    
    tmp1 = -w(i,j,k-1) + w(i,j,k) - w(i,j+1,k-1) + w(i,j+1,k)
    tmp2 = -w(i,j,k) + w(i,j,k+1) - w(i,j+1,k) + w(i,j+1,k+1)
    mwz    = 0.25d0 * (mz(1)  * tmp1 + mz(2)  * tmp2) * dz
    mwzsgs = 0.25d0 * (mtz(1) * tmp1 + mtz(2) * tmp2) * dz
    tyy  = 2.d0 * (2.d0 * mvy - mwz - mux) / 3.d0
    tyx  = muy + mvx
    tyz  = mvz + mwy
    tyysgs = 2.d0 * (2.d0 * mvysgs - mwzsgs - muxsgs) / 3.d0
    tyxsgs = muysgs + mvxsgs
    tyzsgs = mvzsgs + mwysgs
    tyx_tyxsgs = tyx + tyxsgs
    tyy_tyysgs = tyy + tyysgs
    tyz_tyzsgs = tyz + tyzsgs
    utyx = 0.5d0 * (u(i,j,k) + u(i,j+1,k)) * tyx
    vtyy = 0.5d0 * (v(i,j,k) + v(i,j+1,k)) * tyy
    wtyz = 0.5d0 * (w(i,j,k) + w(i,j+1,k)) * tyz
  end subroutine calc_derivative_sgs_y2

  attributes(device) subroutine calc_derivative_z2(nx, ny, nz, i, j, k, dx, dy, dz, mz, mx, my, u, v, w, tzz, tzx, tzy, wtzz, utzx, vtzy)
    integer, intent(in), value :: nx, ny, nz, i, j, k
    real(8), intent(in), value :: dx, dy, dz, mz
    real(8), intent(in)        :: mx(2), my(2)
    real(8), intent(in)        :: u(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), intent(in)        :: v(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), intent(in)        :: w(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), intent(out)       :: tzz, tzx, tzy, wtzz, utzx, vtzy
    real(8) muz, mvz, mwz, mux, mwx, mvy, mwy
    muz = mz * (-u(i,j,k) + u(i,j,k+1)) * dz
    mvz = mz * (-v(i,j,k) + v(i,j,k+1)) * dz
    mwz = mz * (-w(i,j,k) + w(i,j,k+1)) * dz
    mux = 0.25d0 * (mx(1) * (-u(i-1,j,k) + u(i,j,k) - u(i-1,j,k+1) + u(i,j,k+1)) &
                  + mx(2) * (-u(i,j,k) + u(i+1,j,k) - u(i,j,k+1) + u(i+1,j,k+1))) * dx
    mwx = 0.25d0 * (mx(1) * (-w(i-1,j,k) + w(i,j,k) - w(i-1,j,k+1) + w(i,j,k+1)) &
                  + mx(2) * (-w(i,j,k) + w(i+1,j,k) - w(i,j,k+1) + w(i+1,j,k+1))) * dx
    mvy = 0.25d0 * (my(1) * (-v(i,j-1,k) + v(i,j,k) - v(i,j-1,k+1) + v(i,j,k+1)) &
                  + my(2) * (-v(i,j,k) + v(i,j+1,k) - v(i,j,k+1) + v(i,j+1,k+1))) * dy
    mwy = 0.25d0 * (my(1) * (-w(i,j-1,k) + w(i,j,k) - w(i,j-1,k+1) + w(i,j,k+1)) &
                  + my(2) * (-w(i,j,k) + w(i,j+1,k) - w(i,j,k+1) + w(i,j+1,k+1))) * dy
    tzz  = 2.d0 * (2.d0 * mwz - mux - mvy) / 3.d0
    tzx  = mwx + muz
    tzy  = mvz + mwy
    utzx = 0.5d0 * (u(i,j,k) + u(i,j,k+1)) * tzx
    vtzy = 0.5d0 * (v(i,j,k) + v(i,j,k+1)) * tzy
    wtzz = 0.5d0 * (w(i,j,k) + w(i,j,k+1)) * tzz
  end subroutine calc_derivative_z2
  
  attributes(device) subroutine calc_derivative_sgs_z2(nx, ny, nz, i, j, k, dx, dy, dz, mz, mtz, mx, mtx, my, mty, u, v, w, tzz_tzzsgs, tzx_tzxsgs, tzy_tzysgs, wtzz, utzx, vtzy)
    integer, intent(in), value :: nx, ny, nz, i, j, k
    real(8), intent(in), value :: dx, dy, dz, mz, mtz
    real(8), intent(in)        :: mx(2), mtx(2), my(2), mty(2)
    real(8), intent(in)        :: u(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), intent(in)        :: v(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), intent(in)        :: w(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), intent(out)       :: tzz_tzzsgs, tzx_tzxsgs, tzy_tzysgs, wtzz, utzx, vtzy
    real(8) muz, muzsgs, mvz, mvzsgs, mwz, mwzsgs, mux, muxsgs, mwx, mwxsgs, mvy, mvysgs, mwy, mwysgs
    real(8) tzz, tzzsgs, tzx, tzxsgs, tzy, tzysgs, tmp1, tmp2
    tmp1 = (-u(i,j,k) + u(i,j,k+1)) * dz
    muz    = mz  * tmp1
    muzsgs = mtz * tmp1
    tmp1 = (-v(i,j,k) + v(i,j,k+1)) * dz
    mvz    = mz  * tmp1
    mvzsgs = mtz * tmp1
    tmp1 = (-w(i,j,k) + w(i,j,k+1)) * dz
    mwz    = mz  * tmp1
    mwzsgs = mtz * tmp1
    tmp1 = -u(i-1,j,k) + u(i,j,k) - u(i-1,j,k+1) + u(i,j,k+1)
    tmp2 = -u(i,j,k) + u(i+1,j,k) - u(i,j,k+1) + u(i+1,j,k+1)
    mux    = 0.25d0 * (mx(1)  * tmp1 + mx(2)  * tmp2) * dx
    muxsgs = 0.25d0 * (mtx(1) * tmp1 + mtx(2) * tmp2) * dx
    
    tmp1 = -w(i-1,j,k) + w(i,j,k) - w(i-1,j,k+1) + w(i,j,k+1)
    tmp2 = -w(i,j,k) + w(i+1,j,k) - w(i,j,k+1) + w(i+1,j,k+1)
    mwx    = 0.25d0 * (mx(1)  * tmp1 + mx(2)  * tmp2) * dx
    mwxsgs = 0.25d0 * (mtx(1) * tmp1 + mtx(2) * tmp2) * dx
    
    tmp1 = -v(i,j-1,k) + v(i,j,k) - v(i,j-1,k+1) + v(i,j,k+1)
    tmp2 = -v(i,j,k) + v(i,j+1,k) - v(i,j,k+1) + v(i,j+1,k+1)
    mvy    = 0.25d0 * (my(1)  * tmp1 + my(2)  * tmp2) * dy
    mvysgs = 0.25d0 * (mty(1) * tmp1 + mty(2) * tmp2) * dy
    
    tmp1 = -w(i,j-1,k) + w(i,j,k) - w(i,j-1,k+1) + w(i,j,k+1)
    tmp2 = -w(i,j,k) + w(i,j+1,k) - w(i,j,k+1) + w(i,j+1,k+1)
    mwy    = 0.25d0 * (my(1)  * tmp1 + my(2)  * tmp2) * dy
    mwysgs = 0.25d0 * (mty(1) * tmp1 + mty(2) * tmp2) * dy
    tzz = 2.d0 * (2.d0 * mwz - mux - mvy) / 3.d0
    tzx = mwx + muz
    tzy = mvz + mwy
    tzzsgs = 2.d0 * (2.d0 * mwzsgs - muxsgs - mvysgs) / 3.d0
    tzxsgs = mwxsgs + muzsgs
    tzysgs = mvzsgs + mwysgs
    tzx_tzxsgs = tzx + tzxsgs
    tzy_tzysgs = tzy + tzysgs
    tzz_tzzsgs = tzz + tzzsgs
    utzx = 0.5d0 * (u(i,j,k) + u(i,j,k+1)) * tzx
    vtzy = 0.5d0 * (v(i,j,k) + v(i,j,k+1)) * tzy
    wtzz = 0.5d0 * (w(i,j,k) + w(i,j,k+1)) * tzz
  end subroutine calc_derivative_sgs_z2
end module calc_visc_visbal

