module calc_me4_base
  use mod_globals, only : threadsEv, threadsFv, threadsGv
  use calc_visc_base
  implicit none
contains
  attributes(device) subroutine calc_derivative_x4(nx, ny, nz, i, j, k, dx, dy, dz, mu, u, v, w, txx, txy, txz, utxx, vtxy, wtxz)
    integer, intent(in), value :: nx, ny, nz, i, j, k
    real(8), intent(in), value :: dx, dy, dz
    real(8), intent(in)        :: mu(3)
    real(8), intent(in)        :: u(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), intent(in)        :: v(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), intent(in)        :: w(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), intent(out)       :: txx, txy, txz, utxx, vtxy, wtxz
    real(8), dimension(3) :: ux3, vx3, wx3, uy3, uz3, vy3, wz3, u3, v3, w3
    call calc_me4_base_x(nx, ny, nz, i, j, k, dx, dy, dz, u, v, w, ux3, vy3, wz3, uy3, vx3, wx3, uz3, u3, v3, w3)
    call tauxx(mu(:), ux3(:), vy3(:), wz3(:), u3(:), txx, utxx)
    call tauxy(mu(:), uy3(:), vx3(:), v3(:), txy, vtxy)
    call tauxy(mu(:), wx3(:), uz3(:), w3(:), txz, wtxz)
  end subroutine calc_derivative_x4

  attributes(device) subroutine calc_derivative_sgs_x4(nx, ny, nz, i, j, k, dx, dy, dz, mu, mt, u, v, w, txx_txxsgs, txy_txysgs, txz_txzsgs, utxx, vtxy, wtxz)
    integer, intent(in), value :: nx, ny, nz, i, j, k
    real(8), intent(in), value :: dx, dy, dz
    real(8), intent(in)        :: mu(3), mt(3)
    real(8), intent(in)        :: u(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), intent(in)        :: v(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), intent(in)        :: w(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), intent(out)       :: txx_txxsgs, txy_txysgs, txz_txzsgs, utxx, vtxy, wtxz
    real(8), dimension(3) :: ux3, vx3, wx3, uy3, uz3, vy3, wz3, u3, v3, w3
    real(8) txx, txxsgs, txy, txysgs, txz, txzsgs
    call calc_me4_base_x(nx, ny, nz, i, j, k, dx, dy, dz, u, v, w, ux3, vy3, wz3, uy3, vx3, wx3, uz3, u3, v3, w3)
    call tauxx(mu(:), ux3(:), vy3(:), wz3(:), u3(:), txx, utxx)
    call tauxy(mu(:), uy3(:), vx3(:), v3(:), txy, vtxy)
    call tauxy(mu(:), wx3(:), uz3(:), w3(:), txz, wtxz)
    call tauxx(mt(:), ux3(:), vy3(:), wz3(:), u3(:), txxsgs)
    call tauxy(mt(:), uy3(:), vx3(:), v3(:), txysgs)
    call tauxy(mt(:), wx3(:), uz3(:), w3(:), txzsgs)
    txx_txxsgs = txx + txxsgs
    txy_txysgs = txy + txysgs
    txz_txzsgs = txz + txzsgs
  end subroutine calc_derivative_sgs_x4

  attributes(device) subroutine calc_derivative_y4(nx, ny, nz, i, j, k, dx, dy, dz, mu, u, v, w, tyy, tyx, tyz, vtyy, utyx, wtyz)
    integer, intent(in), value :: nx, ny, nz, i, j, k
    real(8), intent(in), value :: dx, dy, dz
    real(8), intent(in)        :: mu(3)
    real(8), intent(in)        :: u(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), intent(in)        :: v(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), intent(in)        :: w(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), intent(out)       :: tyy, tyx, tyz, vtyy, utyx, wtyz
    real(8), dimension(3) :: vy3, wz3, ux3, uy3, vx3, vz3, wy3, u3, v3, w3
    call calc_me4_base_y(nx, ny, nz, i, j, k, dx, dy, dz, u, v, w, vy3, wz3, ux3, uy3, vx3, vz3, wy3, u3, v3, w3)
    call tauxx(mu(:), vy3(:), wz3(:), ux3(:), v3(:), tyy, vtyy)
    call tauxy(mu(:), uy3(:), vx3(:), u3(:), tyx, utyx)
    call tauxy(mu(:), vz3(:), wy3(:), w3(:), tyz, wtyz)
  end subroutine calc_derivative_y4
  
  attributes(device) subroutine calc_derivative_sgs_y4(nx, ny, nz, i, j, k, dx, dy, dz, mu, mt, u, v, w, tyy_tyysgs, tyx_tyxsgs, tyz_tyzsgs, vtyy, utyx, wtyz)
    integer, intent(in), value :: nx, ny, nz, i, j, k
    real(8), intent(in), value :: dx, dy, dz
    real(8), intent(in)        :: mu(3), mt(3)
    real(8), intent(in)        :: u(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), intent(in)        :: v(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), intent(in)        :: w(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), intent(out)       :: tyy_tyysgs, tyx_tyxsgs, tyz_tyzsgs, vtyy, utyx, wtyz
    real(8), dimension(3) :: vy3, wz3, ux3, uy3, vx3, vz3, wy3, u3, v3, w3
    real(8) tyy, tyysgs, tyx, tyxsgs, tyz, tyzsgs
    call calc_me4_base_y(nx, ny, nz, i, j, k, dx, dy, dz, u, v, w, vy3, wz3, ux3, uy3, vx3, vz3, wy3, u3, v3, w3)
    call tauxx(mu(:), vy3(:), wz3(:), ux3(:), v3(:), tyy, vtyy)
    call tauxy(mu(:), uy3(:), vx3(:), u3(:), tyx, utyx)
    call tauxy(mu(:), vz3(:), wy3(:), w3(:), tyz, wtyz)
    call tauxx(mt(:), vy3(:), wz3(:), ux3(:), v3(:), tyysgs)
    call tauxy(mt(:), uy3(:), vx3(:), u3(:), tyxsgs)
    call tauxy(mt(:), vz3(:), wy3(:), w3(:), tyzsgs)
    tyx_tyxsgs = tyx + tyxsgs
    tyy_tyysgs = tyy + tyysgs
    tyz_tyzsgs = tyz + tyzsgs
  end subroutine calc_derivative_sgs_y4
  
  attributes(device) subroutine calc_derivative_z4(nx, ny, nz, i, j, k, dx, dy, dz, mu, u, v, w, tzz, tzx, tzy, wtzz, utzx, vtzy)
    integer, intent(in), value :: nx, ny, nz, i, j, k
    real(8), intent(in), value :: dx, dy, dz
    real(8), intent(in)        :: mu(3)
    real(8), intent(in)        :: u(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), intent(in)        :: v(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), intent(in)        :: w(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), intent(out)       :: tzz, tzx, tzy, wtzz, utzx, vtzy
    real(8), dimension(3) :: wz3, ux3, vy3, wx3, uz3, vz3, wy3, u3, v3, w3
    call calc_me4_base_z(nx, ny, nz, i, j, k, dx, dy, dz, u, v, w, wz3, ux3, vy3, wx3, uz3, vz3, wy3, u3, v3, w3)
    call tauxx(mu(:), wz3(:), ux3(:), vy3(:), w3(:), tzz, wtzz)
    call tauxy(mu(:), wx3(:), uz3(:), u3(:), tzx, utzx)
    call tauxy(mu(:), vz3(:), wy3(:), v3(:), tzy, vtzy)
  end subroutine calc_derivative_z4
  
  attributes(device) subroutine calc_derivative_sgs_z4(nx, ny, nz, i, j, k, dx, dy, dz, mu, mt, u, v, w, tzz_tzzsgs, tzx_tzxsgs, tzy_tzysgs, wtzz, utzx, vtzy)
    integer, intent(in), value :: nx, ny, nz, i, j, k
    real(8), intent(in), value :: dx, dy, dz
    real(8), intent(in)        :: mu(3), mt(3)
    real(8), intent(in)        :: u(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), intent(in)        :: v(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), intent(in)        :: w(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), intent(out)       :: tzz_tzzsgs, tzx_tzxsgs, tzy_tzysgs, wtzz, utzx, vtzy
    real(8), dimension(3) :: wz3, ux3, vy3, wx3, uz3, vz3, wy3, u3, v3, w3
    real(8) tzz, tzzsgs, tzx, tzxsgs, tzy, tzysgs
    call calc_me4_base_z(nx, ny, nz, i, j, k, dx, dy, dz, u, v, w, wz3, ux3, vy3, wx3, uz3, vz3, wy3, u3, v3, w3)
    call tauxx(mu(:), wz3(:), ux3(:), vy3(:), w3(:), tzz, wtzz)
    call tauxy(mu(:), wx3(:), uz3(:), u3(:), tzx, utzx)
    call tauxy(mu(:), vz3(:), wy3(:), v3(:), tzy, vtzy)
    call tauxx(mt(:), wz3(:), ux3(:), vy3(:), w3(:), tzzsgs)
    call tauxy(mt(:), wx3(:), uz3(:), u3(:), tzxsgs)
    call tauxy(mt(:), vz3(:), wy3(:), v3(:), tzysgs)
    tzx_tzxsgs = tzx + tzxsgs
    tzy_tzysgs = tzy + tzysgs
    tzz_tzzsgs = tzz + tzzsgs
  end subroutine calc_derivative_sgs_z4

  attributes(device) subroutine calc_me4_base_x(nx, ny, nz, i, j, k, dx, dy, dz, u, v, w, ux3, vy3, wz3, uy3, vx3, wx3, uz3, u3, v3, w3)
    integer, intent(in), value :: nx, ny, nz, i, j, k
    real(8), intent(in), value :: dx, dy, dz
    real(8), intent(in)        :: u(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), intent(in)        :: v(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), intent(in)        :: w(-1:threadsEv%x+3,-1:threadsEv%y+2,-1:threadsEv%z+2)
    real(8), intent(out), dimension(3) :: ux3, vy3, wz3, uy3, vx3, wx3, uz3, u3, v3, w3
    real(8), dimension(6) :: uy6, uz6, vy6, wz6
    integer off
    ux3(:) = dx6(u(i-2:i+3,j,k), dx)
    vx3(:) = dx6(v(i-2:i+3,j,k), dx)
    wx3(:) = dx6(w(i-2:i+3,j,k), dx)
    do off = -2, 3
      uy6(3+off) = (2.d0 * (-u(i+off,j-1,k) + u(i+off,j+1,k)) - 0.25d0 * (-u(i+off,j-2,k) + u(i+off,j+2,k))) * dy / 3.d0
      uz6(3+off) = (2.d0 * (-u(i+off,j,k-1) + u(i+off,j,k+1)) - 0.25d0 * (-u(i+off,j,k-2) + u(i+off,j,k+2))) * dz / 3.d0
    enddo
    uy3(:) = interpolation6(uy6(:))
    uz3(:) = interpolation6(uz6(:))
    do off = -2, 3
      vy6(3+off) = (2.d0 * (-v(i+off,j-1,k) + v(i+off,j+1,k)) - 0.25d0 * (-v(i+off,j-2,k) + v(i+off,j+2,k))) * dy / 3.d0
    enddo
    vy3(:) = interpolation6(vy6(:))
    do off = -2, 3
      wz6(3+off) = (2.d0 * (-w(i+off,j,k-1) + w(i+off,j,k+1)) - 0.25d0 * (-w(i+off,j,k-2) + w(i+off,j,k+2))) * dz / 3.d0
    enddo
    wz3(:) = interpolation6(wz6(:))
    u3(:) = interpolation6(u(i-2:i+3,j,k))
    v3(:) = interpolation6(v(i-2:i+3,j,k))
    w3(:) = interpolation6(w(i-2:i+3,j,k))
  end subroutine calc_me4_base_x
  
  attributes(device) subroutine calc_me4_base_y(nx, ny, nz, i, j, k, dx, dy, dz, u, v, w, vy3, wz3, ux3, uy3, vx3, vz3, wy3, u3, v3, w3)
    integer, intent(in), value :: nx, ny, nz, i, j, k
    real(8), intent(in), value :: dx, dy, dz
    real(8), intent(in)        :: u(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), intent(in)        :: v(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), intent(in)        :: w(-1:threadsFv%x+2,-1:threadsFv%y+3,-1:threadsFv%z+2)
    real(8), intent(out), dimension(3) :: vy3, wz3, ux3, uy3, vx3, vz3, wy3, u3, v3, w3
    real(8), dimension(6) :: vx6, vz6, ux6, wz6, u6, v6, w6
    integer off
    u6(:)  = u(i,j-2:j+3,k)
    v6(:)  = v(i,j-2:j+3,k)
    w6(:)  = w(i,j-2:j+3,k)
    uy3(:) = dx6(u6(:), dy)
    vy3(:) = dx6(v6(:), dy)
    wy3(:) = dx6(w6(:), dy)
    do off = -2, 3
      vx6(3+off) = (2.d0 * (-v(i-1,j+off,k) + v(i+1,j+off,k)) - 0.25d0 * (-v(i-2,j+off,k) + v(i+2,j+off,k))) * dx / 3.d0
      vz6(3+off) = (2.d0 * (-v(i,j+off,k-1) + v(i,j+off,k+1)) - 0.25d0 * (-v(i,j+off,k-2) + v(i,j+off,k+2))) * dz / 3.d0
    enddo
    vx3(:) = interpolation6(vx6(:))
    vz3(:) = interpolation6(vz6(:))
    do off = -2, 3
      ux6(3+off) = (2.d0 * (-u(i-1,j+off,k) + u(i+1,j+off,k)) - 0.25d0 * (-u(i-2,j+off,k) + u(i+2,j+off,k))) * dx / 3.d0
    enddo
    ux3(:) = interpolation6(ux6(:))
    do off = -2, 3
      wz6(3+off) = (2.d0 * (-w(i,j+off,k-1) + w(i,j+off,k+1)) - 0.25d0 * (-w(i,j+off,k-2) + w(i,j+off,k+2))) * dz / 3.d0
    enddo
    wz3(:) = interpolation6(wz6(:))
    u3(:)  = interpolation6(u6(:))
    v3(:)  = interpolation6(v6(:))
    w3(:)  = interpolation6(w6(:))
  end subroutine calc_me4_base_y

  attributes(device) subroutine calc_me4_base_z(nx, ny, nz, i, j, k, dx, dy, dz, u, v, w, wz3, ux3, vy3, wx3, uz3, vz3, wy3, u3, v3, w3)
    integer, intent(in), value :: nx, ny, nz, i, j, k
    real(8), intent(in), value :: dx, dy, dz
    real(8), intent(in)        :: u(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), intent(in)        :: v(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), intent(in)        :: w(-1:threadsGv%x+2,-1:threadsGv%y+2,-1:threadsGv%z+3)
    real(8), intent(out), dimension(3) :: wz3, ux3, vy3, wx3, uz3, vz3, wy3, u3, v3, w3
    real(8), dimension(6) :: wx6, wy6, ux6, vy6, u6, v6, w6
    integer off
    u6(:)  = u(i,j,k-2:k+3)
    v6(:)  = v(i,j,k-2:k+3)
    w6(:)  = w(i,j,k-2:k+3)
    uz3(:) = dx6(u6(:), dz)
    vz3(:) = dx6(v6(:), dz)
    wz3(:) = dx6(w6(:), dz)
    do off = -2, 3
      wx6(3+off) = (2.d0 * (-w(i-1,j,k+off) + w(i+1,j,k+off)) - 0.25d0 * (-w(i-2,j,k+off) + w(i+2,j,k+off))) * dx / 3.d0
      wy6(3+off) = (2.d0 * (-w(i,j-1,k+off) + w(i,j+1,k+off)) - 0.25d0 * (-w(i,j-2,k+off) + w(i,j+2,k+off))) * dy / 3.d0
    enddo
    wx3(:) = interpolation6(wx6(:))
    wy3(:) = interpolation6(wy6(:))
    do off = -2, 3
      ux6(3+off) = (2.d0 * (-u(i-1,j,k+off) + u(i+1,j,k+off)) - 0.25d0 * (-u(i-2,j,k+off) + u(i+2,j,k+off))) * dx / 3.d0
    enddo
    ux3(:) = interpolation6(ux6(:))
    do off = -2, 3
      vy6(3+off) = (2.d0 * (-v(i,j-1,k+off) + v(i,j+1,k+off)) - 0.25d0 * (-v(i,j-2,k+off) + v(i,j+2,k+off))) * dy / 3.d0
    enddo
    vy3(:) = interpolation6(vy6(:))
    u3(:)  = interpolation6(u6(:))
    v3(:)  = interpolation6(v6(:))
    w3(:)  = interpolation6(w6(:))
  end subroutine calc_me4_base_z
end module calc_me4_base

