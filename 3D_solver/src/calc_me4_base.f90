module calc_me4_base
  use mod_globals, only : gamma, R, Pr, Prt
  use mod_constant, only : Cp, gamma_1
  use calc_sutherland, only : mu6, mu2, mu23, mu32
  use calc_visc_base
  implicit none
contains
  !dir$ inline
  attributes(device) function heat_conduction6(mu, T, dx) result(ans)
    real(8), intent(in), device :: mu(3), T(6)
    real(8), intent(in), value  :: dx
    real(8), dimension(3) ::  kTx
    real(8) :: ans
    kTx(:) = Cp * mu(:) * dx6(T(:), dx) / Pr
    ans = flux4(kTx(:))
  end function heat_conduction6
  
  attributes(device) subroutine calc_me4_base_x(nx, ny, nz, i, j, k, dx, dy, dz, Q, txx, txy, txz, utxx, vtxy, wtxz, kTx)
    integer, intent(in), value  :: nx, ny, nz, i, j, k
    real(8), intent(in), device :: dx(nx-1), dy(ny-1), dz(nz-1), Q(5,nx,ny,nz)
    real(8), intent(out)        :: txx, txy, txz, utxx, vtxy, wtxz, kTx
    real(8), dimension(6,5), device :: tmp4
    real(8), dimension(6), device   :: T6
    real(8), dimension(3), device   :: u3, v3, w3, ux3, vx3, wx3, uy3, vy3, uz3, wz3, mu
    ! dudx & dudy
    tmp4(:,:) = Q(2,i-2:i+3,j-2:j+2,k)
    uy3(:)    = dy65(tmp4(:,:), dy(j))
    ux3(:)    = dx6(tmp4(:,3), dx(i))
    u3(:)     = interpolation6(tmp4(:,3))
    ! dudz
    tmp4(:,:) = Q(2,i-2:i+3,j,k-2:k+2)
    uz3(:)    = dy65(tmp4(:,:), dz(k))
    ! dvdx & dvdy
    tmp4(:,:) = Q(3,i-2:i+3,j-2:j+2,k)
    vy3(:)    = dy65(tmp4(:,:), dy(j))
    vx3(:)    = dx6(tmp4(:,3), dx(i))
    v3(:)     = interpolation6(tmp4(:,3))
    ! dwdx & dwdz
    tmp4(:,:) = Q(4,i-2:i+3,j,k-2:k+2)
    wz3(:)    = dy65(tmp4(:,:), dz(k))
    wx3(:)    = dx6(tmp4(:,3), dx(i))
    w3(:)     = interpolation6(tmp4(:,3))
    ! dTdx
    T6(:)     = Q(5,i-2:i+3,j,k) / (R * Q(1,i-2:i+3,j,k))
    mu(:)     = mu6(T6(:))
    call tauxx_4(mu(:), ux3(:), vy3(:), wz3(:), u3(:), txx, utxx)
    call tauxy_4(mu(:), uy3(:), vx3(:), v3(:), txy, vtxy)
    call tauxy_4(mu(:), wx3(:), uz3(:), w3(:), txz, wtxz)
    kTx = heat_conduction6(mu(:), T6(:), dx(i))
  end subroutine calc_me4_base_x
  
  attributes(device) subroutine calc_me4_base_les_x(nx, ny, nz, i, j, k, dx, dy, dz, Q, mut, qc2, txx, txy, txz, utxx, vtxy, wtxz, kTx, Hsgs)
    integer, intent(in), value  :: nx, ny, nz, i, j, k
    real(8), intent(in), device :: dx(nx-1), dy(ny-1), dz(nz-1), Q(5,nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(out)        :: txx, txy, txz, utxx, vtxy, wtxz, kTx, Hsgs
    real(8), dimension(6,5), device :: tmp4
    real(8), dimension(6), device   :: T6
    real(8), dimension(4), device   :: H
    real(8), dimension(3), device   :: u3, v3, w3, ux3, vx3, wx3, uy3, vy3, uz3, wz3, mu
    real(8) mutx
    ! dudx & dudy
    tmp4(:,:) = Q(2,i-2:i+3,j-2:j+2,k)
    uy3(:)    = dy65(tmp4(:,:), dy(j))
    ux3(:)    = dx6(tmp4(:,3), dx(i))
    u3(:)     = interpolation6(tmp4(:,3))
    ! dudz
    tmp4(:,:) = Q(2,i-2:i+3,j,k-2:k+2)
    uz3(:)    = dy65(tmp4(:,:), dz(k))
    ! dvdx & dvdy
    tmp4(:,:) = Q(3,i-2:i+3,j-2:j+2,k)
    vy3(:)    = dy65(tmp4(:,:), dy(j))
    vx3(:)    = dx6(tmp4(:,3), dx(i))
    v3(:)     = interpolation6(tmp4(:,3))
    ! dwdx & dwdz
    tmp4(:,:) = Q(4,i-2:i+3,j,k-2:k+2)
    wz3(:)    = dy65(tmp4(:,:), dz(k))
    wx3(:)    = dx6(tmp4(:,3), dx(i))
    w3(:)     = interpolation6(tmp4(:,3))
    ! dTdx
    T6(:)     = Q(5,i-2:i+3,j,k) / (R * Q(1,i-2:i+3,j,k))
    mu(:)     = mu6(T6(:))
    call tauxx_4(mu(:), ux3(:), vy3(:), wz3(:), u3(:), txx, utxx)
    call tauxy_4(mu(:), uy3(:), vx3(:), v3(:), txy, vtxy)
    call tauxy_4(mu(:), wx3(:), uz3(:), w3(:), txz, wtxz)
    kTx  = heat_conduction6(mu(:), T6(:), dx(i))
    mutx = 0.0625d0 * (-mut(i-1,j,k) + 9.d0 * (mut(i,j,k) + mut(i+1,j,k)) -mut(i+2,j,k))
    txx  = txx + 2.d0 * mutx * (2.d0 * ux3(2) - vy3(2) - wz3(2)) / 3.d0
    txy  = txy + mutx * (uy3(2) + vx3(2))
    txz  = txz + mutx * (wx3(2) + uz3(2))
    H(:) = (gamma * Q(5,i-1:i+2,j,k) / (Q(1,i-1:i+2,j,k) * gamma_1)) &
           + 0.5d0 * (Q(2,i-1:i+2,j,k)**2 + Q(3,i-1:i+2,j,k)**2 + Q(4,i-1:i+2,j,k)**2) + qc2(i-1:i+2,j,k)
    Hsgs = -mutx * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) / 3.d0) * dx(i) / Prt
  end subroutine calc_me4_base_les_x
  
  attributes(device) subroutine calc_me4_base_y(nx, ny, nz, i, j, k, dx, dy, dz, Q, tyx, tyy, tyz, utyx, vtyy, wtyz, kTy)
    integer, intent(in), value  :: nx, ny, nz, i, j, k
    real(8), intent(in), device :: dx(nx-1), dy(ny-1), dz(nz-1), Q(5,nx,ny,nz)
    real(8), intent(out)        :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy
    ! 4th-order 2
    real(8), dimension(5,6), device :: tmp4
    real(8), dimension(6,5), device :: tmp5
    real(8), dimension(6), device   :: T6, u6, v6, w6
    real(8), dimension(3), device   :: uy3, vy3, wy3, vz3, wz3, ux3, vx3, mu
    ! dudx & dudy
    tmp4(:,:) = Q(2,i-2:i+2,j-2:j+3,k)
    ux3(:)    = dy56(tmp4(:,:), dx(i))
    u6(:)     = tmp4(3,:)
    uy3(:)    = dx6(u6(:), dy(j))
    ! dvdx dvdy
    tmp4(:,:) = Q(3,i-2:i+2,j-2:j+3,k)
    vx3(:)    = dy56(tmp4(:,:), dx(i))
    v6(:)     = tmp4(3,:)
    vy3(:)    = dx6(v6(:), dy(j))
    ! dvdz
    tmp5(:,:) = Q(3,i,j-2:j+3,k-2:k+2)
    vz3(:)    = dy65(tmp5(:,:), dz(k))
    ! dwdz & dwdy
    tmp5(:,:) = Q(4,i,j-2:j+3,k-2:k+2)
    wz3(:)    = dy65(tmp5(:,:), dz(k))
    w6(:)     = tmp5(:,3)
    wy3(:)    = dx6(w6(:), dy(j))
    ! dTdy
    T6(:)     = Q(5,i,j-2:j+3,k) / (R * Q(1,i,j-2:j+3,k))
    mu(:)     = mu6(T6(:))
    call tauxy4(mu(:), uy3(:), vx3(:), u6(:), tyx, utyx)
    call tauxx4(mu(:), vy3(:), wz3(:), ux3(:), v6(:), tyy, vtyy)
    call tauxy4(mu(:), vz3(:), wy3(:), w6(:), tyz, wtyz)
    kTy = heat_conduction6(mu(:), T6(:), dy(j))
  end subroutine calc_me4_base_y
  
  attributes(device) subroutine calc_me4_base_les_y(nx, ny, nz, i, j, k, dy, dx, dz, Q, mut, qc2, tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, Hsgs)
    integer, intent(in), value  :: nx, ny, nz, i, j, k
    real(8), intent(in), device :: dx(nx-1), dy(ny-1), dz(nz-1), Q(5,nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(out)        :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, Hsgs
    ! 4th-order 2
    real(8), dimension(5,6), device :: tmp4
    real(8), dimension(6,5), device :: tmp5
    real(8), dimension(6), device   :: T6, u6, v6, w6
    real(8), dimension(4), device   :: H
    real(8), dimension(3), device   :: uy3, vy3, wy3, vz3, wz3, ux3, vx3, mu
    real(8) muty
    ! dudx & dudy
    tmp4(:,:) = Q(2,i-2:i+2,j-2:j+3,k)
    ux3(:)    = dy56(tmp4(:,:), dx(i))
    u6(:)     = tmp4(3,:)
    uy3(:)    = dx6(u6(:), dy(j))
    ! dvdx dvdy
    tmp4(:,:) = Q(3,i-2:i+2,j-2:j+3,k)
    vx3(:)    = dy56(tmp4(:,:), dx(i))
    v6(:)     = tmp4(3,:)
    vy3(:)    = dx6(v6(:), dy(j))
    ! dvdz
    tmp5(:,:) = Q(3,i,j-2:j+3,k-2:k+2)
    vz3(:)    = dy65(tmp5(:,:), dz(k))
    ! dwdz & dwdy
    tmp5(:,:) = Q(4,i,j-2:j+3,k-2:k+2)
    wz3(:)    = dy65(tmp5(:,:), dz(k))
    w6(:)     = tmp5(:,3)
    wy3(:)    = dx6(w6(:), dy(j))
    ! dTdy
    T6(:)     = Q(5,i,j-2:j+3,k) / (R * Q(1,i,j-2:j+3,k))
    mu(:)     = mu6(T6(:))
    call tauxy4(mu(:), uy3(:), vx3(:), u6(:), tyx, utyx)
    call tauxx4(mu(:), vy3(:), wz3(:), ux3(:), v6(:), tyy, vtyy)
    call tauxy4(mu(:), vz3(:), wy3(:), w6(:), tyz, wtyz)
    kTy = heat_conduction6(mu(:), T6(:), dy(j))
    muty = 0.0625d0 * (-mut(i,j-1,k) + 9.d0 * (mut(i,j,k) + mut(i,j+1,k)) -mut(i,j+2,k))
    tyx  = tyx + muty * (uy3(2) + vx3(2))
    tyy  = tyy + 2.d0 * muty * (2.d0 * vy3(2) - ux3(2) - wz3(2)) / 3.d0
    tyz  = tyz + muty * (vz3(2) + wy3(2))
    H(:) = (gamma * Q(5,i,j-1:j+2,k) / (Q(1,i,j-1:j+2,k) * gamma_1)) &
           + 0.5d0 * (Q(2,i,j-1:j+2,k)**2 + Q(3,i,j-1:j+2,k)**2 + Q(4,i,j-1:j+2,k)**2) + qc2(i,j-1:j+2,k)
    Hsgs = -muty * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) / 3.d0) * dy(j) / Prt
  end subroutine calc_me4_base_les_y
  
  attributes(device) subroutine calc_me4_base_z(nx, ny, nz, i, j, k, dx, dy, dz, Q, tzx, tzy, tzz, utzx, vtzy, wtzz, kTz)
    integer, intent(in), value  :: nx, ny, nz, i, j, k
    real(8), intent(in), device :: dx(nx-1), dy(ny-1), dz(nz-1), Q(5,nx,ny,nz)
    real(8), intent(out)        :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz
    ! 4th-order 2
    real(8), dimension(5,6), device :: tmp4
    real(8), dimension(6), device   :: T6, u6, v6, w6
    real(8), dimension(3), device   :: uz3, vz3, wz3, wx3, ux3, vy3, wy3, mu
    ! dudx & dudz
    tmp4(:,:) = Q(2,i-2:i+2,j,k-2:k+3)
    ux3(:)    = dy56(tmp4(:,:), dx(i))
    u6(:)     = tmp4(3,:)
    uz3(:)    = dx6(u6(:), dz(k))
    ! dvdy & dvdz
    tmp4(:,:) = Q(3,i,j-2:j+2,k-2:k+3)
    vy3(:)    = dy56(tmp4(:,:), dy(j))
    v6(:)     = tmp4(3,:)
    vz3(:)    = dx6(v6(:), dz(k))
    ! dwdx
    tmp4(:,:) = Q(4,i-2:i+2,j,k-2:k+3)
    wx3(:)    = dy56(tmp4(:,:), dx(i))
    ! dwdy & dwdz
    tmp4(:,:) = Q(4,i,j-2:j+2,k-2:k+3)
    wy3(:)    = dy56(tmp4(:,:), dy(j))
    w6(:)     = tmp4(3,:)
    wz3(:)    = dx6(w6(:), dz(k))
    ! dTdz
    T6(:)     = Q(5,i,j,k-2:k+3) / (R * Q(1,i,j,k-2:k+3))
    mu(:)     = mu6(T6(:))
    call tauxy4(mu(:), wx3(:), uz3(:), u6(:), tzx, utzx)
    call tauxy4(mu(:), vz3(:), wy3(:), v6(:), tzy, vtzy)
    call tauxx4(mu(:), wz3(:), ux3(:), vy3(:), w6(:), tzz, wtzz)
    kTz = heat_conduction6(mu(:), T6(:), dz(k))
  end subroutine calc_me4_base_z

  attributes(device) subroutine calc_me4_base_les_z(nx, ny, nz, i, j, k, dx, dy, dz, Q, mut, qc2, tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, Hsgs)
    integer, intent(in), value  :: nx, ny, nz, i, j, k
    real(8), intent(in), device :: dx(nx-1), dy(ny-1), dz(nz-1), Q(5,nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(out)        :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, Hsgs
    ! 4th-order 2
    real(8), dimension(5,6), device :: tmp4
    real(8), dimension(6), device   :: T6, u6, v6, w6
    real(8), dimension(4), device   :: H
    real(8), dimension(3), device   :: uz3, vz3, wz3, wx3, ux3, vy3, wy3, mu
    real(8) mutz
    ! dudx & dudz
    tmp4(:,:) = Q(2,i-2:i+2,j,k-2:k+3)
    ux3(:)    = dy56(tmp4(:,:), dx(i))
    u6(:)     = tmp4(3,:)
    uz3(:)    = dx6(u6(:), dz(k))
    ! dvdy & dvdz
    tmp4(:,:) = Q(3,i,j-2:j+2,k-2:k+3)
    vy3(:)    = dy56(tmp4(:,:), dy(j))
    v6(:)     = tmp4(3,:)
    vz3(:)    = dx6(v6(:), dz(k))
    ! dwdx
    tmp4(:,:) = Q(4,i-2:i+2,j,k-2:k+3)
    wx3(:)    = dy56(tmp4(:,:), dx(i))
    ! dwdy & dwdz
    tmp4(:,:) = Q(4,i,j-2:j+2,k-2:k+3)
    wy3(:)    = dy56(tmp4(:,:), dy(j))
    w6(:)     = tmp4(3,:)
    wz3(:)    = dx6(w6(:), dz(k))
    ! dTdz
    T6(:)     = Q(5,i,j,k-2:k+3) / (R * Q(1,i,j,k-2:k+3))
    mu(:)     = mu6(T6(:))
    call tauxy4(mu(:), wx3(:), uz3(:), u6(:), tzx, utzx)
    call tauxy4(mu(:), vz3(:), wy3(:), v6(:), tzy, vtzy)
    call tauxx4(mu(:), wz3(:), ux3(:), vy3(:), w6(:), tzz, wtzz)
    kTz = heat_conduction6(mu(:), T6(:), dz(k))
    mutz   = 0.0625d0 * (-mut(i,j,k-1) + 9.d0 * (mut(i,j,k) + mut(i,j,k+1)) -mut(i,j,k+2))
    tzx  = tzx + mutz * (wx3(2) + uz3(2))
    tzy  = tzy + mutz * (vz3(2) + wy3(2))
    tzz  = tzz + 2.d0 * mutz * (2.d0 * wz3(2) - ux3(2) - vy3(2)) / 3.d0
    H(:) = (gamma * Q(5,i,j,k-1:k+2) / (Q(1,i,j,k-1:k+2) * gamma_1)) &
           + 0.5d0 * (Q(2,i,j,k-1:k+2)**2 + Q(3,i,j,k-1:k+2)**2 + Q(4,i,j,k-1:k+2)**2) + qc2(i,j,k-1:k+2)
    Hsgs = -mutz * 0.125d0 * (9.d0 * (-H(2) + H(3)) - (-H(1) + H(4)) / 3.d0) * dz(k) / Prt
  end subroutine calc_me4_base_les_z
end module calc_me4_base

