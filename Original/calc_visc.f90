module calc_visc
  use calc_term
  implicit none
contains
  subroutine calc_Ev(id, nx, ny, nz, dX, dY, dZ, mu, kappa, u, v, w, T, Ev)
    integer(kind=2), intent(in) :: id
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in) :: dX, dY, dZ, mu, kappa
    real(8), intent(in), dimension(nx,ny,nz) :: u, v, w, T
    real(8), intent(out), dimension(nx-1,ny-2,nz-2,5) :: Ev
    integer i, j, k
    real(8) ux, uy, uz, vx, vy, wx, wz, txx, txy, txz
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 1, nx-1
          ux = dX * (-u(i,j,k) + u(i+1,j,k))
          vx = dX * (-v(i,j,k) + v(i+1,j,k))
          wx = dX * (-w(i,j,k) + w(i+1,j,k))
          uy = dY * 0.5d0 * (-u(i,j-1,k) + u(i,j+1,k))
          vy = dY * 0.5d0 * (-v(i,j-1,k) + v(i,j+1,k))
          uz = dZ * 0.5d0 * (-u(i,j,k-1) + u(i,j,k+1))
          wz = dZ * 0.5d0 * (-w(i,j,k-1) + w(i,j,k+1))
          txx = 2.d0 * mu * (2.d0 * ux - vy - wz) / 3.d0
          txy = mu * (uy + vx)
          txz = mu * (wx + uz)
          Ev(i,j-1,k-1,1) = 0.d0
          Ev(i,j-1,k-1,2) = txx
          Ev(i,j-1,k-1,3) = txy
          Ev(i,j-1,k-1,4) = txz
          Ev(i,j-1,k-1,5) = txx * 0.5d0 *(u(i,j,k) + u(i+1,j,k)) + txy * 0.5d0 *(v(i,j,k) + v(i+1,j,k)) &
          & + txz * 0.5d0 *(w(i,j,k) + w(i+1,j,k)) + kappa * (-T(i,j,k) + T(i+1,j,k))
        enddo
      enddo
    enddo
  end subroutine calc_Ev
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  subroutine calc_Fv(id, nx, ny, nz, dX, dY, dZ, mu, kappa, u, v, w, T, Fv)
    integer(kind=2), intent(in) :: id
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in) :: dX, dY, dZ, mu, kappa
    real(8), intent(in), dimension(nx,ny,nz) :: u, v, w, T
    real(8), intent(out) :: Fv(nx-2,ny-1,nz-2,5)
    integer i, j, k
    real(8) ux, uy, vx, vy, vz, wy, wz, tyx, tyy, tyz
    do k = 2, nz-1
      do j = 1, ny-1
        do i = 2, nx-1
          uy = dY * (-u(i,j,k) + u(i,j+1,k))
          vy = dY * (-v(i,j,k) + v(i,j+1,k))
          wy = dY * (-w(i,j,k) + w(i,j+1,k))
          ux = dX * 0.5d0 * (-u(i-1,j,k) + u(i+1,j,k))
          vx = dX * 0.5d0 * (-v(i-1,j,k) + v(i+1,j,k))
          vz = dZ * 0.5d0 * (-v(i,j,k-1) + v(i,j,k+1))
          wz = dZ * 0.5d0 * (-w(i,j,k-1) + w(i,j,k+1))
          tyx = mu * (uy + vx)
          tyy = 2.d0 * mu * (2.d0 * vy - wz - ux) / 3.d0
          tyz = mu * (vz + wy)
          Fv(i-1,j,k-1,1) = 0.d0
          Fv(i-1,j,k-1,2) = tyx
          Fv(i-1,j,k-1,3) = tyy
          Fv(i-1,j,k-1,4) = tyz
          Fv(i-1,j,k-1,5) = tyx * 0.5d0 * (u(i,j,k) + u(i,j+1,k)) + tyy * 0.5d0 * (v(i,j,k) + v(i,j+1,k)) &
          & + tyz * 0.5d0 * (w(i,j,k) + w(i,j+1,k)) + kappa * (-T(i,j,k) + T(i,j+1,k))
        enddo
      enddo
    enddo
  end subroutine calc_Fv
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  subroutine calc_Gv(id, nx, ny, nz, dX, dY, dZ, mu, kappa, u, v, w, T, Gv)
    integer(kind=2), intent(in) :: id
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in) :: dX, dY, dZ, mu, kappa
    real(8), intent(in), dimension(nx,ny,nz) :: u, v, w, T
    real(8), intent(out) :: Gv(nx-2,ny-2,nz-1,5)
    integer i, j, k
    real(8) ux, uz, vy, vz, wx, wy, wz, tzx, tzy, tzz
    do k = 1, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          uz = dZ * (-u(i,j,k) + u(i,j,k+1))
          vz = dZ * (-v(i,j,k) + v(i,j,k+1))
          wz = dZ * (-w(i,j,k) + w(i,j,k+1))
          ux = dX * 0.5d0 * (-u(i-1,j,k) + u(i+1,j,k))
          wx = dX * 0.5d0 * (-w(i-1,j,k) + w(i+1,j,k))
          vy = dY * 0.5d0 * (-v(i,j-1,k) + v(i,j+1,k))
          wy = dY * 0.5d0 * (-w(i,j-1,k) + w(i,j+1,k))
          tzx = mu * (wx + uz)
          tzy = mu * (vz + wy)
          tzz = 2.d0 * mu * (2.d0 * wz - ux - vy) / 3.d0
          Gv(i-1,j-1,k,1) = 0.d0
          Gv(i-1,j-1,k,2) = tzx
          Gv(i-1,j-1,k,3) = tzy
          Gv(i-1,j-1,k,4) = tzz
          Gv(i-1,j-1,k,5) = tzx * 0.5d0 * (u(i,j,k) + u(i,j,k+1)) + tzy * 0.5d0 * (v(i,j,k) + v(i,j,k+1)) &
          & + tzz * 0.5d0 * (w(i,j,k) + w(i,j,k+1)) + kappa * (-T(i,j,k) + T(i,j,k+1))
        enddo
      enddo
    enddo
  end subroutine calc_Gv
end module calc_visc
  