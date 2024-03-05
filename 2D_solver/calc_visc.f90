module calc_visc
  implicit none
contains
  attributes(global) subroutine calc_Ev2(id, nx, ny, dX, dY, mu, kappa, u, v, T, Ev)
    integer(kind=2), intent(in), value :: id
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: dX, dY, mu, kappa
    real(8), intent(in), dimension(nx,ny), device :: u, v, T
    real(8), intent(out), dimension(nx-1,ny-2,4), device :: Ev
    integer i, j
    real(8) ux, uy, vx, vy, txx, txy
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1

    ux = dX * (-u(i,j) + u(i+1,j))
    vx = dX * (-v(i,j) + v(i+1,j))
    uy = dY * 0.5d0 * (-u(i,j-1) + u(i,j+1))
    vy = dY * 0.5d0 * (-v(i,j-1) + v(i,j+1))
    txx = 2.d0 * mu * (2.d0 * ux - vy) / 3.d0
    txy = mu * (uy + vx)
    Ev(i,j-1,1) = 0.d0
    Ev(i,j-1,2) = txx
    Ev(i,j-1,3) = txy
    Ev(i,j-1,4) = txx * 0.5d0 *(u(i,j) + u(i+1,j)) + txy * 0.5d0 *(v(i,j) + v(i+1,j)) + kappa * (-T(i,j) + T(i+1,j))
  end subroutine calc_Ev2
  
  attributes(global) subroutine calc_Fv2(id, nx, ny, dX, dY, mu, kappa, u, v, T, Fv)
    integer(kind=2), intent(in), value :: id
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: dX, dY, mu, kappa
    real(8), intent(in), dimension(nx,ny), device :: u, v, T
    real(8), intent(out), device :: Fv(nx-2,ny-1,4)
    integer i, j
    real(8) ux, uy, vx, vy, tyx, tyy
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y

    uy = dY * (-u(i,j) + u(i,j+1))
    vy = dY * (-v(i,j) + v(i,j+1))
    ux = dX * 0.5d0 * (-u(i-1,j) + u(i+1,j))
    vx = dX * 0.5d0 * (-v(i-1,j) + v(i+1,j))
    tyx = mu * (uy + vx)
    tyy = 2.d0 * mu * (2.d0 * vy - ux) / 3.d0
    Fv(i-1,j,1) = 0.d0
    Fv(i-1,j,2) = tyx
    Fv(i-1,j,3) = tyy
    Fv(i-1,j,4) = tyx * 0.5d0 * (u(i,j) + u(i,j+1)) + tyy * 0.5d0 * (v(i,j) + v(i,j+1)) + kappa * (-T(i,j) + T(i,j+1))
  end subroutine calc_Fv2
end module calc_visc
