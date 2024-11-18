module calc_hybrid
  use cudafor
  use mod_globals, only : accuracy, offset, nx, ny, dxi, dyi
  implicit none
contains
  attributes(global) subroutine calc_Ducros(u,v,fd)
    real(8), intent(in), dimension(nx,ny), device :: u, v
    real(8), intent(out), dimension(nx,ny), device :: fd
    real(8) dudx, dudy, dvdx, dvdy
    real(8) div, rot
    real(8) :: eps = 1.d-16
    integer i, j
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    dudx = 0.5d0 * (-u(i-1,j) + u(i+1,j)) * dxi
    dvdx = 0.5d0 * (-v(i-1,j) + v(i+1,j)) * dxi
    dudy = 0.5d0 * (-u(i,j-1) + u(i,j+1)) * dyi
    dvdy = 0.5d0 * (-v(i,j-1) + v(i,j+1)) * dyi
    div = dudx + dvdy
    rot = dvdx - dudy
    fd(i,j) = (div**2) / (div**2 + rot**2 + eps)
    ! boundary
    ! x direction
    if (i == 2) then
      fd(1,j) = fd(2,j)
    elseif (i == nx-1) then
      fd(nx,j) = fd(nx-1,j)
    endif
    ! y direction
    if (j == 2) then
      fd(i,1) = fd(i,2)
    elseif (j == ny-1) then
      fd(i,ny) = fd(i,ny-1)
    endif
  end subroutine calc_Ducros

  attributes(device) function sigmoid(x) result(ans)
    real(8), intent(in), value :: x
    real(8) :: ans
    ans = 0.5d0 * (tanh(10.d0 * (x - 0.5d0)) + 1.d0)
  end function sigmoid

  attributes(device) function wiggle_detector(phi) result(ans)
    real(8), intent(in), device :: phi(4)
    real(8) ans, phi1, phi2
    phi1 = (-phi(1) + phi(2)) * (-phi(2) + phi(3))
    phi2 = (-phi(3) + phi(4)) * (-phi(2) + phi(3))
    ans  = 0.5d0 * (1.d0 - sign(1.d0, min(phi1, phi2)))
  end function wiggle_detector
end module calc_hybrid

