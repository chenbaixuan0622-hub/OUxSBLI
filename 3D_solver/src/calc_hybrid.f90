module calc_hybrid
  use cudafor
  use mod_globals, only : accuracy, offset, gamma
  implicit none
contains
  attributes(global) subroutine calc_Ducros(nx,ny,nz,dx,dy,dz,u,v,w,rho,p,fd)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device      :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1), device      :: dy ! 1 / dy
    real(8), intent(in), dimension(nz-1), device      :: dz ! 1 / dz
    real(8), intent(in), dimension(nx,ny,nz), device  :: u, v, w, rho, p
    real(8), intent(out), dimension(nx,ny,nz), device :: fd
    integer i, j, k
    real(8) dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    real(8) div, rot(3)
    real(8) :: eps = 1.d-16
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1 
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    dudx = (-u(i-1,j,k) + u(i+1,j,k)) * 0.25d0 * (dx(i-1) + dx(i))
    dvdx = (-v(i-1,j,k) + v(i+1,j,k)) * 0.25d0 * (dx(i-1) + dx(i))
    dwdx = (-w(i-1,j,k) + w(i+1,j,k)) * 0.25d0 * (dx(i-1) + dx(i))
    dudy = (-u(i,j-1,k) + u(i,j+1,k)) * 0.25d0 * (dy(j-1) + dy(j))
    dvdy = (-v(i,j-1,k) + v(i,j+1,k)) * 0.25d0 * (dy(j-1) + dy(j))
    dwdy = (-w(i,j-1,k) + w(i,j+1,k)) * 0.25d0 * (dy(j-1) + dy(j))
    dudz = (-u(i,j,k-1) + u(i,j,k+1)) * 0.25d0 * (dz(k-1) + dz(k))
    dvdz = (-v(i,j,k-1) + v(i,j,k+1)) * 0.25d0 * (dz(k-1) + dz(k))
    dwdz = (-w(i,j,k-1) + w(i,j,k+1)) * 0.25d0 * (dz(k-1) + dz(k))
    div = dudx + dvdy + dwdz
    rot(1) = dwdy - dvdz
    rot(2) = dudz - dwdx
    rot(3) = dvdx - dudy
    fd(i,j,k) = (div**2) / (div**2 + (rot(1)**2 + rot(2)**2 + rot(3)**2) + eps)

    fd(i,j,k) = min(1.d0, fd(i,j,k))

    ! boundary
    ! x direction
    if (i == 2) then
      fd(1,j,k) = fd(2,j,k)
    elseif (i == nx-1) then
      fd(nx,j,k) = fd(nx-1,j,k)
    endif
    ! y direction
    if (j == 2) then
      fd(i,1,k) = fd(i,2,k)
    elseif (j == ny-1) then
      fd(i,ny,k) = fd(i,ny-1,k)
    endif
    ! z direction
    if (k == 2) then
      fd(i,j,1) = fd(i,j,2)
    elseif (k == nz-1) then
      fd(i,j,nz) = fd(i,j,nz-1)
    endif
  end subroutine calc_Ducros

  attributes(device) function Albada(rho,p,V) result(phi)
    real(8), intent(in), dimension(4), device   :: rho, p
    real(8), intent(in), dimension(4,3), device :: V
    integer i
    real(8) :: e(4), d1, d2, d3, phim, phip, phi, eps = 1.d-16
    do i = 1, 4
      e(i)  = p(i) / (gamma - 1.d0) + 0.5d0 * rho(i) * (V(i,1)**2 + V(i,2)**2 + V(i,3)**2)
    enddo
    d1   = -e(1) / rho(1) + e(2) / rho(2)
    d2   = -e(2) / rho(2) + e(3) / rho(3)
    d3   = -e(3) / rho(3) + e(4) / rho(4)
    phip = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    phim = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    phi  = min(phim, phip)
  end function Albada

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

