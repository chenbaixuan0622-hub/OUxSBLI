module calc_hybrid
  use cudafor
  use mod_globals, only : dzi
  implicit none
contains
  attributes(global) subroutine calc_Ducros(nx,ny,nz,dx,dy,u,v,w,fd)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device      :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1), device      :: dy ! 1 / dy
    real(8), intent(in), dimension(nx,ny,nz), device  :: u, v, w
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
    dudz = (-u(i,j,k-1) + u(i,j,k+1)) * 0.5d0 * dzi
    dvdz = (-v(i,j,k-1) + v(i,j,k+1)) * 0.5d0 * dzi
    dwdz = (-w(i,j,k-1) + w(i,j,k+1)) * 0.5d0 * dzi
    div = dudx + dvdy + dwdz
    rot(1) = dwdy - dvdz
    rot(2) = dudz - dwdx
    rot(3) = dvdx - dudy
    fd(i,j,k) = (div**2) / (div**2 + (rot(1)**2 + rot(2)**2 + rot(3)**2) + eps)

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

  attributes(global) subroutine calc_E_hybrid(nx,ny,nz,u,v,w,fd,E_upwind,E)
    integer, intent(in), value                                  :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device            :: u, v, w, fd
    real(8), intent(in), dimension(nx-1,ny-2,nz-2,5), device    :: E_upwind
    real(8), intent(inout), dimension(nx-1,ny-2,nz-2,5), device :: E
    ! Albada
    !real(8) d1, d2, d3, phi_p, phi_m, phi, E_tvd(5)
    integer i, j, k
    real(8) fdx
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    !if (2 <= i) then
    !  d1 = -e(i-1,j,k) / rho(i-1,j,k) + e(i,j,k) / rho(i,j,k)
    !else
    !  d1  = 0.d0
    !endif
    !d2 = -e(i,j,k) / rho(i,j,k) + e(i+1,j,k) / rho(i+1,j,k)
    !if (i <= nx-2) then
    !  d3 = -e(i+1,j,k) / rho(i+1,j,k) + e(i+2,j,k) / rho(i+2,j,k)
    !else
    !  d3 = 0.d0
    !endif
    !phi_p = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    !phi_m = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    !phi = min(phi_p, phi_m)
    !E_tvd(:) = phi * E_keep(i-1+1,j-1,k-1,:) + (1.d0 - phi) * E_upwind(i-1+1,j-1,k-1,:)

    fdx = max(fd(i,j,k), fd(i+1,j,k))
    if (0.4d0 <= fdx) then
      fdx = 1.d0
    else
      fdx = 0.d0
    endif
    !fdx = min(1.d0, max(0.d0, fdx))
    E(i,j-1,k-1,:) = (1.d0 - fdx) * E(i,j-1,k-1,:) + fdx * E_upwind(i,j-1,k-1,:)
  end subroutine calc_E_hybrid
  
  attributes(global) subroutine calc_F_hybrid(nx,ny,nz,u,v,w,fd,F_upwind,F)
    integer, intent(in), value                                  :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device            :: u, v, w, fd
    real(8), intent(in), dimension(nx-2,ny-1,nz-2,5), device    :: F_upwind
    real(8), intent(inout), dimension(nx-2,ny-1,nz-2,5), device :: F
    ! Albada
    !real(8) d1, d2, d3, phi_p, phi_m, phi, F_tvd(5)
    integer i, j, k
    real(8) fdy
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    !if (2 <= j) then
    !  d1 = -e(i,j-1,k) / rho(i,j-1,k) + e(i,j,k) / rho(i,j,k)
    !else
    !  d1  = 0.d0
    !endif
    !d2 = -e(i,j,k) / rho(i,j,k) + e(i,j+1,k) / rho(i,j+1,k)
    !if (j <= ny-2) then
    !  d3 = -e(i,j+1,k) / rho(i,j+1,k) + e(i,j+2,k) / rho(i,j+2,k)
    !else
    !  d3 = 0.d0
    !endif
    !phi_p = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    !phi_m = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    !phi = min(phi_p, phi_m)
    !F_tvd(:) = phi * F_keep(i-1,j-1+1,k-1,:) + (1.d0 - phi) * F_upwind(i-1,j-1+1,k-1,:)

    fdy = max(fd(i,j,k), fd(i,j+1,k))
    if (0.4d0 <= fdy) then
      fdy = 1.d0
    else
      fdy = 0.d0
    endif
    !fdy = min(1.d0, max(fdy, 0.d0))
    F(i-1,j,k-1,:) = (1.d0 - fdy) * F(i-1,j,k-1,:) + fdy * F_upwind(i-1,j,k-1,:)
  end subroutine calc_F_hybrid

  attributes(global) subroutine calc_G_hybrid(nx,ny,nz,u,v,w,fd,G_upwind,G)
    integer, intent(in), value                                  :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device            :: u, v, w, fd
    real(8), intent(in), dimension(nx-2,ny-2,nz-1,5), device    :: G_upwind
    real(8), intent(inout), dimension(nx-2,ny-2,nz-1,5), device :: G
    ! Albada
    !real(8) d1, d2, d3, phi_p, phi_m, phi, G_tvd(5)
    integer i, j, k
    real(8) fdz
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    !if (2 <= k) then
    !  d1 = -e(i,j,k-1) / rho(i,j,k-1) + e(i,j,k) / rho(i,j,k)
    !else
    !  d1  = 0.d0
    !endif
    !d2 = -e(i,j,k) / rho(i,j,k) + e(i,j,k+1) / rho(i,j,k+1)
    !if (k <= nz-2) then
    !  d3 = -e(i,j,k+1) / rho(i,j,k+1) + e(i,j,k+2) / rho(i,j,k+2)
    !else
    !  d3 = 0.d0
    !endif
    !phi_p = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    !phi_m = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    !phi = min(phi_p, phi_m)
    !G_tvd(:) = phi * G_keep(i-1,j-1,k-1+1,:) + (1.d0 - phi) * G_upwind(i-1,j-1,k-1+1,:)
    
    fdz = max(fd(i,j,k), fd(i,j,k+1))
    if (0.4d0 <= fdz) then
      fdz = 1.d0
    else
      fdz = 0.d0
    endif
    !fdz = min(1.d0, max(0.d0, fdz))
    G(i-1,j-1,k,:) = (1.d0 - fdz) * G(i-1,j-1,k,:) + fdz * G_upwind(i-1,j-1,k,:)
  end subroutine calc_G_hybrid
end module calc_hybrid

