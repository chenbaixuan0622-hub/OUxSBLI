module calc_hybrid
  use cudafor
  use mod_globals, only : accuracy, offset, nx, ny, nz, dxi, dyi, dzi
  implicit none
contains
  attributes(global) subroutine calc_E_hybrid(rho,u,v,w,e,E_keep,E_upwind,E_hybrid)
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w, e
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E_keep, E_upwind
    real(8), intent(out), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E_hybrid
    real(8) d1, d2, d3, phi_p, phi_m, phi, E_tvd(5), div, rot(3), fd
    real(8) :: eps = 1.d-16
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset - 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    if (2 <= i) then
      d1 = -e(i-1,j,k) / rho(i-1,j,k) + e(i,j,k) / rho(i,j,k)
    else
      d1  = 0.d0
    endif
    d2 = -e(i,j,k) / rho(i,j,k) + e(i+1,j,k) / rho(i+1,j,k)
    if (i <= nx-2) then
      d3 = -e(i+1,j,k) / rho(i+1,j,k) + e(i+2,j,k) / rho(i+2,j,k)
    else
      d3 = 0.d0
    endif
    phi_p = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    phi_m = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    phi = min(phi_p, phi_m)
    E_tvd(:) = phi * E_keep(i-offset+1,j-offset,k-offset,:) + (1.d0 - phi) * E_upwind(i-offset+1,j-offset,k-offset,:)
    
    div = dxi * (-u(i,j,k) + u(i+1,j,k)) + 0.5d0 * dyi * (-v(i,j-1,k) + v(i,j+1,k)) + 0.5d0 * dzi * (-w(i,j,k-1) + w(i,j,k+1))
    
    rot(1) = 0.5d0 * dyi * (-w(i,j-1,k) + w(i,j+1,k)) - 0.5d0 * dzi * (-v(i,j,k-1) + v(i,j,k+1))
    rot(2) = 0.5d0 * dzi * (-u(i,j,k-1) + u(i,j,k+1)) - dxi * (-w(i,j,k) + w(i+1,j,k))
    rot(3) = dxi * (-v(i,j,k) + v(i+1,j,k)) - 0.5d0 * dyi * (-u(i,j-1,k) + u(i,j+1,k))
    
    fd = div**2 / (div**2 + (rot(1)**2 + rot(2)**2 + rot(3)**2) + eps)
    E_hybrid(i-offset+1,j-offset,k-offset,:) = (1.d0 - fd) * E_keep(i-offset+1,j-offset,k-offset,:) + fd * E_tvd(:)
  end subroutine calc_E_hybrid
  
  attributes(global) subroutine calc_F_hybrid(rho,u,v,w,e,F_keep,F_upwind,F_hybrid)
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w, e
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F_keep, F_upwind
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F_hybrid
    real(8) d1, d2, d3, phi_p, phi_m, phi, F_tvd(5), div, rot(3), fd
    real(8) :: eps = 1.d-16
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset - 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    if (2 <= j) then
      d1 = -e(i,j-1,k) / rho(i,j-1,k) + e(i,j,k) / rho(i,j,k)
    else
      d1  = 0.d0
    endif
    d2 = -e(i,j,k) / rho(i,j,k) + e(i,j+1,k) / rho(i,j+1,k)
    if (j <= ny-2) then
      d3 = -e(i,j+1,k) / rho(i,j+1,k) + e(i,j+2,k) / rho(i,j+2,k)
    else
      d3 = 0.d0
    endif
    phi_p = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    phi_m = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    phi = min(phi_p, phi_m)
    F_tvd(:) = phi * F_keep(i-offset,j-offset+1,k-offset,:) + (1.d0 - phi) * F_upwind(i-offset,j-offset+1,k-offset,:)
    
    div = 0.5d0 * dxi * (-u(i-1,j,k) + u(i+1,j,k)) + dyi * (-v(i,j,k) + v(i,j+1,k)) + 0.5d0 * dzi * (-w(i,j,k-1) + w(i,j,k+1))
    
    rot(1) = dyi * (-w(i,j,k) + w(i,j+1,k)) - 0.5d0 * dzi * (-v(i,j,k-1) + v(i,j,k+1))
    rot(2) = 0.5d0 * dzi * (-u(i,j,k-1) + u(i,j,k+1)) - 0.5d0 * dxi * (-w(i-1,j,k) + w(i+1,j,k))
    rot(3) = 0.5d0 * dxi * (-v(i-1,j,k) + v(i+1,j,k)) - dyi * (-u(i,j,k) + u(i,j+1,k))
    
    fd = div**2 / (div**2 + (rot(1)**2 + rot(2)**2 + rot(3)**2) + eps)
    F_hybrid(i-offset,j-offset+1,k-offset,:) = (1.d0 - fd) * F_keep(i-offset,j-offset+1,k-offset,:) + fd * F_tvd(:)
  end subroutine calc_F_hybrid

  attributes(global) subroutine calc_G_hybrid(rho,u,v,w,e,G_keep,G_upwind,G_hybrid)
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w, e
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G_keep, G_upwind
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G_hybrid
    real(8) d1, d2, d3, phi_p, phi_m, phi, G_tvd(5), div, rot(3), fd
    real(8) :: eps = 1.d-16
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset - 1
    if (2 <= k) then
      d1 = -e(i,j,k-1) / rho(i,j,k-1) + e(i,j,k) / rho(i,j,k)
    else
      d1  = 0.d0
    endif
    d2 = -e(i,j,k) / rho(i,j,k) + e(i,j,k+1) / rho(i,j,k+1)
    if (k <= nz-2) then
      d3 = -e(i,j,k+1) / rho(i,j,k+1) + e(i,j,k+2) / rho(i,j,k+2)
    else
      d3 = 0.d0
    endif
    phi_p = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    phi_m = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    phi = min(phi_p, phi_m)
    G_tvd(:) = phi * G_keep(i-offset,j-offset,k-offset+1,:) + (1.d0 - phi) * G_upwind(i-offset,j-offset,k-offset+1,:)
    
    div = 0.5d0 * dxi * (-u(i-1,j,k) + u(i+1,j,k)) + 0.5d0 * dyi * (-v(i,j-1,k) + v(i,j+1,k)) + dzi * (-w(i,j,k) + w(i,j,k+1))
    
    rot(1) = 0.5d0 * dyi * (-w(i,j-1,k) + w(i,j+1,k)) - dzi * (-v(i,j,k) + v(i,j,k+1))
    rot(2) = dzi * (-u(i,j,k) + u(i,j,k+1)) - 0.5d0 * dxi * (-w(i-1,j,k) + w(i+1,j,k))
    rot(3) = 0.5d0 * dxi * (-v(i-1,j,k) + v(i+1,j,k)) - 0.5d0 * dyi * (-u(i,j-1,k) + u(i,j+1,k))
    
    fd = div**2 / (div**2 + (rot(1)**2 + rot(2)**2 + rot(3)**2) + eps)
    G_hybrid(i-offset,j-offset,k-offset+1,:) = (1.d0 - fd) * G_keep(i-offset,j-offset,k-offset+1,:) + fd * G_tvd(:)
  end subroutine calc_G_hybrid
end module calc_hybrid

