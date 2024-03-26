module calc_hybrid
  use cudafor
  use mod_globals, only : accuracy, offset, nx, ny, nz, dxi, dyi, dzi
  implicit none
contains
  attributes(global) subroutine calc_E_tvd(Q,E_keep,E_upwind,E_tvd)
    real(8), intent(in), dimension(nx,ny,nz,5) :: Q
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E_keep, E_upwind
    real(8), intent(out), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E_tvd
    real(8) d1, d2, d3, phi_p, phi_m, phi
    real(8) :: eps = 1.d-16
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    if (2 <= i) then
      d1 = -Q(i-1,j,k,5) / Q(i-1,j,k,1) + Q(i,j,k,5) / Q(i,j,k,1)
    else
      d1  = 0.d0
    endif
    d2 = -Q(i,j,k,5) / Q(i,j,k,1) + Q(i+1,j,k,5) / Q(i+1,j,k,1)
    if (i <= nx-2) then
      d3 = -Q(i+1,j,k,5) / Q(i+1,j,k,1) + Q(i+2,j,k,5) / Q(i+2,j,k,1)
    else
      d3 = 0.d0
    endif
    phi_p = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    phi_m = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    phi = min(phi_p, phi_m)
    E_tvd(i,j-offset,k-offset,:) = phi * E_keep(i,j-offset,k-offset,:) + (1.d0 - phi) * E_upwind(i,j-offset,k-offset,:)
  end subroutine calc_E_tvd
  
  attributes(global) subroutine calc_F_tvd(Q,F_keep,F_upwind,F_tvd)
    real(8), intent(in), dimension(nx,ny,nz,5) :: Q
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F_keep, F_upwind
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F_tvd
    real(8) d1, d2, d3, phi_p, phi_m, phi
    real(8) :: eps = 1.d-16
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    if (2 <= j) then
      d1 = -Q(i,j-1,k,5) / Q(i,j-1,k,1) + Q(i,j,k,5) / Q(i,j,k,1)
    else
      d1  = 0.d0
    endif
    d2 = -Q(i,j,k,5) / Q(i,j,k,1) + Q(i,j+1,k,5) / Q(i,j+1,k,1)
    if (j <= ny-2) then
      d3 = -Q(i,j+1,k,5) / Q(i,j+1,k,1) + Q(i,j+2,k,5) / Q(i,j+2,k,1)
    else
      d3 = 0.d0
    endif
    phi_p = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    phi_m = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    phi = min(phi_p, phi_m)
    F_tvd(i-offset,j,k-offset,:) = phi * F_keep(i-offset,j,k-offset,:) + (1.d0 - phi) * F_upwind(i-offset,j,k-offset,:)
  end subroutine calc_F_tvd

  attributes(global) subroutine calc_G_tvd(Q,G_keep,G_upwind,G_tvd)
    real(8), intent(in), dimension(nx,ny,nz,5) :: Q
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G_keep, G_upwind
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G_tvd
    real(8) d1, d2, d3, phi_p, phi_m, phi
    real(8) :: eps = 1.d-16
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    if (2 <= k) then
      d1 = -Q(i,j,k-1,5) / Q(i,j,k-1,1) + Q(i,j,k,5) / Q(i,j,k,1)
    else
      d1  = 0.d0
    endif
    d2 = -Q(i,j,k,5) / Q(i,j,k,1) + Q(i,j,k+1,5) / Q(i,j,k+1,1)
    if (k <= nz-2) then
      d3 = -Q(i,j,k+1,5) / Q(i,j,k+1,1) + Q(i,j,k+2,5) / Q(i,j,k+2,1)
    else
      d3 = 0.d0
    endif
    phi_p = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    phi_m = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    phi = min(phi_p, phi_m)
    G_tvd(i-offset,j-offset,k,:) = phi * G_keep(i-offset,j-offset,k,:) + (1.d0 - phi) * G_upwind(i-offset,j-offset,k,:)
  end subroutine calc_G_tvd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(global) subroutine calc_E_hybrid(u,v,w,E_keep,E_tvd,E_hybrid)
    real(8), intent(in), dimension(nx,ny,nz), device :: u, v, w
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E_keep, E_tvd
    real(8), intent(out), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E_hybrid
    real(8) div, rot(3), fd
    real(8) :: eps = 1.d-6
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    div = dxi * (-u(i,j,k) + u(i+1,j,k)) + 0.5d0 * dyi * (-v(i,j-1,k) + v(i,j+1,k)) + 0.5d0 * dzi * (-w(i,j,k-1) + w(i,j,k+1))
    
    rot(1) = 0.5d0 * dyi * (-w(i,j-1,k) + w(i,j+1,k)) - 0.5d0 * dzi * (-v(i,j,k-1) + v(i,j,k+1))
    rot(2) = 0.5d0 * dzi * (-u(i,j,k-1) + u(i,j,k+1)) - dxi * (-w(i,j,k) + w(i+1,j,k))
    rot(3) = dxi * (-v(i,j,k) + v(i+1,j,k)) - 0.5d0 * dyi * (-u(i,j-1,k) + u(i,j+1,k))
    
    fd = div**2 / (div**2 + (rot(1)**2 + rot(2)**2 + rot(3)**2) + eps)
    E_hybrid(i-offset,j-offset,k-offset,:) = (1.d0 - fd) * E_keep(i-offset,j-offset,k-offset,:) + fd * E_tvd(i-offset,j-offset,k-offset,:)
  end subroutine calc_E_hybrid
  
  attributes(global) subroutine calc_F_hybrid(u,v,w,F_keep,F_tvd,F_hybrid)
    real(8), intent(in), dimension(nx,ny,nz), device :: u, v, w
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F_keep, F_tvd
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F_hybrid
    real(8) div, rot(3), fd
    real(8) :: eps = 1.d-6
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    div = 0.5d0 * dxi * (-u(i-1,j,k) + u(i+1,j,k)) + dyi * (-v(i,j,k) + v(i,j+1,k)) + 0.5d0 * dzi * (-w(i,j,k-1) + w(i,j,k+1))
    
    rot(1) = dyi * (-w(i,j,k) + w(i,j+1,k)) - 0.5d0 * dzi * (-v(i,j,k-1) + v(i,j,k+1))
    rot(2) = 0.5d0 * dzi * (-u(i,j,k-1) + u(i,j,k+1)) - 0.5d0 * dxi * (-w(i-1,j,k) + w(i+1,j,k))
    rot(3) = 0.5d0 * dxi * (-v(i-1,j,k) + v(i+1,j,k)) - dyi * (-u(i,j,k) + u(i,j+1,k))
    
    fd = div**2 / (div**2 + (rot(1)**2 + rot(2)**2 + rot(3)**2) + eps)
    F_hybrid(i-offset,j-offset,k-offset,:) = (1.d0 - fd) * F_keep(i-offset,j-offset,k-offset,:) + fd * F_tvd(i-offset,j-offset,k-offset,:)
  end subroutine calc_F_hybrid
  
  attributes(global) subroutine calc_G_hybrid(u,v,w,G_keep,G_tvd,G_hybrid)
    real(8), intent(in), dimension(nx,ny,nz), device :: u, v, w
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G_keep, G_tvd
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G_hybrid
    real(8) div, rot(3), fd
    real(8) :: eps = 1.d-6
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    div = 0.5d0 * dxi * (-u(i-1,j,k) + u(i+1,j,k)) + 0.5d0 * dyi * (-v(i,j-1,k) + v(i,j+1,k)) + dzi * (-w(i,j,k) + w(i,j,k+1))
    
    rot(1) = 0.5d0 * dyi * (-w(i,j-1,k) + w(i,j+1,k)) - dzi * (-v(i,j,k) + v(i,j,k+1))
    rot(2) = dzi * (-u(i,j,k) + u(i,j,k+1)) - 0.5d0 * dxi * (-w(i-1,j,k) + w(i+1,j,k))
    rot(3) = 0.5d0 * dxi * (-v(i-1,j,k) + v(i+1,j,k)) - 0.5d0 * dyi * (-u(i,j-1,k) + u(i,j+1,k))
    
    fd = div**2 / (div**2 + (rot(1)**2 + rot(2)**2 + rot(3)**2) + eps)
    G_hybrid(i-offset,j-offset,k-offset,:) = (1.d0 - fd) * G_keep(i-offset,j-offset,k-offset,:) + fd * G_tvd(i-offset,j-offset,k-offset,:)
  end subroutine calc_G_hybrid
end module calc_hybrid

