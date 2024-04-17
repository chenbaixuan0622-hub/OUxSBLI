module calc_hybrid
  use cudafor
  use mod_globals, only : accuracy, offset, nx, ny, dxi, dyi
  implicit none
contains
  attributes(global) subroutine calc_E_hybrid(rho,u,v,e,E_keep,E_upwind,E_hybrid)
    real(8), intent(in), dimension(nx,ny) :: rho, u, v, e
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device :: E_keep, E_upwind
    real(8), intent(out), dimension(nx-accuracy+1,ny-accuracy,4), device :: E_hybrid
    real(8) d1, d2, d3, phi_p, phi_m, phi, E_tvd(4), div, rot, fd
    real(8) :: eps = 1.d-16
    integer i, j
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset - 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    if (2 <= i) then
      d1 = -e(i-1,j) / rho(i-1,j) + e(i,j) / rho(i,j)
    else
      d1  = 0.d0
    endif
    d2 = -e(i,j) / rho(i,j) + e(i+1,j) / rho(i+1,j)
    if (i <= nx-2) then
      d3 = -e(i+1,j) / rho(i+1,j) + e(i+2,j) / rho(i+2,j)
    else
      d3 = 0.d0
    endif
    phi_p = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    phi_m = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    phi = min(phi_p, phi_m)
    !phi = 0.d0
    E_tvd(:) = phi * E_keep(i-offset+1,j-offset,:) + (1.d0 - phi) * E_upwind(i-offset+1,j-offset,:)
    
    div = dxi * (-u(i,j) + u(i+1,j)) + 0.5d0 * dyi * (-v(i,j-1) + v(i,j+1))
    rot = dxi * (-v(i,j) + v(i+1,j)) - 0.5d0 * dyi * (-u(i,j-1) + u(i,j+1))
    !fd = div**2 / (div**2 + rot**2 + eps)
    fd = 1.d0
    E_hybrid(i-offset+1,j-offset,:) = (1.d0 - fd) * E_keep(i-offset+1,j-offset,:) + fd * E_tvd(:)
  end subroutine calc_E_hybrid
  
  attributes(global) subroutine calc_F_hybrid(rho,u,v,e,F_keep,F_upwind,F_hybrid)
    real(8), intent(in), dimension(nx,ny) :: rho, u, v, e
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device :: F_keep, F_upwind
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy+1,4), device :: F_hybrid
    real(8) d1, d2, d3, phi_p, phi_m, phi, F_tvd(4), div, rot, fd
    real(8) :: eps = 1.d-16
    integer i, j
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset - 1
    if (2 <= j) then
      d1 = -e(i,j-1) / rho(i,j-1) + e(i,j) / rho(i,j)
    else
      d1  = 0.d0
    endif
    d2 = -e(i,j) / rho(i,j) + e(i,j+1) / rho(i,j+1)
    if (j <= ny-2) then
      d3 = -e(i,j+1) / rho(i,j+1) + e(i,j+2) / rho(i,j+2)
    else
      d3 = 0.d0
    endif
    phi_p = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    phi_m = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    phi = min(phi_p, phi_m)
    !phi = 0.d0
    F_tvd(:) = phi * F_keep(i-offset,j-offset+1,:) + (1.d0 - phi) * F_upwind(i-offset,j-offset+1,:)
    
    div = 0.5d0 * dxi * (-u(i-1,j) + u(i+1,j)) + dyi * (-v(i,j) + v(i,j+1))
    rot = 0.5d0 * dxi * (-v(i-1,j) + v(i+1,j)) - dyi * (-u(i,j) + u(i,j+1))
    !fd = div**2 / (div**2 + rot**2 + eps)
    fd = 1.d0
    F_hybrid(i-offset,j-offset+1,:) = (1.d0 - fd) * F_keep(i-offset,j-offset+1,:) + fd * F_tvd(:)
  end subroutine calc_F_hybrid
end module calc_hybrid

