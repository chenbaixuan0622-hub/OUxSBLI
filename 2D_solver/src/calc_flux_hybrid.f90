module calc_flux_hybrid
  use cudafor
  use mod_globals, only : accuracy, offset, nx, ny, dxi, dyi
  implicit none
contains
  attributes(global) subroutine calc_E_tvd(Q,E_keep,E_upwind,E_tvd)
    real(8), intent(in), dimension(nx,ny,4) :: Q
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device :: E_keep, E_upwind
    real(8), intent(out), dimension(nx-accuracy+1,ny-accuracy,4), device :: E_tvd
    real(8) d1, d2, d3, phi_p, phi_m, phi
    real(8) :: eps = 1.d-16
    integer i, j
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    if (2 <= i) then
      d1 = -Q(i-1,j,4) / Q(i-1,j,1) + Q(i,j,4) / Q(i,j,1)
    else
      d1  = 0.d0
    endif
    d2 = -Q(i,j,4) / Q(i,j,1) + Q(i+1,j,4) / Q(i+1,j,1)
    if (i <= nx-2) then
      d3 = -Q(i+1,j,4) / Q(i+1,j,1) + Q(i+2,j,4) / Q(i+2,j,1)
    else
      d3 = 0.d0
    endif
    phi_p = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    phi_m = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    phi = min(phi_p, phi_m)
    E_tvd(i,j-offset,:) = phi * E_keep(i,j-offset,:) + (1.d0 - phi) * E_upwind(i,j-offset,:)
  end subroutine calc_E_tvd
  
  attributes(global) subroutine calc_F_tvd(Q,F_keep,F_upwind,F_tvd)
    real(8), intent(in), dimension(nx,ny,4) :: Q
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device :: F_keep, F_upwind
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy+1,4), device :: F_tvd
    real(8) d1, d2, d3, phi_p, phi_m, phi
    real(8) :: eps = 1.d-16
    integer i, j
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    if (2 <= j) then
      d1 = -Q(i,j-1,4) / Q(i,j-1,1) + Q(i,j,4) / Q(i,j,1)
    else
      d1  = 0.d0
    endif
    d2 = -Q(i,j,4) / Q(i,j,1) + Q(i,j+1,4) / Q(i,j+1,1)
    if (j <= ny-2) then
      d3 = -Q(i,j+1,4) / Q(i,j+1,1) + Q(i,j+2,4) / Q(i,j+2,1)
    else
      d3 = 0.d0
    endif
    phi_p = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    phi_m = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    phi = min(phi_p, phi_m)
    F_tvd(i-offset,j,:) = phi * F_keep(i-offset,j,:) + (1.d0 - phi) * F_upwind(i-offset,j,:)
  end subroutine calc_F_tvd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(global) subroutine calc_E_hybrid(u,v,E_keep,E_tvd,E_hybrid)
    real(8), intent(in), dimension(nx,ny), device :: u, v
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device :: E_keep, E_tvd
    real(8), intent(out), dimension(nx-accuracy+1,ny-accuracy,4), device :: E_hybrid
    real(8) div, rot, fd
    real(8) :: eps = 1.d-6
    integer i, j
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    div = dxi * (-u(i,j) + u(i+1,j)) + 0.5d0 * dyi * (-v(i,j-1) + v(i,j+1))
    rot = dxi * (-v(i,j) + v(i+1,j)) - 0.5d0 * dyi * (-u(i,j-1) + u(i,j+1))
    fd = div**2 / (div**2 + rot**2 + eps)
    E_hybrid(i-offset,j-offset,:) = (1.d0 - fd) * E_keep(i-offset,j-offset,:) + fd * E_tvd(i-offset,j-offset,:)
  end subroutine calc_E_hybrid
  
  attributes(global) subroutine calc_F_hybrid(u,v,F_keep,F_tvd,F_hybrid)
    real(8), intent(in), dimension(nx,ny), device :: u, v
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device :: F_keep, F_tvd
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy+1,4), device :: F_hybrid
    real(8) div, rot, fd
    real(8) :: eps = 1.d-6
    integer i, j
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    div = 0.5d0 * dxi * (-u(i-1,j) + u(i+1,j)) + dyi * (-v(i,j) + v(i,j+1))
    rot = 0.5d0 * dxi * (-v(i-1,j) + v(i+1,j)) - dyi * (-u(i,j) + u(i,j+1))
    fd = div**2 / (div**2 + rot**2 + eps)
    F_hybrid(i-offset,j-offset,:) = (1.d0 - fd) * F_keep(i-offset,j-offset,:) + fd * F_tvd(i-offset,j-offset,:)
  end subroutine calc_F_hybrid
end module calc_flux_hybrid

