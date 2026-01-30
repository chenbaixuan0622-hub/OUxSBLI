module calc_igr
  use cusparse
  use mod_globals, only : blocks, threads, dx, over_dx, alpha
  implicit none
contains
  attributes(global) subroutine set_band(nx, Q, dl, d, du, b)
    integer, intent(in), value   :: nx
    real(8), intent(in), device  :: Q(3,nx)
    real(4), intent(out), device :: dl(nx-2), d(nx-2), du(nx-2), b(nx-2)
    real(4) dudx, RHS, a0, alpha_a0
    real(4) over_rho_im1, over_rho_i, over_rho_ip1
    real(4) over_rhox1, over_rhox2
    integer i, ii
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    if (i < 2 .or. i > nx-1) return
    ii = i - 1
    dudx  = 0.5e0 * real(-Q(2,i-1) + Q(2,i+1)) * over_dx
    RHS   = 2.e0 * (dudx**2)
    over_rho_im1 = 1.e0 / real(Q(1,i-1))
    over_rho_i   = 1.e0 / real(Q(1,i))
    over_rho_ip1 = 1.e0 / real(Q(1,i+1))
    over_rhox1   = 0.5e0 * (over_rho_im1 + over_rho_i)
    over_rhox2   = 0.5e0 * (over_rho_i   + over_rho_ip1)
    a0 = over_rho_i + alpha * (over_rhox1 + over_rhox2) * over_dx**2
    alpha_a0 = alpha / a0
    ! band coefficients
    d(ii)  = 1.0e0
    dl(ii) = - over_rhox1 * over_dx**2 * alpha_a0
    du(ii) = - over_rhox2 * over_dx**2 * alpha_a0
    b(ii)  = RHS * alpha_a0
    ! Neumann boundary condition
    if (ii == 1) then
      d(ii)  = d(ii) + dl(ii)
      dl(ii) = 0.e0
    endif
    if (ii == nx-2) then
      d(ii)  = d(ii) + du(ii)
      du(ii) = 0.e0
    endif
  end subroutine set_band


  subroutine calc_sigma(nx, Q, sigma)
    integer, intent(in), value     :: nx
    real(8), intent(in), device    :: Q(3,nx) ! rho, u, p
    real(4), intent(inout), device :: sigma(nx)
    real(4), allocatable, device :: dl(:), d(:), du(:), b(:), buffer(:)
    integer(8) buffersize
    type(cusparseHandle) handle
    integer n, istat
    n = nx - 2
    allocate(dl(n), d(n), du(n), b(n))
    call set_band<<<blocks,threads>>>(nx, Q, dl, d, du, b)
    istat = cusparseCreate(handle)
    istat = cusparseSgtsv2_bufferSizeExt(handle, n, 1, dl, d, du, b, n, bufferSize)
    allocate(buffer(bufferSize))
    istat = cusparseSgtsv2(handle, n, 1, dl, d, du, b, n, buffer)
    istat = cusparseDestroy(handle)
    sigma(2:nx-1) = b
    sigma(1)  = sigma(2)
    sigma(nx) = sigma(nx-1)
    deallocate(dl, d, du, b, buffer)
  end subroutine calc_sigma
end module calc_igr

