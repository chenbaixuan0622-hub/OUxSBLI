module calc_forcing
  use cudafor
  use cufft
  implicit none
  integer plan_fwd, plan_inv
contains
  subroutine pack_uv(accuracy, nx, ny, Q, u, v)
    integer, intent(in), value   :: accuracy, nx, ny
    real(8), intent(in), device  :: Q(4,nx,ny)
    real(8), intent(out), device :: u(nx-accuracy,ny-accuracy)
    real(8), intent(out), device :: v(nx-accuracy,ny-accuracy)
    integer i, j, offset
    offset = accuracy / 2
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny-accuracy
      do i = 1, nx-accuracy
        u(i,j) = Q(2,i+offset,j+offset)
        v(i,j) = Q(3,i+offset,j+offset)
    enddo;enddo
  end subroutine pack_uv


  subroutine unpack_F(accuracy, nx, ny, F_x, F_y)
    integer, intent(in), value     :: accuracy, nx, ny
    real(8), intent(inout), device :: F_x(nx-accuracy,ny-accuracy)
    real(8), intent(inout), device :: F_y(nx-accuracy,ny-accuracy)
    real(8) norm
    integer i, j
    norm = 1.d0 / dble((nx-accuracy) * (ny-accuracy))
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny-accuracy
      do i = 1, nx-accuracy
        F_x(i,j) = F_x(i,j) * norm
        F_y(i,j) = F_y(i,j) * norm
    enddo;enddo
  end subroutine unpack_F


  subroutine calc_linear_forcing(accuracy, nx, ny, kmin, kmax, alpha, u_h, v_h, F_x, F_y)
    integer, intent(in), value      :: accuracy, nx, ny, kmin, kmax
    real(8), intent(in), value      :: alpha
    complex(8), intent(in), device  :: u_h((nx-accuracy)/2+1,ny-accuracy)
    complex(8), intent(in), device  :: v_h((nx-accuracy)/2+1,ny-accuracy)
    complex(8), intent(out), device :: F_x((nx-accuracy)/2+1,ny-accuracy)
    complex(8), intent(out), device :: F_y((nx-accuracy)/2+1,ny-accuracy)
    real(8) kx, ky, k2, kabs, dot_prod
    integer nx_a, ny_a, i, j
    nx_a = nx - accuracy
    ny_a = ny - accuracy
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny_a
      do i = 1, nx_a / 2 + 1
        kx = dble(i - 1)
        if (j <= ny_a/2 + 1) then
          ky = dble(j - 1)
        else
          ky = dble(j - 1 - ny_a)
        endif
        k2 = kx*kx + ky*ky
        kabs = sqrt(k2)
        if (kabs >= kmin .and. kabs <= kmax .and. k2 > 0.d0) then
          dot_prod = (kx * u_h(i,j) + ky * v_h(i,j)) / k2
          F_x(i,j) = alpha * (u_h(i,j) - dot_prod * kx)
          F_y(i,j) = alpha * (v_h(i,j) - dot_prod * ky)
        else
          F_x(i,j) = (0.d0, 0.d0)
          F_y(i,j) = (0.d0, 0.d0)
        endif
    enddo;enddo
  end subroutine calc_linear_forcing


  subroutine calc_forcing_main(accuracy, plan_fwd, plan_inv, nx, ny, kmin, kmax, alpha, Q, F_x, F_y)
    integer, intent(in), value   :: accuracy, plan_fwd, plan_inv, nx, ny, kmin, kmax
    real(8), intent(in), value   :: alpha
    real(8), intent(in), device  :: Q(4,nx,ny)
    real(8), intent(out), device :: F_x(nx-accuracy,ny-accuracy)
    real(8), intent(out), device :: F_y(nx-accuracy,ny-accuracy)
    real(8), device :: u(nx-accuracy,ny-accuracy), v(nx-accuracy,ny-accuracy)
    complex(8), device ::  u_c((nx-accuracy)/2+1,ny-accuracy),  v_c((nx-accuracy)/2+1,ny-accuracy)
    complex(8), device :: F_xc((nx-accuracy)/2+1,ny-accuracy), F_yc((nx-accuracy)/2+1,ny-accuracy)
    integer istat
    call pack_uv(accuracy, nx, ny, Q, u, v)
    istat = cufftExecD2Z(plan_fwd, u, u_c)
    istat = cufftExecD2Z(plan_fwd, v, v_c)
    call calc_linear_forcing(accuracy, nx, ny, kmin, kmax, alpha, u_c, v_c, F_xc, F_yc) 
    istat = cufftExecZ2D(plan_inv, F_xc, F_x)
    istat = cufftExecZ2D(plan_inv, F_yc, F_y)
    call unpack_F(accuracy, nx, ny, F_x, F_y)
  end subroutine calc_forcing_main
end module calc_forcing

