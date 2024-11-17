module calc_steps
  use cudafor
  use mod_globals, only : accuracy, offset, dt
  implicit none
contains
  subroutine calc_step(nx,ny,coef1,coef2,dx,dy,E,F,Q,Q2,Rs)
    integer, intent(in), value                                           :: nx, ny
    real(8), intent(in), value                                           :: coef1, coef2
    real(8), intent(in), dimension(nx-1), device                         :: dx
    real(8), intent(in), dimension(ny-1), device                         :: dy
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device  :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device  :: F
    real(8), intent(in), dimension(nx,ny,4), device                      :: Q
    real(8), intent(out), dimension(nx,ny,4), device                     :: Q2
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,4), device :: Rs
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          R = dt / dy(j-offset) * (-E(i-offset,j-offset,k) + E(i-offset+1,j-offset,k)) &
          & + dt / dx(i-offset) * (-F(i-offset,j-offset,k) + F(i-offset,j-offset+1,k))
          Q2(i,j,k) = Q(i,j,k) - coef1 * R
          Rs(i-offset,j-offset,k) = Rs(i-offset,j-offset,k) + coef2 * R
    enddo;enddo;enddo
  end subroutine calc_step
  
  subroutine calc_step4(nx,ny,dx,dy,E,F,Rs,Q)
    integer, intent(in), value                                           :: nx, ny
    real(8), intent(in), dimension(nx-1), device                         :: dx
    real(8), intent(in), dimension(ny-1), device                         :: dy
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device  :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device  :: F
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,4), device :: Rs
    real(8), intent(inout), dimension(nx,ny,4), device                   :: Q
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          R = dt / dy(j-offset) * (-E(i-offset,j-offset,k) + E(i-offset+1,j-offset,k)) &
          & + dt / dx(i-offset) * (-F(i-offset,j-offset,k) + F(i-offset,j-offset+1,k))
          Rs(i-offset,j-offset,k) = Rs(i-offset,j-offset,k) + R
          Q(i,j,k) = Q(i,j,k) - Rs(i-offset,j-offset,k) / 6.d0
          Rs(i-offset,j-offset,k) = 0.d0
    enddo;enddo;enddo
  end subroutine calc_step4
end module calc_steps

