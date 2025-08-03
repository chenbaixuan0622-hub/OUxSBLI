module calc_steps
  use cudafor
  use mod_globals, only : accuracy, offset, dt
  implicit none
  interface
    subroutine calc_step(nx, ny, coef1, coef2, dx, dy, E, F, Q, Q2, Rs)
      integer, intent(in), value                                                     :: nx, ny
      real(8), intent(in), value                                                     :: coef1, coef2
      real(8), intent(in), dimension(nx-1), device                                   :: dx
      real(8), intent(in), dimension(ny-1), device                                   :: dy
      real(8), intent(in), dimension(4,nx-accuracy+1,ny-accuracy), device            :: E
      real(8), intent(in), dimension(4,nx-accuracy,ny-accuracy+1), device            :: F
      real(8), intent(in), dimension(4,nx,ny), device                                :: Q
      real(8), intent(out), dimension(4,nx,ny), device                               :: Q2
      real(8), intent(inout), dimension(4,nx-accuracy,ny-accuracy), optional, device :: Rs
    end subroutine
  end interface
  
  interface
    subroutine calc_step_forcing(nx, ny, coef1, coef2, dx, dy, E, F, fx, fy, Q, Q2, Rs)
      integer, intent(in), value                                                     :: nx, ny
      real(8), intent(in), value                                                     :: coef1, coef2
      real(8), intent(in), dimension(nx-1), device                                   :: dx
      real(8), intent(in), dimension(ny-1), device                                   :: dy
      real(8), intent(in), dimension(4,nx-accuracy+1,ny-accuracy), device            :: E
      real(8), intent(in), dimension(4,nx-accuracy,ny-accuracy+1), device            :: F
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy), device                :: fx
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy), device                :: fy
      real(8), intent(in), dimension(4,nx,ny), device                                :: Q
      real(8), intent(out), dimension(4,nx,ny), device                               :: Q2
      real(8), intent(inout), dimension(4,nx-accuracy,ny-accuracy), optional, device :: Rs
    end subroutine
  end interface
contains
  subroutine calc_R(nx, ny, dx, dy, E, F, R)
    integer, intent(in), value                                          :: nx, ny
    real(8), intent(in), dimension(nx-1), device                        :: dx
    real(8), intent(in), dimension(ny-1), device                        :: dy
    real(8), intent(in), dimension(4,nx-accuracy+1,ny-accuracy), device :: E
    real(8), intent(in), dimension(4,nx-accuracy,ny-accuracy+1), device :: F
    real(8), intent(out), dimension(4,nx-accuracy,ny-accuracy),device   :: R
    integer i, j, k
    !$cuf kernel do(3) <<<*,*>>>
    do j = 1, ny-2*offset
      do i = 1, nx-2*offset
        do k = 1, 4
          R(k,i,j) = &
          &   dt / dy(j) * (-E(k,i,j) + E(k,i+1,j)) &
          & + dt / dx(i) * (-F(k,i,j) + F(k,i,j+1))
    enddo;enddo;enddo
  end subroutine calc_R

  subroutine calc_R_forcing(nx, ny, dx, dy, E, F, fx, fy, R)
    integer, intent(in), value                                          :: nx, ny
    real(8), intent(in), dimension(nx-1), device                        :: dx
    real(8), intent(in), dimension(ny-1), device                        :: dy
    real(8), intent(in), dimension(4,nx-accuracy+1,ny-accuracy), device :: E
    real(8), intent(in), dimension(4,nx-accuracy,ny-accuracy+1), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy), device     :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy), device     :: fy
    real(8), intent(out), dimension(4,nx-accuracy,ny-accuracy),device   :: R
    integer i, j, k
    !$cuf kernel do(3) <<<*,*>>>
    do j = 1, ny-2*offset
      do i = 1, nx-2*offset
        do k = 1, 4
          R(k,i,j) = &
          &   dt / dy(j) * (-E(k,i,j) + E(k,i+1,j) + fx(i,j) / dx(i)) &
          & + dt / dx(i) * (-F(k,i,j) + F(k,i,j+1) + fy(i,j) / dy(j))
    enddo;enddo;enddo
  end subroutine calc_R_forcing

  subroutine calc_step(nx, ny, coef1, coef2, dx, dy, E, F, Q, Q2, Rs)
    integer, intent(in), value                                                     :: nx, ny
    real(8), intent(in), value                                                     :: coef1, coef2
    real(8), intent(in), dimension(nx-1), device                                   :: dx
    real(8), intent(in), dimension(ny-1), device                                   :: dy
    real(8), intent(in), dimension(4,nx-accuracy+1,ny-accuracy), device            :: E
    real(8), intent(in), dimension(4,nx-accuracy,ny-accuracy+1), device            :: F
    real(8), intent(in), dimension(4,nx,ny), device                                :: Q
    real(8), intent(out), dimension(4,nx,ny), device                               :: Q2
    real(8), intent(inout), dimension(4,nx-accuracy,ny-accuracy), optional, device :: Rs
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do j = 1, ny-2*offset
      do i = 1, nx-2*offset
        do k = 1, 4
          R = dt / dy(j) * (-E(k,i,j) + E(k,i+1,j)) &
          & + dt / dx(i) * (-F(k,i,j) + F(k,i,j+1))
          Q2(k,i+offset,j+offset) = Q(k,i+offset,j+offset) - coef1 * R
          if (present(Rs)) then
            Rs(k,i,j) = Rs(k,i,j) + coef2 * R
          endif
    enddo;enddo;enddo
  end subroutine calc_step
  
  subroutine calc_step_forcing(nx, ny, coef1, coef2, dx, dy, E, F, fx, fy, Q, Q2, Rs)
    integer, intent(in), value                                                     :: nx, ny
    real(8), intent(in), value                                                     :: coef1, coef2
    real(8), intent(in), dimension(nx-1), device                                   :: dx
    real(8), intent(in), dimension(ny-1), device                                   :: dy
    real(8), intent(in), dimension(4,nx-accuracy+1,ny-accuracy), device            :: E
    real(8), intent(in), dimension(4,nx-accuracy,ny-accuracy+1), device            :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy), device                :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy), device                :: fy
    real(8), intent(in), dimension(4,nx,ny), device                                :: Q
    real(8), intent(out), dimension(4,nx,ny), device                               :: Q2
    real(8), intent(inout), dimension(4,nx-accuracy,ny-accuracy), optional, device :: Rs
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do j = 1, ny-2*offset
      do i = 1, nx-2*offset
        do k = 1, 4
          R = dt / dy(j) * (-E(k,i,j) + E(k,i+1,j) + fx(i,j) / dx(i)) &
          & + dt / dx(i) * (-F(k,i,j) + F(k,i,j+1) + fy(i,j) / dy(j))
          Q2(k,i+offset,j+offset) = Q(k,i+offset,j+offset) - coef1 * R
          if (present(Rs)) then
            Rs(k,i,j) = Rs(k,i,j) + coef2 * R
          endif
    enddo;enddo;enddo
  end subroutine calc_step_forcing
  
  subroutine calc_step2_3(nx, ny, coef1, coef2, coef3, coef4, dx, dy, E, F, Q, Q2)
    integer, intent(in), value                                          :: nx, ny
    real(8), intent(in), value                                          :: coef1, coef2, coef3, coef4
    real(8), intent(in), dimension(nx-1), device                        :: dx
    real(8), intent(in), dimension(ny-1), device                        :: dy
    real(8), intent(in), dimension(4,nx-accuracy+1,ny-accuracy), device :: E
    real(8), intent(in), dimension(4,nx-accuracy,ny-accuracy+1), device :: F
    real(8), intent(in), dimension(4,nx,ny), device                     :: Q
    real(8), intent(inout), dimension(4,nx,ny), device                  :: Q2
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do j = 1, ny-2*offset
      do i = 1, nx-2*offset
        do k = 1, 4
          R = dt / dy(j) * (-E(k,i,j) + E(k,i+1,j)) &
          & + dt / dx(i) * (-F(k,i,j) + F(k,i,j+1))
          Q2(k,i+offset,j+offset) = (coef1 * Q(k,i+offset,j+offset) + coef2 * Q2(k,i+offset,j+offset) - coef3 * R) / coef4
    enddo;enddo;enddo
  end subroutine calc_step2_3
  
  subroutine calc_step2_3_forcing(nx, ny, coef1, coef2, coef3, coef4, dx, dy, E, F, fx, fy, Q, Q2)
    integer, intent(in), value                                          :: nx, ny
    real(8), intent(in), value                                          :: coef1, coef2, coef3, coef4
    real(8), intent(in), dimension(nx-1), device                        :: dx
    real(8), intent(in), dimension(ny-1), device                        :: dy
    real(8), intent(in), dimension(4,nx-accuracy+1,ny-accuracy), device :: E
    real(8), intent(in), dimension(4,nx-accuracy,ny-accuracy+1), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy), device     :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy), device     :: fy
    real(8), intent(in), dimension(4,nx,ny), device                     :: Q
    real(8), intent(inout), dimension(4,nx,ny), device                  :: Q2
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do j = 1, ny-2*offset
      do i = 1, nx-2*offset
        do k = 1, 4
          R = dt / dy(j) * (-E(k,i,j) + E(k,i+1,j) + fx(i,j) / dx(i)) &
          & + dt / dx(i) * (-F(k,i,j) + F(k,i,j+1) + fy(i,j) / dy(j))
          Q2(k,i+offset,j+offset) = (coef1 * Q(k,i+offset,j+offset) + coef2 * Q2(k,i+offset,j+offset) - coef3 * R) / coef4
    enddo;enddo;enddo
  end subroutine calc_step2_3_forcing
  
  subroutine calc_step4(nx, ny, dx, dy, E, F, Rs, Q)
    integer, intent(in), value                                           :: nx, ny
    real(8), intent(in), dimension(nx-1), device                         :: dx
    real(8), intent(in), dimension(ny-1), device                         :: dy
    real(8), intent(in), dimension(4,nx-accuracy+1,ny-accuracy), device  :: E
    real(8), intent(in), dimension(4,nx-accuracy,ny-accuracy+1), device  :: F
    real(8), intent(inout), dimension(4,nx-accuracy,ny-accuracy), device :: Rs
    real(8), intent(inout), dimension(4,nx,ny), device                   :: Q
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do j = 1, ny-2*offset
      do i = 1, nx-2*offset
        do k = 1, 4
          R = dt / dy(j) * (-E(k,i,j) + E(k,i+1,j)) &
          & + dt / dx(i) * (-F(k,i,j) + F(k,i,j+1))
          Rs(k,i,j) = Rs(k,i,j) + R
          Q(k,i+offset,j+offset) = Q(k,i+offset,j+offset) - Rs(k,i,j) / 6.d0
          Rs(k,i,j) = 0.d0
    enddo;enddo;enddo
  end subroutine calc_step4

  subroutine calc_step4_forcing(nx, ny, dx, dy, E, F, fx, fy, Rs, Q)
    integer, intent(in), value                                           :: nx, ny
    real(8), intent(in), dimension(nx-1), device                         :: dx
    real(8), intent(in), dimension(ny-1), device                         :: dy
    real(8), intent(in), dimension(4,nx-accuracy+1,ny-accuracy), device  :: E
    real(8), intent(in), dimension(4,nx-accuracy,ny-accuracy+1), device  :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy), device      :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy), device      :: fy
    real(8), intent(inout), dimension(4,nx-accuracy,ny-accuracy), device :: Rs
    real(8), intent(inout), dimension(4,nx,ny), device                   :: Q
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do j = 1, ny-2*offset
      do i = 1, nx-2*offset
        do k = 1, 4
          R = dt / dy(j) * (-E(k,i,j) + E(k,i+1,j) + fx(i,j) / dx(i)) &
          & + dt / dx(i) * (-F(k,i,j) + F(k,i,j+1) + fy(i,j) / dy(j))
          Rs(k,i,j) = Rs(k,i,j) + R
          Q(k,i+offset,j+offset) = Q(k,i+offset,j+offset) - Rs(k,i,j) / 6.d0
          Rs(k,i,j) = 0.d0
    enddo;enddo;enddo
  end subroutine calc_step4_forcing

  subroutine calc_error(nx, ny, R1, R2, R1_new, R2_new, err)
    integer, intent(in), value                                        :: nx, ny
    real(8), intent(in), dimension(4,nx-accuracy,ny-accuracy), device :: R1, R2, R1_new, R2_new
    real(8), intent(out)                                              :: err
    integer i, j, k
    err = 0.d0
    !$cuf kernel do(3) <<<*,*>>>
    do j = 1, ny-2*offset
      do i = 1, nx-2*offset
        do k = 1, 4
          err = err + sqrt((R1(k,i,j) - R1_new(k,i,j)**2)) + sqrt((R2(k,i,j) - R2_new(k,i,j))**2)
    enddo;enddo;enddo
    err = err / (dble(nx - accuracy) * dble(ny * accuracy) * 4.d0)
  end subroutine calc_error

  subroutine calc_Gauss_step(nx, ny, a1, a2, dx, dy, R1, R2, Q, Q2)
    integer, intent(in), value                                        :: nx, ny
    real(8), intent(in), value                                        :: a1, a2
    real(8), intent(in), dimension(nx-1), device                      :: dx
    real(8), intent(in), dimension(ny-1), device                      :: dy
    real(8), intent(in), dimension(4,nx-accuracy,ny-accuracy), device :: R1
    real(8), intent(in), dimension(4,nx-accuracy,ny-accuracy), device :: R2
    real(8), intent(in), dimension(4,nx,ny), device                   :: Q
    real(8), intent(out), dimension(4,nx,ny), device                  :: Q2
    integer i, j, k
    !$cuf kernel do(3) <<<*,*>>>
    do j = 1+offset, ny-offset
      do i = 1+offset, nx-offset
        do k = 1, 4
          Q2(k,i,j) = Q(k,i,j) - (a1 * R1(k,i-offset,j-offset) + a2 * R2(k,i-offset,j-offset))
    enddo;enddo;enddo
  end subroutine calc_Gauss_step
  
  subroutine calc_Gauss_step_Q(nx, ny, a1, a2, dx, dy, R1, R2, Q)
    integer, intent(in), value                                        :: nx, ny
    real(8), intent(in), value                                        :: a1, a2
    real(8), intent(in), dimension(nx-1), device                      :: dx
    real(8), intent(in), dimension(ny-1), device                      :: dy
    real(8), intent(in), dimension(4,nx-accuracy,ny-accuracy), device :: R1
    real(8), intent(in), dimension(4,nx-accuracy,ny-accuracy), device :: R2
    real(8), intent(inout), dimension(4,nx,ny), device                :: Q
    integer i, j, k
    !$cuf kernel do(3) <<<*,*>>>
    do j = 1+offset, ny-offset
      do i = 1+offset, nx-offset
        do k = 1, 4
          Q(k,i,j) = Q(k,i,j) - (a1 * R1(k,i-offset,j-offset) + a2 * R2(k,i-offset,j-offset))
    enddo;enddo;enddo
  end subroutine calc_Gauss_step_Q
end module calc_steps

