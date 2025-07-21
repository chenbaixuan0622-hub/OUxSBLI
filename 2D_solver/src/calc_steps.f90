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
      real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device            :: E
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device            :: F
      real(8), intent(in), dimension(nx,ny,4), device                                :: Q
      real(8), intent(out), dimension(nx,ny,4), device                               :: Q2
      real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,4), optional, device :: Rs
    end subroutine
  end interface
  
  interface
    subroutine calc_step_forcing(nx, ny, coef1, coef2, dx, dy, E, F, fx, fy, Q, Q2, Rs)
      integer, intent(in), value                                                     :: nx, ny
      real(8), intent(in), value                                                     :: coef1, coef2
      real(8), intent(in), dimension(nx-1), device                                   :: dx
      real(8), intent(in), dimension(ny-1), device                                   :: dy
      real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device            :: E
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device            :: F
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device              :: fx
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device              :: fy
      real(8), intent(in), dimension(nx,ny,4), device                                :: Q
      real(8), intent(out), dimension(nx,ny,4), device                               :: Q2
      real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,4), optional, device :: Rs
    end subroutine
  end interface
  
  interface
    subroutine calc_step2(nx, ny, coef1, coef2, coef3, coef4, dx, dy, E, F, Q, Q2, Rs)
      integer, intent(in), value                                                     :: nx, ny
      real(8), intent(in), value                                                     :: coef1, coef2, coef3, coef4
      real(8), intent(in), dimension(nx-1), device                                   :: dx
      real(8), intent(in), dimension(ny-1), device                                   :: dy
      real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device            :: E
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device            :: F
      real(8), intent(in), dimension(nx,ny,4), device                                :: Q
      real(8), intent(inout), dimension(nx,ny,4), device                             :: Q2
      real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,4), optional, device :: Rs
    end subroutine
  end interface
  
  interface
    subroutine calc_step2_forcing(nx, ny, coef1, coef2, coef3, coef4, dx, dy, E, F, fx, fy, Q, Q2, Rs)
      integer, intent(in), value                                                     :: nx, ny
      real(8), intent(in), value                                                     :: coef1, coef2, coef3, coef4
      real(8), intent(in), dimension(nx-1), device                                   :: dx
      real(8), intent(in), dimension(ny-1), device                                   :: dy
      real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device            :: E
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device            :: F
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device              :: fx
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device              :: fy
      real(8), intent(in), dimension(nx,ny,4), device                                :: Q
      real(8), intent(inout), dimension(nx,ny,4), device                             :: Q2
      real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,4), optional, device :: Rs
    end subroutine
  end interface
contains
  subroutine calc_R(nx, ny, dx, dy, E, F, R)
    integer, intent(in), value                                          :: nx, ny
    real(8), intent(in), dimension(nx-1), device                        :: dx
    real(8), intent(in), dimension(ny-1), device                        :: dy
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy,4),device   :: R
    integer i, j, k
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          R(i,j,k) = &
          &   dt / dy(j) * (-E(i,j,k) + E(i+1,j,k)) &
          & + dt / dx(i) * (-F(i,j,k) + F(i,j+1,k))
    enddo;enddo;enddo
  end subroutine calc_R

  subroutine calc_R_forcing(nx, ny, dx, dy, E, F, fx, fy, R)
    integer, intent(in), value                                          :: nx, ny
    real(8), intent(in), dimension(nx-1), device                        :: dx
    real(8), intent(in), dimension(ny-1), device                        :: dy
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device   :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device   :: fy
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy,4),device   :: R
    integer i, j, k
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          R(i,j,k) = &
          &   dt / dy(j) * (-E(i,j,k) + E(i+1,j,k) + fx(i,j,k) / dx(i)) &
          & + dt / dx(i) * (-F(i,j,k) + F(i,j+1,k) + fy(i,j,k) / dy(j))
    enddo;enddo;enddo
  end subroutine calc_R_forcing

  subroutine calc_step(nx, ny, coef1, coef2, dx, dy, E, F, Q, Q2, Rs)
    integer, intent(in), value                                                     :: nx, ny
    real(8), intent(in), value                                                     :: coef1, coef2
    real(8), intent(in), dimension(nx-1), device                                   :: dx
    real(8), intent(in), dimension(ny-1), device                                   :: dy
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device            :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device            :: F
    real(8), intent(in), dimension(nx,ny,4), device                                :: Q
    real(8), intent(out), dimension(nx,ny,4), device                               :: Q2
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,4), optional, device :: Rs
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          R = dt / dy(j) * (-E(i,j,k) + E(i+1,j,k)) &
          & + dt / dx(i) * (-F(i,j,k) + F(i,j+1,k))
          Q2(i+offset,j+offset,k) = Q(i+offset,j+offset,k) - coef1 * R
          if (present(Rs)) then
            Rs(i,j,k) = Rs(i,j,k) + coef2 * R
          endif
    enddo;enddo;enddo
  end subroutine calc_step
  
  subroutine calc_step_forcing(nx, ny, coef1, coef2, dx, dy, E, F, fx, fy, Q, Q2, Rs)
    integer, intent(in), value                                                     :: nx, ny
    real(8), intent(in), value                                                     :: coef1, coef2
    real(8), intent(in), dimension(nx-1), device                                   :: dx
    real(8), intent(in), dimension(ny-1), device                                   :: dy
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device            :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device            :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device              :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device              :: fy
    real(8), intent(in), dimension(nx,ny,4), device                                :: Q
    real(8), intent(out), dimension(nx,ny,4), device                               :: Q2
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,4), optional, device :: Rs
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          R = dt / dy(j) * (-E(i,j,k) + E(i+1,j,k) + fx(i,j,k) / dx(i)) &
          & + dt / dx(i) * (-F(i,j,k) + F(i,j+1,k) + fy(i,j,k) / dy(j))
          Q2(i+offset,j+offset,k) = Q(i+offset,j+offset,k) - coef1 * R
          if (present(Rs)) then
            Rs(i,j,k) = Rs(i,j,k) + coef2 * R
          endif
    enddo;enddo;enddo
  end subroutine calc_step_forcing
  
  subroutine calc_step2(nx, ny, coef1, coef2, coef3, coef4, dx, dy, E, F, Q, Q2, Rs)
    integer, intent(in), value                                                     :: nx, ny
    real(8), intent(in), value                                                     :: coef1, coef2, coef3, coef4
    real(8), intent(in), dimension(nx-1), device                                   :: dx
    real(8), intent(in), dimension(ny-1), device                                   :: dy
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device            :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device            :: F
    real(8), intent(in), dimension(nx,ny,4), device                                :: Q
    real(8), intent(inout), dimension(nx,ny,4), device                             :: Q2
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,4), optional, device :: Rs
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          R = dt / dy(j) * (-E(i,j,k) + E(i+1,j,k)) &
          & + dt / dx(i) * (-F(i,j,k) + F(i,j+1,k))
          Q2(i+offset,j+offset,k) = (coef1 * Q(i+offset,j+offset,k) + coef2 * Q2(i+offset,j+offset,k) - coef3 * R) / coef4
          if (present(Rs)) then
            Rs(i,j,k) = R
          endif
    enddo;enddo;enddo
  end subroutine calc_step2
  
  subroutine calc_step2_forcing(nx, ny, coef1, coef2, coef3, coef4, dx, dy, E, F, fx, fy, Q, Q2, Rs)
    integer, intent(in), value                                                     :: nx, ny
    real(8), intent(in), value                                                     :: coef1, coef2, coef3, coef4
    real(8), intent(in), dimension(nx-1), device                                   :: dx
    real(8), intent(in), dimension(ny-1), device                                   :: dy
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device            :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device            :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device              :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device              :: fy
    real(8), intent(in), dimension(nx,ny,4), device                                :: Q
    real(8), intent(inout), dimension(nx,ny,4), device                             :: Q2
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,4), optional, device :: Rs
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          R = dt / dy(j) * (-E(i,j,k) + E(i+1,j,k) + fx(i,j,k) / dx(i)) &
          & + dt / dx(i) * (-F(i,j,k) + F(i,j+1,k) + fy(i,j,k) / dy(j))
          Q2(i+offset,j+offset,k) = (coef1 * Q(i+offset,j+offset,k) + coef2 * Q2(i+offset,j+offset,k) - coef3 * R) / coef4
          if (present(Rs)) then
            Rs(i,j,k) = R
          endif
    enddo;enddo;enddo
  end subroutine calc_step2_forcing
  
  subroutine calc_step3(nx, ny, dx, dy, E, F, Q3, Q)
    integer, intent(in), value                                          :: nx, ny
    real(8), intent(in), dimension(nx-1), device                        :: dx
    real(8), intent(in), dimension(ny-1), device                        :: dy
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    real(8), intent(in), dimension(nx,ny,4), device                     :: Q3
    real(8), intent(inout), dimension(nx,ny,4), device                  :: Q
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          R = dt / dy(j) * (-E(i,j,k) + E(i+1,j,k)) &
          & + dt / dx(i) * (-F(i,j,k) + F(i,j+1,k))
          Q(i+offset,j+offset,k) = (Q(i+offset,j+offset,k) + 2.0d0 * Q3(i+offset,j+offset,k) - 2.d0 * R) / 3.d0
    enddo;enddo;enddo
  end subroutine calc_step3
 
  subroutine calc_step3_forcing(nx, ny, dx, dy, E, F, fx, fy, Q3, Q)
    integer, intent(in), value                                          :: nx, ny
    real(8), intent(in), dimension(nx-1), device                        :: dx
    real(8), intent(in), dimension(ny-1), device                        :: dy
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device   :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device   :: fy
    real(8), intent(in), dimension(nx,ny,4), device                     :: Q3
    real(8), intent(inout), dimension(nx,ny,4), device                  :: Q
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          R = dt / dy(j) * (-E(i,j,k) + E(i+1,j,k) + fx(i,j,k) / dx(i)) &
          & + dt / dx(i) * (-F(i,j,k) + F(i,j+1,k) + fy(i,j,k) / dy(j))
          Q(i+offset,j+offset,k) = (Q(i+offset,j+offset,k) + 2.0d0 * Q3(i+offset,j+offset,k) - 2.d0 * R) / 3.d0
    enddo;enddo;enddo
  end subroutine calc_step3_forcing
 
  subroutine calc_step4(nx, ny, dx, dy, E, F, Rs, Q)
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
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          R = dt / dy(j) * (-E(i,j,k) + E(i+1,j,k)) &
          & + dt / dx(i) * (-F(i,j,k) + F(i,j+1,k))
          Rs(i,j,k) = Rs(i,j,k) + R
          Q(i+offset,j+offset,k) = Q(i+offset,j+offset,k) - Rs(i,j,k) / 6.d0
          Rs(i,j,k) = 0.d0
    enddo;enddo;enddo
  end subroutine calc_step4

  subroutine calc_step4_forcing(nx, ny, dx, dy, E, F, fx, fy, Rs, Q)
    integer, intent(in), value                                           :: nx, ny
    real(8), intent(in), dimension(nx-1), device                         :: dx
    real(8), intent(in), dimension(ny-1), device                         :: dy
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device  :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device  :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device    :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device    :: fy
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,4), device :: Rs
    real(8), intent(inout), dimension(nx,ny,4), device                   :: Q
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          R = dt / dy(j) * (-E(i,j,k) + E(i+1,j,k) + fx(i,j,k) / dx(i)) &
          & + dt / dx(i) * (-F(i,j,k) + F(i,j+1,k) + fy(i,j,k) / dy(j))
          Rs(i,j,k) = Rs(i,j,k) + R
          Q(i+offset,j+offset,k) = Q(i+offset,j+offset,k) - Rs(i,j,k) / 6.d0
          Rs(i,j,k) = 0.d0
    enddo;enddo;enddo
  end subroutine calc_step4_forcing

  subroutine calc_error(nx, ny, R1, R2, R1_new, R2_new, err)
    integer, intent(in), value                                        :: nx, ny
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device :: R1, R2, R1_new, R2_new
    real(8), intent(out)                                              :: err
    integer i, j, k
    err = 0.d0
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          err = err + sqrt((R1(i,j,k) - R1_new(i,j,k)**2)) + sqrt((R2(i,j,k) - R2_new(i,j,k))**2)
    enddo;enddo;enddo
    err = err / (dble(nx - accuracy) * dble(ny * accuracy) * 4.d0)
  end subroutine calc_error

  subroutine calc_Gauss_step(nx, ny, a1, a2, dx, dy, R1, R2, Q, Q2)
    integer, intent(in), value                                        :: nx, ny
    real(8), intent(in), value                                        :: a1, a2
    real(8), intent(in), dimension(nx-1), device                      :: dx
    real(8), intent(in), dimension(ny-1), device                      :: dy
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device :: R1
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device :: R2
    real(8), intent(in), dimension(nx,ny,4), device                   :: Q
    real(8), intent(out), dimension(nx,ny,4), device                  :: Q2
    integer i, j, k
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          Q2(i,j,k) = Q(i,j,k) - (a1 * R1(i-offset,j-offset,k) + a2 * R2(i-offset,j-offset,k))
    enddo;enddo;enddo
  end subroutine calc_Gauss_step
  
  subroutine calc_Gauss_step_Q(nx, ny, a1, a2, dx, dy, R1, R2, Q)
    integer, intent(in), value                                        :: nx, ny
    real(8), intent(in), value                                        :: a1, a2
    real(8), intent(in), dimension(nx-1), device                      :: dx
    real(8), intent(in), dimension(ny-1), device                      :: dy
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device :: R1
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,4), device :: R2
    real(8), intent(inout), dimension(nx,ny,4), device                :: Q
    integer i, j, k
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          Q(i,j,k) = Q(i,j,k) - (a1 * R1(i-offset,j-offset,k) + a2 * R2(i-offset,j-offset,k))
    enddo;enddo;enddo
  end subroutine calc_Gauss_step_Q
end module calc_steps

