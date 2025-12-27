module calc_steps
  use cudafor
  use mod_globals, only : dt
  use mod_constant, only : one_sixth
  implicit none
contains
  subroutine calc_R(nx, ny, dx, dy, E, F, R)
    integer, intent(in), value                          :: nx, ny
    real(8), intent(in), dimension(nx-1), device        :: dx
    real(8), intent(in), dimension(ny-1), device        :: dy
    real(8), intent(in), dimension(4,nx-1,ny-2), device :: E
    real(8), intent(in), dimension(4,nx-2,ny-1), device :: F
    real(8), intent(out), dimension(4,nx-2,ny-2),device :: R
    real(8) dtdx, dtdy
    integer i, j, l
    !$cuf kernel do(2) <<<*,(32,4)>>>
    do j = 1, ny-2
      do i = 1, nx-2
        dtdx = dt * dx(i)
        dtdy = dt * dy(j)
        do l = 1, 4
          R(l,i,j) = &
          &   dtdy * (-E(l,i,j) + E(l,i+1,j)) &
          & + dtdx * (-F(l,i,j) + F(l,i,j+1))
    enddo;enddo;enddo
  end subroutine calc_R


  subroutine calc_step1(nx, ny, coef, dx, dy, E, F, Q, Q2)
    integer, intent(in), value                          :: nx, ny
    real(8), intent(in), value                          :: coef
    real(8), intent(in), dimension(nx-1), device        :: dx
    real(8), intent(in), dimension(ny-1), device        :: dy
    real(8), intent(in), dimension(4,nx-1,ny-2), device :: E
    real(8), intent(in), dimension(4,nx-2,ny-1), device :: F
    real(8), intent(in), dimension(4,nx,ny), device     :: Q
    real(8), intent(out), dimension(4,nx,ny), device    :: Q2
    real(8) R, dtdx, dtdy
    integer i, j, l
    !$cuf kernel do(2) <<<*,(32,4)>>>
    do j = 1, ny-2
      do i = 1, nx-2
        dtdx = dt * dx(i)
        dtdy = dt * dy(j)
        do l = 1, 4
            R = dtdy * (-E(l,i,j) + E(l,i+1,j)) &
            & + dtdx * (-F(l,i,j) + F(l,i,j+1))
            Q2(l,i+1,j+1) = Q(l,i+1,j+1) - coef * R
    enddo;enddo;enddo
  end subroutine calc_step1


  subroutine calc_step(nx, ny, coef1, coef2, dx, dy, E, F, Q, Q2, Rs)
    integer, intent(in), value                             :: nx, ny
    real(8), intent(in), value                             :: coef1, coef2
    real(8), intent(in), dimension(nx-1), device           :: dx
    real(8), intent(in), dimension(ny-1), device           :: dy
    real(8), intent(in), dimension(4,nx-1,ny-2), device    :: E
    real(8), intent(in), dimension(4,nx-2,ny-1), device    :: F
    real(8), intent(in), dimension(4,nx,ny), device        :: Q
    real(8), intent(out), dimension(4,nx,ny), device       :: Q2
    real(8), intent(inout), dimension(4,nx-2,ny-2), device :: Rs
    real(8) R, dtdx, dtdy
    integer i, j, l
    !$cuf kernel do(2) <<<*,(32,4)>>>
    do j = 1, ny-2
      do i = 1, nx-2
        dtdx = dt * dx(i)
        dtdy = dt * dy(j)
        do l = 1, 4
            R = dtdy * (-E(l,i,j) + E(l,i+1,j)) &
            & + dtdx * (-F(l,i,j) + F(l,i,j+1))
            Q2(l,i+1,j+1) = Q(l,i+1,j+1) - coef1 * R
            Rs(l,i,j) = Rs(l,i,j) + coef2 * R
    enddo;enddo;enddo
  end subroutine calc_step
  
 
  subroutine calc_step2_3(nx, ny, coef1, coef2, coef3, coef4, dx, dy, E, F, Qin, Qout)
    integer, intent(in), value                          :: nx, ny
    real(8), intent(in), value                          :: coef1, coef2, coef3, coef4
    real(8), intent(in), dimension(nx-1), device        :: dx
    real(8), intent(in), dimension(ny-1), device        :: dy
    real(8), intent(in), dimension(4,nx-1,ny-2), device :: E
    real(8), intent(in), dimension(4,nx-2,ny-1), device :: F
    real(8), intent(in), dimension(4,nx,ny), device     :: Qin
    real(8), intent(inout), dimension(4,nx,ny), device  :: Qout
    real(8) R, dtdx, dtdy
    integer i, j, l
    !$cuf kernel do(2) <<<*,(32,4)>>>
    do j = 1, ny-2
      do i = 1, nx-2
        dtdx = dt * dx(i)
        dtdy = dt * dy(j)
        do l = 1, 4
            R = dtdy * (-E(l,i,j) + E(l,i+1,j)) &
            & + dtdx * (-F(l,i,j) + F(l,i,j+1))
            Qout(l,i+1,j+1) = (coef1 * Qin(l,i+1,j+1) + coef2 * Qout(l,i+1,j+1) - coef3 * R) / coef4
    enddo;enddo;enddo
  end subroutine calc_step2_3
  
 
  subroutine calc_step4(nx, ny, dx, dy, E, F, Rs, Q)
    integer, intent(in), value                             :: nx, ny
    real(8), intent(in), dimension(nx-1), device           :: dx
    real(8), intent(in), dimension(ny-1), device           :: dy
    real(8), intent(in), dimension(4,nx-1,ny-2), device    :: E
    real(8), intent(in), dimension(4,nx-2,ny-1), device    :: F
    real(8), intent(inout), dimension(4,nx-2,ny-2), device :: Rs
    real(8), intent(inout), dimension(4,nx,ny), device     :: Q
    real(8) R, dtdx, dtdy
    integer i, j, l
    !$cuf kernel do(2) <<<*,(32,4)>>>
    do j = 1, ny-2
      do i = 1, nx-2
        dtdx = dt * dx(i)
        dtdy = dt * dy(j)
        do l = 1, 4
            R = dtdy * (-E(l,i,j) + E(l,i+1,j)) &
            & + dtdx * (-F(l,i,j) + F(l,i,j+1))
            Rs(l,i,j) = Rs(l,i,j) + R
            Q(l,i+1,j+1) = Q(l,i+1,j+1) - Rs(l,i,j) * one_sixth
            Rs(l,i,j) = 0.d0
    enddo;enddo;enddo
  end subroutine calc_step4


  subroutine calc_error(nx, ny, R1, R2, R1_new, R2_new, err)
    integer, intent(in), value                          :: nx, ny
    real(8), intent(in), dimension(4,nx-2,ny-2), device :: R1, R2, R1_new, R2_new
    real(8), intent(out)                                :: err
    integer i, j, l
    err = 0.d0
    !$cuf kernel do(2) <<<*,(32,4,1)>>>
    do j = 1, ny-2
      do i = 1, nx-2
        do l = 1, 4
          err = err + sqrt((R1(l,i,j) - R1_new(l,i,j)**2)) + sqrt((R2(l,i,j) - R2_new(l,i,j))**2)
    enddo;enddo;enddo
    err = err / (dble(nx - 2) * dble(ny * 2) * 4.d0)
  end subroutine calc_error


  subroutine calc_Gauss_step(nx, ny, a1, a2, dx, dy, R1, R2, Q, Q2)
    integer, intent(in), value                          :: nx, ny
    real(8), intent(in), value                          :: a1, a2
    real(8), intent(in), dimension(nx-1), device        :: dx
    real(8), intent(in), dimension(ny-1), device        :: dy
    real(8), intent(in), dimension(4,nx-2,ny-2), device :: R1
    real(8), intent(in), dimension(4,nx-2,ny-2), device :: R2
    real(8), intent(in), dimension(4,nx,ny), device     :: Q
    real(8), intent(out), dimension(4,nx,ny), device    :: Q2
    integer i, j, l
    !$cuf kernel do(2) <<<*,(32,4)>>>
    do j = 2, ny-1
      do i = 2, nx-1
        do l = 1, 4
            Q2(l,i,j) = Q(l,i,j) - (a1 * R1(l,i-1,j-1) + a2 * R2(l,i-1,j-1))
    enddo;enddo;enddo
  end subroutine calc_Gauss_step


  subroutine calc_Gauss_step_Q(nx, ny, a1, a2, dx, dy, R1, R2, Q)
    integer, intent(in), value                          :: nx, ny
    real(8), intent(in), value                          :: a1, a2
    real(8), intent(in), dimension(nx-1), device        :: dx
    real(8), intent(in), dimension(ny-1), device        :: dy
    real(8), intent(in), dimension(4,nx-2,ny-2), device :: R1
    real(8), intent(in), dimension(4,nx-2,ny-2), device :: R2
    real(8), intent(inout), dimension(4,nx,ny), device  :: Q
    integer i, j, l
    !$cuf kernel do(2) <<<*,(32,4)>>>
    do j = 2, ny-1
      do i = 2, nx-1
        do l = 1, 4
          Q(l,i,j) = Q(l,i,j) - (a1 * R1(l,i-1,j-1) + a2 * R2(l,i-1,j-1))
    enddo;enddo;enddo
  end subroutine calc_Gauss_step_Q
end module calc_steps

