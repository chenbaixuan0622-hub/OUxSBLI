module calc_steps
  use cudafor
  use mod_globals, only : accuracy, offset, dt
  implicit none
  interface
    subroutine calc_step(nx,coef1,coef2,dx,E,Q,Q2,Rs)
      integer, intent(in), value                                         :: nx
      real(8), intent(in), value                                         :: coef1, coef2
      real(8), intent(in), dimension(nx-1), device                       :: dx
      real(8), intent(in), dimension(nx-accuracy+1,3), device            :: E
      real(8), intent(in), dimension(nx,3), device                       :: Q
      real(8), intent(out), dimension(nx,3), device                      :: Q2
      real(8), intent(inout), dimension(nx-accuracy,3), optional, device :: Rs
    end subroutine
  end interface
  interface
    subroutine calc_step2(nx,coef1,coef2,coef3,coef4,dx,E,Q,Q2,Q3,Rs)
      integer, intent(in), value                                         :: nx
      real(8), intent(in), value                                         :: coef1, coef2, coef3, coef4
      real(8), intent(in), dimension(nx-1), device                       :: dx
      real(8), intent(in), dimension(nx-accuracy+1,3), device            :: E
      real(8), intent(in), dimension(nx,3), device                       :: Q, Q2
      real(8), intent(out), dimension(nx,3), device                      :: Q3
      real(8), intent(inout), dimension(nx-accuracy,3), optional, device :: Rs
    end subroutine
  end interface
contains
  subroutine calc_step(nx,coef1,coef2,dx,E,Q,Q2,Rs)
    integer, intent(in), value                                         :: nx
    real(8), intent(in), value                                         :: coef1, coef2
    real(8), intent(in), dimension(nx-1), device                       :: dx ! 1 / dx
    real(8), intent(in), dimension(nx-accuracy+1,3), device            :: E
    real(8), intent(in), dimension(nx,3), device                       :: Q
    real(8), intent(out), dimension(nx,3), device                      :: Q2
    real(8), intent(inout), dimension(nx-accuracy,3), optional, device :: Rs
    integer i, j
    real(8) R
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, 3
      do i = 1+offset, nx-offset
        R = dt * dx(i-offset) * (-E(i-offset,j) + E(i-offset+1,j))
        Q2(i,j) = Q(i,j) - coef1 * R
        if (present(Rs)) then
          Rs(i-offset,j) = Rs(i-offset,j) + coef2 * R
        endif
    enddo;enddo
  end subroutine calc_step
  
  subroutine calc_step2(nx,coef1,coef2,coef3,coef4,dx,E,Q,Q2,Q3,Rs)
    integer, intent(in), value                                         :: nx
    real(8), intent(in), value                                         :: coef1, coef2, coef3, coef4
    real(8), intent(in), dimension(nx-1), device                       :: dx ! 1 / dx
    real(8), intent(in), dimension(nx-accuracy+1,3), device            :: E
    real(8), intent(in), dimension(nx,3), device                       :: Q, Q2
    real(8), intent(out), dimension(nx,3), device                      :: Q3
    real(8), intent(inout), dimension(nx-accuracy,3), optional, device :: Rs
    integer i, j
    real(8) R
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, 3 
      do i = 1+offset, nx-offset
        R = dt * dx(i-offset) * (-E(i-offset,j) + E(i-offset+1,j))
        Q3(i,j) = (coef1 * Q(i,j) + coef2 * Q2(i,j) - coef3 * R) / coef4
        if (present(Rs)) then
          Rs(i-offset,j) = R
        endif
    enddo;enddo
  end subroutine calc_step2
  
  subroutine calc_step3(nx,dx,E,Q3,Q)
    integer, intent(in), value                              :: nx
    real(8), intent(in), dimension(nx-1), device            :: dx ! 1 / dx
    real(8), intent(in), dimension(nx-accuracy+1,3), device :: E
    real(8), intent(in), dimension(nx,3), device            :: Q3
    real(8), intent(inout), dimension(nx,3), device         :: Q
    integer i, j
    real(8) R
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, 3
      do i = 1+offset, nx-offset
        R = dt * dx(i-offset) * (-E(i-offset,j) + E(i-offset+1,j))
        Q(i,j) = (Q(i,j) + 2.0d0 * Q3(i,j) - 2.d0 * R) / 3.d0
    enddo;enddo
  end subroutine calc_step3
 
  subroutine calc_step4(nx,dx,E,Rs,Q)
    integer, intent(in), value                               :: nx
    real(8), intent(in), dimension(nx-1), device             :: dx ! 1 / dx
    real(8), intent(in), dimension(nx-accuracy+1,3), device  :: E
    real(8), intent(inout), dimension(nx-accuracy,3), device :: Rs
    real(8), intent(inout), dimension(nx,3), device          :: Q
    integer i, j
    real(8) R
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, 3
      do i = 1+offset, nx-offset
        R = dt * dx(i-offset) * (-E(i-offset,j) + E(i-offset+1,j))
        Rs(i-offset,j) = Rs(i-offset,j) + R
        Q(i,j) = Q(i,j) - Rs(i-offset,j) / 6.d0
        Rs(i-offset,j) = 0.d0
    enddo;enddo
  end subroutine calc_step4
end module calc_steps

