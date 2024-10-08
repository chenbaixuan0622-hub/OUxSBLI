module calc_steps
  use cudafor
  use mod_globals, only : accuracy, offset, dt
  implicit none
  interface
    subroutine calc_step(nx,ny,nz,coef1,coef2,dx,dy,dz,E,F,G,Q,Q2,Rs)
      integer, intent(in), value                                                                  :: nx, ny, nz
      real(8), intent(in), value                                                                  :: coef1, coef2
      real(8), intent(in), dimension(nx-1), device                                                :: dx
      real(8), intent(in), dimension(ny-1), device                                                :: dy
      real(8), intent(in), dimension(nz-1), device                                                :: dz
      real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device             :: E
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device             :: F
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device             :: G
      real(8), intent(in), dimension(nx,ny,nz,5), device                                          :: Q
      real(8), intent(out), dimension(nx,ny,nz,5), device                                         :: Q2
      real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), optional, device  :: Rs
    end subroutine
  end interface
  interface
    subroutine calc_step2(nx,ny,nz,coef1,coef2,coef3,coef4,dx,dy,dz,E,F,G,Q,Q2,Q3,Rs)
      integer, intent(in), value                                                                  :: nx, ny, nz
      real(8), intent(in), value                                                                  :: coef1, coef2, coef3, coef4
      real(8), intent(in), dimension(nx-1), device                                                :: dx
      real(8), intent(in), dimension(ny-1), device                                                :: dy
      real(8), intent(in), dimension(nz-1), device                                                :: dz
      real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device             :: E
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device             :: F
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device             :: G
      real(8), intent(in), dimension(nx,ny,nz,5), device                                          :: Q, Q2
      real(8), intent(out), dimension(nx,ny,nz,5), device                                         :: Q3
      real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), optional, device  :: Rs
    end subroutine
  end interface
contains
  subroutine calc_step(nx,ny,nz,coef1,coef2,dx,dy,dz,E,F,G,Q,Q2,Rs)
    integer, intent(in), value                                                                  :: nx, ny, nz
    real(8), intent(in), value                                                                  :: coef1, coef2
    real(8), intent(in), dimension(nx-1), device                                                :: dx
    real(8), intent(in), dimension(ny-1), device                                                :: dy
    real(8), intent(in), dimension(nz-1), device                                                :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device             :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device             :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device             :: G
    real(8), intent(in), dimension(nx,ny,nz,5), device                                          :: Q
    real(8), intent(out), dimension(nx,ny,nz,5), device                                         :: Q2
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), optional, device  :: Rs
    integer i, j, k, l
    real(8) R
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1+offset, nz-offset
        do j = 1+offset, ny-offset
          do i = 1+offset, nx-offset
            R = dt / (dy(j-offset) * dz(k-offset)) * (-E(i-offset,j-offset,k-offset,l) + E(i-offset+1,j-offset,k-offset,l)) &
            & + dt / (dz(k-offset) * dx(i-offset)) * (-F(i-offset,j-offset,k-offset,l) + F(i-offset,j-offset+1,k-offset,l)) &
            & + dt / (dx(i-offset) * dy(j-offset)) * (-G(i-offset,j-offset,k-offset,l) + G(i-offset,j-offset,k-offset+1,l))
            Q2(i,j,k,l) = Q(i,j,k,l) - coef1 * R
            if (present(Rs)) then
              Rs(i-offset,j-offset,k-offset,l) = Rs(i-offset,j-offset,k-offset,l) + coef2 * R
            endif
    enddo;enddo;enddo;enddo
  end subroutine calc_step
  
  subroutine calc_step2(nx,ny,nz,coef1,coef2,coef3,coef4,dx,dy,dz,E,F,G,Q,Q2,Q3,Rs)
    integer, intent(in), value                                                                  :: nx, ny, nz
    real(8), intent(in), value                                                                  :: coef1, coef2, coef3, coef4
    real(8), intent(in), dimension(nx-1), device                                                :: dx
    real(8), intent(in), dimension(ny-1), device                                                :: dy
    real(8), intent(in), dimension(nz-1), device                                                :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device             :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device             :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device             :: G
    real(8), intent(in), dimension(nx,ny,nz,5), device                                          :: Q, Q2
    real(8), intent(out), dimension(nx,ny,nz,5), device                                         :: Q3
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), optional, device  :: Rs
    integer i, j, k, l
    real(8) R
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1+offset, nz-offset
        do j = 1+offset, ny-offset
          do i = 1+offset, nx-offset
            R = dt / (dy(j-offset) * dz(k-offset)) * (-E(i-offset,j-offset,k-offset,l) + E(i-offset+1,j-offset,k-offset,l)) &
            & + dt / (dz(k-offset) * dx(i-offset)) * (-F(i-offset,j-offset,k-offset,l) + F(i-offset,j-offset+1,k-offset,l)) &
            & + dt / (dx(i-offset) * dy(j-offset)) * (-G(i-offset,j-offset,k-offset,l) + G(i-offset,j-offset,k-offset+1,l))
            Q3(i,j,k,l) = (coef1 * Q(i,j,k,l) + coef2 * Q2(i,j,k,l) - coef3 * R) / coef4
            if (present(Rs)) then
              Rs(i-offset,j-offset,k-offset,l) = R
            endif
    enddo;enddo;enddo;enddo
  end subroutine calc_step2
  
  subroutine calc_step3(nx,ny,nz,dx,dy,dz,E,F,G,Q3,Q)
    integer, intent(in), value                                                      :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device                                    :: dx
    real(8), intent(in), dimension(ny-1), device                                    :: dy
    real(8), intent(in), dimension(nz-1), device                                    :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    real(8), intent(in), dimension(nx,ny,nz,5), device                              :: Q3
    real(8), intent(inout), dimension(nx,ny,nz,5), device                           :: Q
    integer i, j, k, l
    real(8) R
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1+offset, nz-offset
        do j = 1+offset, ny-offset
          do i = 1+offset, nx-offset
            R = dt / (dy(j-offset) * dz(k-offset)) * (-E(i-offset,j-offset,k-offset,l) + E(i-offset+1,j-offset,k-offset,l)) &
            & + dt / (dz(k-offset) * dx(i-offset)) * (-F(i-offset,j-offset,k-offset,l) + F(i-offset,j-offset+1,k-offset,l)) &
            & + dt / (dx(i-offset) * dy(j-offset)) * (-G(i-offset,j-offset,k-offset,l) + G(i-offset,j-offset,k-offset+1,l))
            Q(i,j,k,l) = (Q(i,j,k,l) + 2.0d0 * Q3(i,j,k,l) - 2.d0 * R) / 3.d0
    enddo;enddo;enddo;enddo
  end subroutine calc_step3
 
  subroutine calc_step4(nx,ny,nz,dx,dy,dz,E,F,G,Rs,Q)
    integer, intent(in), value                                                        :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device                                      :: dx
    real(8), intent(in), dimension(ny-1), device                                      :: dy
    real(8), intent(in), dimension(nz-1), device                                      :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device   :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device   :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device   :: G
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device  :: Rs
    real(8), intent(inout), dimension(nx,ny,nz,5), device                             :: Q
    integer i, j, k, l
    real(8) R
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1+offset, nz-offset
        do j = 1+offset, ny-offset
          do i = 1+offset, nx-offset
            R = dt / (dy(j-offset) * dz(k-offset)) * (-E(i-offset,j-offset,k-offset,l) + E(i-offset+1,j-offset,k-offset,l)) &
            & + dt / (dz(k-offset) * dx(i-offset)) * (-F(i-offset,j-offset,k-offset,l) + F(i-offset,j-offset+1,k-offset,l)) &
            & + dt / (dx(i-offset) * dy(j-offset)) * (-G(i-offset,j-offset,k-offset,l) + G(i-offset,j-offset,k-offset+1,l))
            Rs(i-offset,j-offset,k-offset,l) = Rs(i-offset,j-offset,k-offset,l) + R
            Q(i,j,k,l) = Q(i,j,k,l) - Rs(i-offset,j-offset,k-offset,l) / 6.d0
            Rs(i-offset,j-offset,k-offset,l) = 0.d0
    enddo;enddo;enddo;enddo
  end subroutine calc_step4
end module calc_steps

