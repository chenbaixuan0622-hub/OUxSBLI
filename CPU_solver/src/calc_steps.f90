module calc_steps
  use mod_globals, only : accuracy, offset, dt
  implicit none
contains
  subroutine calc_step(nx,ny,nz,coef1,coef2,dx,dy,dz,E,F,G,Q,Q2,Rs)
    integer, intent(in)                                                      :: nx, ny, nz
    real(8), intent(in)                                                      :: coef1, coef2
    real(8), intent(in), dimension(nx-1)                                     :: dx
    real(8), intent(in), dimension(ny-1)                                     :: dy
    real(8), intent(in), dimension(nz-1)                                     :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5)  :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5)  :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5)  :: G
    real(8), intent(in), dimension(nx,ny,nz,5)                               :: Q
    real(8), intent(out), dimension(nx,ny,nz,5)                              :: Q2
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5) :: Rs
    integer i, j, k, l
    real(8) R
    do l = 1, 5
      do k = 1+offset, nz-offset
        do j = 1+offset, ny-offset
          do i = 1+offset, nx-offset
            R = dt * dy(j-offset) * dz(k-offset) * (-E(i-offset,j-offset,k-offset,l) + E(i-offset+1,j-offset,k-offset,l)) &
            & + dt * dz(k-offset) * dx(i-offset) * (-F(i-offset,j-offset,k-offset,l) + F(i-offset,j-offset+1,k-offset,l)) &
            & + dt * dx(i-offset) * dy(j-offset) * (-G(i-offset,j-offset,k-offset,l) + G(i-offset,j-offset,k-offset+1,l))
            Q2(i,j,k,l) = Q(i,j,k,l) - coef1 * R
            Rs(i-offset,j-offset,k-offset,l) = Rs(i-offset,j-offset,k-offset,l) + coef2 * R
    enddo;enddo;enddo;enddo
  end subroutine calc_step
  
  subroutine calc_step4(nx,ny,nz,dx,dy,dz,E,F,G,Rs,Q)
    integer, intent(in)                                                      :: nx, ny, nz
    real(8), intent(in), dimension(nx-1)                                     :: dx
    real(8), intent(in), dimension(ny-1)                                     :: dy
    real(8), intent(in), dimension(nz-1)                                     :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5)  :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5)  :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5)  :: G
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5) :: Rs
    real(8), intent(inout), dimension(nx,ny,nz,5)                            :: Q
    integer i, j, k, l
    real(8) R
    do l = 1, 5
      do k = 1+offset, nz-offset
        do j = 1+offset, ny-offset
          do i = 1+offset, nx-offset
            R = dt * dy(j-offset) * dz(k-offset) * (-E(i-offset,j-offset,k-offset,l) + E(i-offset+1,j-offset,k-offset,l)) &
            & + dt * dz(k-offset) * dx(i-offset) * (-F(i-offset,j-offset,k-offset,l) + F(i-offset,j-offset+1,k-offset,l)) &
            & + dt * dx(i-offset) * dy(j-offset) * (-G(i-offset,j-offset,k-offset,l) + G(i-offset,j-offset,k-offset+1,l))
            Rs(i-offset,j-offset,k-offset,l) = Rs(i-offset,j-offset,k-offset,l) + R
            Q(i,j,k,l) = Q(i,j,k,l) - Rs(i-offset,j-offset,k-offset,l) / 6.d0
            Rs(i-offset,j-offset,k-offset,l) = 0.d0
    enddo;enddo;enddo;enddo
  end subroutine calc_step4
end module calc_steps

