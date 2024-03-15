module calc_steps
  use mod_globals, only : accuracy
  implicit none
contains
  attributes(device) subroutine calc_R(dtdx,dtdy,E1,E2,F1,F2,R)
    real(8), intent(in), value :: dtdx, dtdy, E1, E2, F1, F2
    real(8), intent(out) :: R
    R = dtdx * (-E1 + E2) + dtdy * (-F1 + F2)
  end subroutine calc_R

  attributes(global) subroutine calc_step1(nx,ny,dtdx,dtdy,E,F,Q,Q2)
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: dtdx, dtdy
    real(8), intent(in), device :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(in), device :: F(nx-accuracy,ny-accuracy+1,4)
    real(8), intent(in), dimension(nx,ny,4), device :: Q
    real(8), intent(out), dimension(nx,ny,4), device :: Q2
    integer i, j, k
    integer :: offset = accuracy / 2
    real(8) R
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    do k = 1, 4
      call calc_R(dtdx,dtdy,E(i-offset,j-offset,k),E(i-offset+1,j-offset,k), &
              & F(i-offset,j-offset,k),F(i-offset,j-offset+1,k),R)
      Q2(i,j,k) = Q(i,j,k) - R
    enddo 
  end subroutine calc_step1
  
  attributes(global) subroutine calc_step2(nx,ny,dtdx,dtdy,E,F,Q,Q2,Q3)
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: dtdx, dtdy
    real(8), intent(in), device :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(in), device :: F(nx-accuracy,ny-accuracy+1,4)
    real(8), intent(in), dimension(nx,ny,4), device :: Q
    real(8), intent(in), dimension(nx,ny,4), device :: Q2
    real(8), intent(out), dimension(nx,ny,4), device :: Q3
    integer i, j, k
    integer :: offset = accuracy / 2
    real(8) R
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    do k = 1, 4
      call calc_R(dtdx,dtdy,E(i-offset,j-offset,k),E(i-offset+1,j-offset,k), &
              & F(i-offset,j-offset,k),F(i-offset,j-offset+1,k),R)
      Q3(i,j,k) = 0.75d0 * Q(i,j,k) + 0.25d0 * Q2(i,j,k) - 0.25d0 * R
    enddo 
  end subroutine calc_step2
  
  attributes(global) subroutine calc_step3(nx,ny,dtdx,dtdy,E,F,Q3,Q)
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: dtdx, dtdy
    real(8), intent(in), device :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(in), device :: F(nx-accuracy,ny-accuracy+1,4)
    real(8), intent(in), dimension(nx,ny,4), device :: Q3
    real(8), intent(inout), dimension(nx,ny,4), device :: Q
    integer i, j, k
    integer :: offset = accuracy / 2
    real(8) R
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    do k = 1, 4
      call calc_R(dtdx,dtdy,E(i-offset,j-offset,k),E(i-offset+1,j-offset,k), &
              & F(i-offset,j-offset,k),F(i-offset,j-offset+1,k),R)
      Q(i,j,k) = (Q(i,j,k) + 2.0d0 * Q3(i,j,k) - 2.d0 * R) / 3.d0
    enddo 
  end subroutine calc_step3
end module calc_steps

