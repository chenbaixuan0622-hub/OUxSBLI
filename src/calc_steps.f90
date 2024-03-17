module calc_steps
  use mod_globals, only : accuracy, offset, nx, ny, nz, dtdx, dtdy, dtdz
  implicit none
  interface calc_R
    module procedure calc_R_2D, calc_R_3D
  end interface

  interface calc_step1
    module procedure calc_step1_2D, calc_step1_3D
  end interface

   
  interface calc_step2
    module procedure calc_step2_2D, calc_step2_3D
  end interface

  interface calc_step3
    module procedure calc_step3_2D, calc_step3_3D
  end interface

contains
  attributes(device) subroutine calc_R_2D(E1,E2,F1,F2,R)
    real(8), intent(in), value :: E1, E2, F1, F2
    real(8), intent(out) :: R
    R = dtdx * (-E1 + E2) + dtdy * (-F1 + F2)
  end subroutine calc_R_2D
  
  attributes(device) subroutine calc_R_3D(E1,E2,F1,F2,G1,G2,R)
    real(8), intent(in), value :: E1, E2, F1, F2, G1, G2
    real(8), intent(out) :: R
    R = dtdx * (-E1 + E2) + dtdy * (-F1 + F2) + dtdz * (-G1 + G2)
  end subroutine calc_R_3D

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(global) subroutine calc_step1_2D(E,F,Q,Q2)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    real(8), intent(in), dimension(nx,ny,4), device :: Q
    real(8), intent(out), dimension(nx,ny,4), device :: Q2
    integer i, j, k
    real(8) R
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    do k = 1, 4
      call calc_R(E(i-offset,j-offset,k),E(i-offset+1,j-offset,k), &
              & F(i-offset,j-offset,k),F(i-offset,j-offset+1,k),R)
      Q2(i,j,k) = Q(i,j,k) - R
    enddo 
  end subroutine calc_step1_2D

  attributes(global) subroutine calc_step1_3D(E,F,G,Q,Q2)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q
    real(8), intent(out), dimension(nx,ny,nz,5), device :: Q2
    integer i, j, k, l
    real(8) R
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    do l = 1, 5
      call calc_R(E(i-offset,j-offset,k-offset,l),E(i-offset+1,j-offset,k-offset,l), &
              & F(i-offset,j-offset,k-offset,l),F(i-offset,j-offset+1,k-offset,l), &
              & G(i-offset,j-offset,k-offset,l),G(i-offset,j-offset,k-offset+1,l),R)
      Q2(i,j,k,l) = Q(i,j,k,l) - R
    enddo 
  end subroutine calc_step1_3D
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(global) subroutine calc_step2_2D(E,F,Q,Q2,Q3)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    real(8), intent(in), dimension(nx,ny,4), device :: Q
    real(8), intent(in), dimension(nx,ny,4), device :: Q2
    real(8), intent(out), dimension(nx,ny,4), device :: Q3
    integer i, j, k
    real(8) R
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    do k = 1, 4
      call calc_R(E(i-offset,j-offset,k),E(i-offset+1,j-offset,k), &
              & F(i-offset,j-offset,k),F(i-offset,j-offset+1,k),R)
      Q3(i,j,k) = 0.75d0 * Q(i,j,k) + 0.25d0 * Q2(i,j,k) - 0.25d0 * R
    enddo 
  end subroutine calc_step2_2D
  
  attributes(global) subroutine calc_step2_3D(E,F,G,Q,Q2,Q3)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q2
    real(8), intent(out), dimension(nx,ny,nz,5), device :: Q3
    integer i, j, k, l
    real(8) R
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    do l = 1, 5
      call calc_R(E(i-offset,j-offset,k-offset,l),E(i-offset+1,j-offset,k-offset,l), &
              & F(i-offset,j-offset,k-offset,l),F(i-offset,j-offset+1,k-offset,l), &
              & G(i-offset,j-offset,k-offset,l),G(i-offset,j-offset,k-offset+1,l),R)
      Q3(i,j,k,l) = 0.75d0 * Q(i,j,k,l) + 0.25d0 * Q2(i,j,k,l) - 0.25d0 * R
    enddo
  end subroutine calc_step2_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(global) subroutine calc_step3_2D(E,F,Q3,Q)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    real(8), intent(in), dimension(nx,ny,4), device :: Q3
    real(8), intent(inout), dimension(nx,ny,4), device :: Q
    integer i, j, k
    real(8) R
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    do k = 1, 4
      call calc_R(E(i-offset,j-offset,k),E(i-offset+1,j-offset,k), &
              & F(i-offset,j-offset,k),F(i-offset,j-offset+1,k),R)
      Q(i,j,k) = (Q(i,j,k) + 2.0d0 * Q3(i,j,k) - 2.d0 * R) / 3.d0
    enddo 
  end subroutine calc_step3_2D
  
  attributes(global) subroutine calc_step3_3D(E,F,G,Q3,Q)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q3
    real(8), intent(inout), dimension(nx,ny,nz,5), device :: Q
    integer i, j, k, l
    real(8) R
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    do l = 1, 5
      call calc_R(E(i-offset,j-offset,k-offset,l),E(i-offset+1,j-offset,k-offset,l), &
              & F(i-offset,j-offset,k-offset,l),F(i-offset,j-offset+1,k-offset,l), &
              & G(i-offset,j-offset,k-offset,l),G(i-offset,j-offset,k-offset+1,l),R)
      Q(i,j,k,l) = (Q(i,j,k,l) + 2.0d0 * Q3(i,j,k,l) - 2.d0 * R) / 3.d0
    enddo
  end subroutine calc_step3_3D
end module calc_steps

