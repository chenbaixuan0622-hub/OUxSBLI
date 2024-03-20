module calc_steps
  use cudafor
  use mod_globals, only : accuracy, offset, nx, ny, nz, dtdx, dtdy, dtdz
  implicit none
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

  subroutine calc_step1_2D(E,F,Q,Q2)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    real(8), intent(in), dimension(nx,ny,4), device :: Q
    real(8), intent(out), dimension(nx,ny,4), device :: Q2
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1+offset, ny-offset
        do i =1+offset, nx-offset
          R = dtdx * (-E(i-offset,j-offset,k) + E(i-offset+1,j-offset,k)) &
          & + dtdy * (-F(i-offset,j-offset,k) + F(i-offset,j-offset+1,k))
          Q2(i,j,k) = Q(i,j,k) - R
        enddo
      enddo
    enddo 
  end subroutine calc_step1_2D

  subroutine calc_step1_3D(E,F,G,Q,Q2)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q
    real(8), intent(out), dimension(nx,ny,nz,5), device :: Q2
    integer i, j, k, l
    real(8) R
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1+offset, nz-offset
        do j = 1+offset, ny-offset
          do i = 1+offset, nx-offset
            R = dtdx * (-E(i-offset,j-offset,k-offset,l) + E(i-offset+1,j-offset,k-offset,l)) &
            & + dtdy * (-F(i-offset,j-offset,k-offset,l) + F(i-offset,j-offset+1,k-offset,l)) &
            & + dtdz * (-G(i-offset,j-offset,k-offset,l) + G(i-offset,j-offset,k-offset+1,l))
            Q2(i,j,k,l) = Q(i,j,k,l) - R
          enddo
        enddo
      enddo
    enddo 
  end subroutine calc_step1_3D
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_step2_2D(E,F,Q,Q2,Q3)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    real(8), intent(in), dimension(nx,ny,4), device :: Q
    real(8), intent(in), dimension(nx,ny,4), device :: Q2
    real(8), intent(out), dimension(nx,ny,4), device :: Q3
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          R = dtdx * (-E(i-offset,j-offset,k) + E(i-offset+1,j-offset,k)) &
          & + dtdy * (-F(i-offset,j-offset,k) + F(i-offset,j-offset+1,k))
          Q3(i,j,k) = 0.75d0 * Q(i,j,k) + 0.25d0 * Q2(i,j,k) - 0.25d0 * R
        enddo
      enddo
    enddo 
  end subroutine calc_step2_2D
  
  subroutine calc_step2_3D(E,F,G,Q,Q2,Q3)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q2
    real(8), intent(out), dimension(nx,ny,nz,5), device :: Q3
    integer i, j, k, l
    real(8) R
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1+offset, nz-offset
        do j = 1+offset, ny-offset
          do i = 1+offset, nx-offset
            R = dtdx * (-E(i-offset,j-offset,k-offset,l) + E(i-offset+1,j-offset,k-offset,l)) &
            & + dtdy * (-F(i-offset,j-offset,k-offset,l) + F(i-offset,j-offset+1,k-offset,l)) &
            & + dtdz * (-G(i-offset,j-offset,k-offset,l) + G(i-offset,j-offset,k-offset+1,l))
            Q3(i,j,k,l) = 0.75d0 * Q(i,j,k,l) + 0.25d0 * Q2(i,j,k,l) - 0.25d0 * R
          enddo
        enddo
      enddo
    enddo
  end subroutine calc_step2_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_step3_2D(E,F,Q3,Q)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    real(8), intent(in), dimension(nx,ny,4), device :: Q3
    real(8), intent(inout), dimension(nx,ny,4), device :: Q
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          R = dtdx * (-E(i-offset,j-offset,k) + E(i-offset+1,j-offset,k)) &
          & + dtdy * (-F(i-offset,j-offset,k) + F(i-offset,j-offset+1,k))
          Q(i,j,k) = (Q(i,j,k) + 2.0d0 * Q3(i,j,k) - 2.d0 * R) / 3.d0
        enddo
      enddo
    enddo 
  end subroutine calc_step3_2D
  
  subroutine calc_step3_3D(E,F,G,Q3,Q)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q3
    real(8), intent(inout), dimension(nx,ny,nz,5), device :: Q
    integer i, j, k, l
    real(8) R
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1+offset, nz-offset
        do j = 1+offset, ny-offset
          do i = 1+offset, nx-offset
            R = dtdx * (-E(i-offset,j-offset,k-offset,l) + E(i-offset+1,j-offset,k-offset,l)) &
            & + dtdy * (-F(i-offset,j-offset,k-offset,l) + F(i-offset,j-offset+1,k-offset,l)) &
            & + dtdz * (-G(i-offset,j-offset,k-offset,l) + G(i-offset,j-offset,k-offset+1,l))
            Q(i,j,k,l) = (Q(i,j,k,l) + 2.0d0 * Q3(i,j,k,l) - 2.d0 * R) / 3.d0
          enddo
        enddo
      enddo
    enddo
  end subroutine calc_step3_3D
end module calc_steps

