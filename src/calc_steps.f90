module calc_steps
  use cudafor
  use mod_globals, only : accuracy, offset, nx, ny, nz, dtdx, dtdy, dtdz
  implicit none
  ! 3rd-order TVD RungeKutta
  interface calc_step1
    module procedure calc_step1_2D, calc_step1_3D
  end interface
   
  interface calc_step2
    module procedure calc_step2_2D, calc_step2_3D
  end interface

  interface calc_step3
    module procedure calc_step3_2D, calc_step3_3D
  end interface

  ! 4th-order traditional RungeKutta
  interface calc_step
    module procedure calc_step_2D, calc_step_3D, calc_step_simple
  end interface

  interface calc_step4
    module procedure calc_step4_2D, calc_step4_3D
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
    real(8), intent(in), dimension(nx,ny,4), device :: Q, Q2
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
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q, Q2
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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine calc_step_2D(coef1,coef2,E,F,Rs,Q,Q2)
    real(8), intent(in) :: coef1, coef2
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,4), device :: Rs
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
          Rs(i-offset,j-offset,k) = Rs(i-offset,j-offset,k) + coef2 * R
          Q2(i,j,k) = Q(i,j,k) - coef1 * R
        enddo
      enddo
    enddo 
  end subroutine calc_step_2D

  subroutine calc_step_3D(coef1,coef2,E,F,G,Rs,Q,Q2)
    real(8), intent(in) :: coef1, coef2
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device :: Rs
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
            Rs(i-offset,j-offset,k-offset,l) = Rs(i-offset,j-offset,k-offset,l) + coef2 * R
            Q2(i,j,k,l) = Q(i,j,k,l) - coef1 * R
          enddo
        enddo
      enddo
    enddo
  end subroutine calc_step_3D
 
  subroutine calc_step_simple(coef,E,F,G,Q)
    real(8), intent(in) :: coef
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
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
            Q(i,j,k,l) = Q(i,j,k,l) - coef * R
          enddo
        enddo
      enddo
    enddo
  end subroutine calc_step_simple

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_step4_2D(E,F,Rs,Q)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,4), device :: Rs
    real(8), intent(inout), dimension(nx,ny,4), device :: Q
    integer i, j, k
    real(8) R
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, 4
      do j = 1+offset, ny-offset
        do i =1+offset, nx-offset
          R = dtdx * (-E(i-offset,j-offset,k) + E(i-offset+1,j-offset,k)) &
          & + dtdy * (-F(i-offset,j-offset,k) + F(i-offset,j-offset+1,k))
          Rs(i-offset,j-offset,k) = Rs(i-offset,j-offset,k) + R
          Q(i,j,k) = Q(i,j,k) - Rs(i-offset,j-offset,k) / 6.d0
          Rs(i-offset,j-offset,k) = 0.d0
        enddo
      enddo
    enddo 
  end subroutine calc_step4_2D

  subroutine calc_step4_3D(E,F,G,Rs,Q)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device :: Rs
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
            Rs(i-offset,j-offset,k-offset,l) = Rs(i-offset,j-offset,k-offset,l) + R
            Q(i,j,k,l) = Q(i,j,k,l) - Rs(i-offset,j-offset,k-offset,l) / 6.d0
            Rs(i-offset,j-offset,k-offset,l) = 0.d0
          enddo
        enddo
      enddo
    enddo
  end subroutine calc_step4_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_step5(E,F,G,Q,Q4,Q5,R4)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q, Q4
    real(8), intent(out), dimension(nx,ny,nz,5), device :: Q5
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device :: R4
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
            Q5(i,j,k,l) = (9.d0 * Q(i,j,k,l) + 6.0d0 * Q4(i,j,k,l) - R) / 15.d0
            R4(i-offset,j-offset,k-offset,l) = R
          enddo
        enddo
      enddo
    enddo
  end subroutine calc_step5

  subroutine calc_step10(E,F,G,R4,Q4,Q9,Q)
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device :: R4
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q4, Q9
    real(8), intent(inout), dimension(nx,ny,nz,5), device :: Q
    integer i, j, k, l
    real(8) R9
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1+offset, nz-offset
        do j = 1+offset, ny-offset
          do i = 1+offset, nx-offset
            R9 = dtdx * (-E(i-offset,j-offset,k-offset,l) + E(i-offset+1,j-offset,k-offset,l)) &
            & + dtdy * (-F(i-offset,j-offset,k-offset,l) + F(i-offset,j-offset+1,k-offset,l)) &
            & + dtdz * (-G(i-offset,j-offset,k-offset,l) + G(i-offset,j-offset,k-offset+1,l))
            Q(i,j,k,l) = (2.d0 * Q(i,j,k,l) + 18.0d0 * Q4(i,j,k,l) + 30.d0 * Q9(i,j,k,l) &
            & - 3.d0 * R4(i-offset,j-offset,k-offset,l) - 5.d0 * R9) / 50.d0
          enddo
        enddo
      enddo
    enddo
  end subroutine calc_step10
end module calc_steps

