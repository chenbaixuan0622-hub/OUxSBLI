module calc_steps
  use mod_globals, only : accuracy
  implicit none
  interface calc_step1
    module procedure calc_step1_Euler, calc_step1_NS
  end interface

  interface calc_step2
    module procedure calc_step2_Euler, calc_step2_NS
  end interface

  interface calc_step3
    module procedure calc_step3_Euler, calc_step3_NS
  end interface
contains
  subroutine calc_step1_Euler(nx,ny,dtdx,dtdy,E,F,Q,Q2)
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: dtdx, dtdy
    real(8), intent(in), device :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(in), device :: F(nx-accuracy,ny-accuracy+1,4)
    real(8), intent(in), dimension(nx,ny,4), device :: Q
    real(8), intent(out), dimension(nx,ny,4), device :: Q2
    integer i, j, offset
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    !$acc kernels deviceptr(E,F,Q,Q2)
    offset = accuracy / 2
    !$acc loop collapse(2)
    do j = 1+offset, ny-offset
      do i = 1+offset, nx-offset
        Q2(i,j,:) = Q(i,j,:) &
                & - dtdx * (-E(i-offset,j-offset,:) + E(i-offset+1,j-offset,:)) &
                & - dtdy * (-F(i-offset,j-offset,:) + F(i-offset,j-offset+1,:)) 
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_step1_Euler

  subroutine calc_step1_NS(nx,ny,dtdx,dtdy,E,F,Ev,Fv,Q,Q2)
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: dtdx, dtdy
    real(8), intent(in), device :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(in), device :: F(nx-accuracy,ny-accuracy+1,4)
    real(8), intent(in), device :: Ev(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(in), device :: Fv(nx-accuracy,ny-accuracy+1,4)
    real(8), intent(in), dimension(nx,ny,4), device :: Q
    real(8), intent(out), dimension(nx,ny,4), device :: Q2
    integer i, j, offset
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    !$acc kernels deviceptr(E,F,Q,Q2)
    offset = accuracy / 2
    !$acc loop collapse(2)
    do j = 1+offset, ny-offset
      do i = 1+offset, nx-offset
        Q2(i,j,:) = Q(i,j,:) &
                & - dtdx * (-E(i-offset,j-offset,:) + E(i-offset+1,j-offset,:) & 
                & + Ev(i-offset,j-offset,:) - Ev(i-offset+1,j-offset,:)) &
                & - dtdy * (-F(i-offset,j-offset,:) + F(i-offset,j-offset+1,:) &
                & + Fv(i-offset,j-offset,:) - Fv(i-offset,j-offset+1,:))
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_step1_NS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine calc_step2_Euler(nx,ny,dtdx,dtdy,E,F,Q,Q2,Q3)
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: dtdx, dtdy
    real(8), intent(in), device :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(in), device :: F(nx-accuracy,ny-accuracy+1,4)
    real(8), intent(in), dimension(nx,ny,4), device :: Q
    real(8), intent(in), dimension(nx,ny,4), device :: Q2
    real(8), intent(out), dimension(nx,ny,4), device :: Q3
    integer i, j, offset
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    !$acc kernels deviceptr(E,F,Q,Q2,Q3)
    offset = accuracy / 2
    !$acc loop collapse(2)
    do j = 1+offset, ny-offset
      do i = 1+offset, nx-offset
        Q3(i,j,:) = 0.75d0 * Q(i,j,:) + 0.25d0 * Q2(i,j,:) - 0.25d0 * ( &
                & + dtdx * (-E(i-offset,j-offset,:) + E(i-offset+1,j-offset,:)) &
                & + dtdy * (-F(i-offset,j-offset,:) + F(i-offset,j-offset+1,:))) 
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_step2_Euler

  subroutine calc_step2_NS(nx,ny,dtdx,dtdy,E,F,Ev,Fv,Q,Q2,Q3)
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: dtdx, dtdy
    real(8), intent(in), device :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(in), device :: F(nx-accuracy,ny-accuracy+1,4)
    real(8), intent(in), device :: Ev(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(in), device :: Fv(nx-accuracy,ny-accuracy+1,4)
    real(8), intent(in), dimension(nx,ny,4), device :: Q
    real(8), intent(in), dimension(nx,ny,4), device :: Q2
    real(8), intent(out), dimension(nx,ny,4), device :: Q3
    integer i, j, offset
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    !$acc kernels deviceptr(E,F,Q,Q2,Q3)
    offset = accuracy / 2
    !$acc loop collapse(2)
    do j = 1+offset, ny-offset
      do i = 1+offset, nx-offset
        Q3(i,j,:) = 0.75d0 * Q(i,j,:) + 0.25d0 * Q2(i,j,:) - 0.25d0 * ( &
                & + dtdx * (-E(i-offset,j-offset,:) + E(i-offset+1,j-offset,:) &
                & + Ev(i-offset,j-offset,:) - Ev(i-offset+1,j-offset,:)) &
                & + dtdy * (-F(i-offset,j-offset,:) + F(i-offset,j-offset+1,:) &
                & + Fv(i-offset,j-offset,:) - Fv(i-offset,j-offset+1,:))) 
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_step2_NS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine calc_step3_Euler(nx,ny,dtdx,dtdy,E,F,Q3,Q)
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: dtdx, dtdy
    real(8), intent(in), device :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(in), device :: F(nx-accuracy,ny-accuracy+1,4)
    real(8), intent(in), dimension(nx,ny,4), device :: Q3
    real(8), intent(inout), dimension(nx,ny,4), device :: Q
    integer i, j, offset
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    !$acc kernels deviceptr(E,F,Q3,Q)
    offset = accuracy / 2
    !$acc loop collapse(2)
    do j = 1+offset, ny-offset
      do i = 1+offset, nx-offset
        Q(i,j,:) = (Q(i,j,:) + 2.0d0 * Q3(i,j,:) - 2.d0 * ( &
                & + dtdx * (-E(i-offset,j-offset,:) + E(i-offset+1,j-offset,:)) &
                & + dtdy * (-F(i-offset,j-offset,:) + F(i-offset,j-offset+1,:)))) / 3.d0 
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_step3_Euler

  subroutine calc_step3_NS(nx,ny,dtdx,dtdy,E,F,Ev,Fv,Q3,Q)
    integer, intent(in), value :: nx, ny
    real(8), intent(in), value :: dtdx, dtdy
    real(8), intent(in), device :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(in), device :: F(nx-accuracy,ny-accuracy+1,4)
    real(8), intent(in), device :: Ev(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(in), device :: Fv(nx-accuracy,ny-accuracy+1,4)
    real(8), intent(in), dimension(nx,ny,4), device :: Q3
    real(8), intent(inout), dimension(nx,ny,4), device :: Q
    integer i, j, offset
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    !$acc kernels deviceptr(E,F,Q3,Q)
    offset = accuracy / 2
    !$acc loop collapse(2)
    do j = 1+offset, ny-offset
      do i = 1+offset, nx-offset
        Q(i,j,:) = (Q(i,j,:) + 2.0d0 * Q3(i,j,:) - 2.d0 * ( &
                & + dtdx * (-E(i-offset,j-offset,:) + E(i-offset+1,j-offset,:) & 
                & + Ev(i-offset,j-offset,:) - Ev(i-offset+1,j-offset,:)) &
                & + dtdy * (-F(i-offset,j-offset,:) + F(i-offset,j-offset+1,:) &
                & + Fv(i-offset,j-offset,:) - Fv(i-offset,j-offset+1,:)))) / 3.d0
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_step3_NS
end module calc_steps

