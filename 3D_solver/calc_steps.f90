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
  subroutine calc_step1_Euler(nx,ny,nz,dtdx,dtdy,dtdz,E,F,G,Q,Q2)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: dtdx, dtdy, dtdz
    real(8), intent(in), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(in), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(in), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q
    real(8), intent(out), dimension(nx,ny,nz,5), device :: Q2
    integer i, j, k, offset
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    !$acc kernels deviceptr(E,F,G,Q,Q2)
    offset = accuracy / 2
    !$acc loop collapse(3)
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          Q2(i,j,k,:) = Q(i,j,k,:) &
          & - dtdx * (-E(i-offset,j-offset,k-offset,:) + E(i-offset+1,j-offset,k-offset,:)) &
          & - dtdy * (-F(i-offset,j-offset,k-offset,:) + F(i-offset,j-offset+1,k-offset,:)) &
          & - dtdz * (-G(i-offset,j-offset,k-offset,:) + G(i-offset,j-offset,k-offset+1,:))
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_step1_Euler

  subroutine calc_step1_NS(nx,ny,nz,dtdx,dtdy,dtdz,E,F,G,Ev,Fv,Gv,Q,Q2)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: dtdx, dtdy, dtdz
    real(8), intent(in), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(in), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(in), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), intent(in), device :: Ev(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(in), device :: Fv(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(in), device :: Gv(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q
    real(8), intent(out), dimension(nx,ny,nz,5), device :: Q2
    integer i, j, k, offset
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    !$acc kernels deviceptr(E,F,G,Q,Q2)
    offset = accuracy / 2
    !$acc loop collapse(3)
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          Q2(i,j,k,:) = Q(i,j,k,:) &
          & - dtdx * (-E(i-offset,j-offset,k-offset,:) + E(i-offset+1,j-offset,k-offset,:) & 
          & + Ev(i-offset,j-offset,k-offset,:) - Ev(i-offset+1,j-offset,k-offset,:)) &
          & - dtdy * (-F(i-offset,j-offset,k-offset,:) + F(i-offset,j-offset+1,k-offset,:) &
          & + Fv(i-offset,j-offset,k-offset,:) - Fv(i-offset,j-offset+1,k-offset,:)) &
          & - dtdz * (-G(i-offset,j-offset,k-offset,:) + G(i-offset,j-offset,k-offset+1,:) &
          & + Gv(i-offset,j-offset,k-offset,:) - Gv(i-offset,j-offset,k-offset+1,:))
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_step1_NS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine calc_step2_Euler(nx,ny,nz,dtdx,dtdy,dtdz,E,F,G,Q,Q2,Q3)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: dtdx, dtdy, dtdz
    real(8), intent(in), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(in), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(in), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q2
    real(8), intent(out), dimension(nx,ny,nz,5), device :: Q3
    integer i, j, k, offset
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    !$acc kernels deviceptr(E,F,G,Q,Q2,Q3)
    offset = accuracy / 2
    !$acc loop collapse(3)
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          Q3(i,j,k,:) = 0.75d0 * Q(i,j,k,:) + 0.25d0 * Q2(i,j,k,:) - 0.25d0 * ( &
          & + dtdx * (-E(i-offset,j-offset,k-offset,:) + E(i-offset+1,j-offset,k-offset,:)) &
          & + dtdy * (-F(i-offset,j-offset,k-offset,:) + F(i-offset,j-offset+1,k-offset,:)) &
          & + dtdz * (-G(i-offset,j-offset,k-offset,:) + G(i-offset,j-offset,k-offset+1,:)))
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_step2_Euler

  subroutine calc_step2_NS(nx,ny,nz,dtdx,dtdy,dtdz,E,F,G,Ev,Fv,Gv,Q,Q2,Q3)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: dtdx, dtdy, dtdz
    real(8), intent(in), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(in), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(in), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), intent(in), device :: Ev(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(in), device :: Fv(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(in), device :: Gv(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q2
    real(8), intent(out), dimension(nx,ny,nz,5), device :: Q3
    integer i, j, k, offset
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    !$acc kernels deviceptr(E,F,G,Q,Q2,Q3)
    offset = accuracy / 2
    !$acc loop collapse(3)
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          Q3(i,j,k,:) = 0.75d0 * Q(i,j,k,:) + 0.25d0 * Q2(i,j,k,:) - 0.25d0 * ( &
          & + dtdx * (-E(i-offset,j-offset,k-offset,:) + E(i-offset+1,j-offset,k-offset,:) &
          & + Ev(i-offset,j-offset,k-offset,:) - Ev(i-offset+1,j-offset,k-offset,:)) &
          & + dtdy * (-F(i-offset,j-offset,k-offset,:) + F(i-offset,j-offset+1,k-offset,:) &
          & + Fv(i-offset,j-offset,k-offset,:) - Fv(i-offset,j-offset+1,k-offset,:)) &
          & + dtdz * (-G(i-offset,j-offset,k-offset,:) + G(i-offset,j-offset,k-offset+1,:) &
          & + Gv(i-offset,j-offset,k-offset,:) - Gv(i-offset,j-offset,k-offset+1,:)))
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_step2_NS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine calc_step3_Euler(nx,ny,nz,dtdx,dtdy,dtdz,E,F,G,Q3,Q)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: dtdx, dtdy, dtdz
    real(8), intent(in), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(in), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(in), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q3
    real(8), intent(inout), dimension(nx,ny,nz,5), device :: Q
    integer i, j, k, offset
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    !$acc kernels deviceptr(E,F,G,Q3,Q)
    offset = accuracy / 2
    !$acc loop collapse(3)
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          Q(i,j,k,:) = (Q(i,j,k,:) + 2.0d0 * Q3(i,j,k,:) - 2.d0 * ( &
          & + dtdx * (-E(i-offset,j-offset,k-offset,:) + E(i-offset+1,j-offset,k-offset,:)) &
          & + dtdy * (-F(i-offset,j-offset,k-offset,:) + F(i-offset,j-offset+1,k-offset,:)) &
          & + dtdz * (-G(i-offset,j-offset,k-offset,:) + G(i-offset,j-offset,k-offset+1,:)))) / 3.d0
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_step3_Euler

  subroutine calc_step3_NS(nx,ny,nz,dtdx,dtdy,dtdz,E,F,G,Ev,Fv,Gv,Q3,Q)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: dtdx, dtdy, dtdz
    real(8), intent(in), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(in), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(in), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), intent(in), device :: Ev(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(in), device :: Fv(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(in), device :: Gv(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q3
    real(8), intent(inout), dimension(nx,ny,nz,5), device :: Q
    integer i, j, k, offset
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    !$acc kernels deviceptr(E,F,G,Q3,Q)
    offset = accuracy / 2
    !$acc loop collapse(3)
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          Q(i,j,k,:) = (Q(i,j,k,:) + 2.0d0 * Q3(i,j,k,:) - 2.d0 * ( &
          & + dtdx * (-E(i-offset,j-offset,k-offset,:) + E(i-offset+1,j-offset,k-offset,:) & 
          & + Ev(i-offset,j-offset,k-offset,:) - Ev(i-offset+1,j-offset,k-offset,:)) &
          & + dtdy * (-F(i-offset,j-offset,k-offset,:) + F(i-offset,j-offset+1,k-offset,:) &
          & + Fv(i-offset,j-offset,k-offset,:) - Fv(i-offset,j-offset+1,k-offset,:)) &
          & + dtdz * (-G(i-offset,j-offset,k-offset,:) + G(i-offset,j-offset,k-offset+1,:) &
          & + Gv(i-offset,j-offset,k-offset,:) - Gv(i-offset,j-offset,k-offset+1,:)))) / 3.d0
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_step3_NS
end module calc_steps

