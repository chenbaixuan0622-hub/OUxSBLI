module calc_steps
  use mod_globals, only : accuracy
  implicit none
contains
  function R(dX,dY,dZ,E,F,G) result(ans)
    !$acc routine
    real(8), intent(in), value :: dX, dY, dZ
    real(8), intent(in), dimension(2,5), device :: E, F, G
    real(8) :: ans(5)
    ans(:) = dX * (-E(1,:) + E(2,:)) +dY * (-F(1,:) + F(2,:)) + dZ * (-G(1,:) + G(2,:))
  end function R

  subroutine calc_step1(nx,ny,nz,dX,dY,dZ,dt,E,F,G,Q,Q2)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: dX, dY, dZ, dt
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
          Q2(i,j,k,:) = Q(i,j,k,:)&
          -dt*R(dX,dY,dZ,E(i-offset:i-offset+1,j-offset,k-offset,:),&
          F(i-offset,j-offset:j-offset+1,k-offset,:),G(i-offset,j-offset,k-offset:k-offset+1,:))
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_step1
  
  subroutine calc_step2(nx,ny,nz,dX,dY,dZ,dt,E,F,G,Q,Q2,Q3)
    integer, intent(in), value :: nx, ny,nz
    real(8), intent(in), value :: dX, dY, dZ, dt
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
          Q3(i,j,k,:) = 0.75d0 * Q(i,k,j,:) + 0.25d0 * Q2(i,j,k,:)&
          -0.25d0*dt*R(dX,dY,dZ,E(i-offset:i-offset+1,j-offset,k-offset,:),&
          F(i-offset,j-offset:j-offset+1,k-offset,:),G(i-offset,j-offset,k-offset:k-offset+1,:))
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_step2
  
  subroutine calc_step3(nx,ny,nz,dX,dY,dZ,dt,E,F,G,Q3,Q)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: dX, dY, dZ, dt
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
          Q(i,j,k,:) = (Q(i,j,k,:) + 2.0d0 * Q3(i,j,k,:)&
          -2.0d0*dt*R(dX,dY,dZ,E(i-offset:i-offset+1,j-offset,k-offset,:),&
          F(i-offset,j-offset:j-offset+1,k-offset,:),G(i-offset,j-offset,k-offset:k-offset+1,:)))/3.d0
        enddo
      enddo
    enddo
    !$acc end kernels
  end subroutine calc_step3
end module calc_steps