module calc_steps
  use mod_globals, only : accuracy
  implicit none
contains
  attributes(global) subroutine calc_step1(nx,ny,nz,dX,dY,dZ,dt,E,F,G,Q,Q2)
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
    offset = accuracy / 2
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    Q2(i,j,k,:) = Q(i,j,k,:) - dt * ( &
    & -dX * (-E(i-offset,j-offset,k-offset,:) + E(i-offset+1,j-offset,k-offset,:)) &
    & -dY * (-F(i-offset,j-offset,k-offset,:) + F(i-offset,j-offset+1,k-offset,:)) &
    & -dZ * (-G(i-offset,j-offset,k-offset,:) + G(i-offset,j-offset,k-offset+1,:)))
  end subroutine calc_step1
  
  attributes(global) subroutine calc_step2(nx,ny,nz,dX,dY,dZ,dt,E,F,G,Q,Q2,Q3)
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
    offset = accuracy / 2
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    Q3(i,j,k,:) = 0.75d0 * Q(i,k,j,:) + 0.25d0 * Q2(i,j,k,:) - 0.25d0 * dt * ( &
    & -dX * (-E(i-offset,j-offset,k-offset,:) + E(i-offset+1,j-offset,k-offset,:)) &
    & -dY * (-F(i-offset,j-offset,k-offset,:) + F(i-offset,j-offset+1,k-offset,:)) &
    & -dZ * (-G(i-offset,j-offset,k-offset,:) + G(i-offset,j-offset,k-offset+1,:)))
  end subroutine calc_step2
  
  attributes(global) subroutine calc_step3(nx,ny,nz,dX,dY,dZ,dt,E,F,G,Q3,Q)
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
    offset = accuracy / 2
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    Q(i,j,k,:) = (Q(i,j,k,:) + 2.0d0 * Q3(i,j,k,:) - 2.d0 * dt * ( &
    & -dX * (-E(i-offset,j-offset,k-offset,:) + E(i-offset+1,j-offset,k-offset,:)) &
    & -dY * (-F(i-offset,j-offset,k-offset,:) + F(i-offset,j-offset+1,k-offset,:)) &
    & -dZ * (-G(i-offset,j-offset,k-offset,:) + G(i-offset,j-offset,k-offset+1,:)))) / 3.d0
  end subroutine calc_step3
end module calc_steps