module set
  use mod_globals, only : nx, ny, nz, dy, gamma
  implicit none
contains
  subroutine set_init(Q)
    real(8), intent(out), dimension(nx,ny,nz,5) :: Q
    integer j
    real(8) y, v
    real(8) :: delta = 0.1d0
    do j = 1, ny
      y = dy * dble(j - 1)
      v = 2.d0 * (y / delta)**(1.d0/7.d0)
      Q(:,j,:,1) = 1.4d0
      Q(:,j,:,2) = 1.4d0 * v
      Q(:,j,:,3) = 0.d0
      Q(:,j,:,4) = 0.d0
      Q(:,j,:,5) = 1.d0 / (gamma - 1.d0) + 0.5d0 * 1.4d0 * v**2
    enddo
  end subroutine set_init
  
  subroutine set_bc(id_accuracy,Q)
    integer(kind=2), intent(in), value :: id_accuracy
    real(8), intent(inout), device :: Q(nx,ny,nz,5)
    real(8), y, v, p_out
    real(8) :: delta = 0.1d0
    integer i, j, k
    !$cuf kernel do <<<*,*>>>
    do k = 1, nz
      do j = 2, ny-1
        ! outlet
        Q(nx,j,k,1) = Q(nx-1,j,k,1)
        Q(nx,j,k,2) = Q(nx-1,j,k,2)
        Q(nx,j,k,3) = Q(nx-1,j,k,3)
        Q(nx,j,k,4) = Q(nx-1,j,k,4)
        Q(nx,j,k,5) = Q(nx-1,j,k,5)
        ! inlet
        Q(1,j,k,1) = Q(nx,j,k,1)
        y = dy * dble(j - 1)
        v = 2.d0 * (y / delta)**(1.d0 / 7.d0)
        Q(1,j,k,2) = 0.5d0 * (Q(nx,j,k,2) + Q(nx,j,k,1) * v)
        Q(1,j,k,3) = Q(nx,j,k,3)
        Q(1,j,k,4) = Q(nx,j,k,4)
        p_out = (gamma - 1.d0) * (Q(nx,j,k,5) - 0.5d0 * &
        & (Q(nx,j,k,2)**2 + Q(nx,j,k,3)**2 + Q(nx,j,k,4)**2) / Q(nx,j,k,1))
        Q(1,j,k,5) = p_out / (gamma - 1.d0) + 0.5d0 * &
        & Q(nx,j,k,1) * (v**2 + (Q(nx,j,k,3) / Q(nx,j,k,1))**2 + (Q(nx,j,k,4) / Q(nx,j,k,1))**2)
      enddo
    enddo

    !$cuf kernel do <<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        ! bottom
        Q(i,2,k,1) = Q(i,3,k,1)
        Q(i,2,k,2) = 0.d0
        Q(i,2,k,3) = 0.d0
        Q(i,2,k,4) = 0.d0
        Q(i,2,k,5) = Q(i,3,k,5)
        ! imaginary
        Q(i,1,k,1) = Q(i,3,k,1)
        Q(i,1,k,2) = Q(i,3,k,2)
        Q(i,1,k,3) = -Q(i,3,k,3)
        Q(i,1,k,4) = Q(i,3,k,4)
        Q(i,1,k,5) = Q(i,3,k,5)
        ! top
        Q(i,ny,k,:) = Q(i,ny-1,k,:)
      enddo
    enddo

    !$cuf kernel do <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        ! cyclic
        Q(i,j,1,:) = Q(i,j,nz-1,:)
        Q(i,j,nz,:) = Q(i,j,2,:)
      enddo
    enddo
  end subroutine set_bc
end module set

