module set
  use mod_globals, only : accuracy, offset, nx, ny, nz, dx, dy, dz, gamma, RHO0, L0, V0, p0
  implicit none

  interface set_bc
    module procedure set_bc2, set_bc4
  end interface
  
contains
  function linspace(x1, x2, n) result(x)
    real(8), intent(in) :: x1, x2
    integer, intent(in) :: n
    integer i
    real(8) x(n)
    x = x1 + (x1 + x2) * (/ (dble(i - 1) / dble(n - 1), i = 1, n) /)
  end function linspace

  subroutine set_grid(x,y,z)
    real(8), intent(out), dimension(nx,ny,nz) :: x, y, z
    integer i, j, k
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          x(i,j,k) = dble(i-1) * dx
          y(i,j,k) = dble(j-1) * dy
          z(i,j,k) = dble(k-1) * dz
        enddo
      enddo
    enddo
  end subroutine set_grid
  
  subroutine set_init(Q,Vin)
    real(8), intent(out), dimension(nx,ny,nz,5) :: Q
    real(8), intent(in), dimension(ny,2) :: Vin
    integer i, j, k
    real(8) :: pi = 2.d0 * acos(0.d0)
    real(8) x(nx-accuracy), y(ny-accuracy), z(nz-accuracy)
    x = linspace(0.d0, 2.d0 * pi, nx-accuracy)
    y = linspace(0.d0, 2.d0 * pi, ny-accuracy)
    z = linspace(0.d0, 2.d0 * pi, nz-accuracy)
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          ! rho
          Q(i,j,k,1) = RHO0
          ! rho u
          Q(i,j,k,2) = RHO0 * V0 * sin(x(i-offset)) * cos(y(j-offset)) * cos(z(k-offset))
          ! rho v
          Q(i,j,k,3) = - RHO0 * V0 * cos(x(i-offset)) * sin(y(j-offset)) * cos(z(k-offset))
          ! rho w
          Q(i,j,k,4) = 0.d0
          ! p / (gamma - 1) + 0.5 * (rhou ** 2 + rhov ** 2 ) / rho
          Q(i,j,k,5) = (p0+RHO0*(V0**2)*(cos(2.d0*x(i-offset)/L0)+cos(2.d0*y(j-offset)/L0))*(cos(2.d0*z(k-offset)/L0)+2.d0)/16.d0)&
          / (gamma - 1.d0) + 0.5d0 * (Q(i,j,k,2) ** 2 + Q(i,j,k,3) ** 2) / Q(i,j,k,1)
        enddo
      enddo
    enddo
    if (accuracy == 2) then
      call set_bc_init2(Q)
    else
      call set_bc_init4(Q)
    endif
  end subroutine set_init
  
  subroutine set_bc_init2(Q)
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer i, j, k
    do k = 2, nz-1
      do j = 2, ny-1
        Q(1,j,k,:) = Q(nx-1,j,k,:)
        Q(nx,j,k,:) = Q(2,j,k,:)
      enddo
    enddo

    do k = 2, nz-1
      do i = 2, nx-1
        Q(i,1,k,:) = Q(i,ny-1,k,:)
        Q(i,ny,k,:) = Q(i,2,k,:)
      enddo
    enddo

    do k = 2, nz-1
      Q(1,1,k,:) = Q(nx-1,ny-1,k,:)
      Q(nx,1,k,:) = Q(2,ny-1,k,:)
      Q(1,ny,k,:) = Q(nx-1,2,k,:)
      Q(nx,ny,k,:) = Q(2,2,k,:)
    enddo

    do j = 1, ny
      do i = 1, nx
        Q(i,j,1,:) = Q(i,j,nz-1,:)
        Q(i,j,nz,:) = Q(i,j,2,:)
      enddo
    enddo
  end subroutine set_bc_init2

  subroutine set_bc_init4(Q)
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer i, j, k

    do k = 3, nz-2
      do j = 3, ny-2
        Q(1:2,j,k,:) = Q(nx-3:nx-2,j,k,:)
        Q(nx-1:nx,j,k,:) = Q(3:4,j,k,:)
      enddo
    enddo

    do k = 3, nz-2
      do i = 3, nx-2
        Q(i,1:2,k,:) = Q(i,ny-3:ny-2,k,:)
        Q(i,ny-1:ny,k,:) = Q(i,3:4,k,:)
      enddo
    enddo

    do k = 3, nz-2
      Q(1:2,1:2,k,:) = Q(nx-3:nx-2,ny-3:ny-2,k,:)
      Q(nx-1:nx,1:2,k,:) = Q(3:4,ny-3:ny-2,k,:)
      Q(1:2,ny-1:ny,k,:) = Q(nx-3:nx-2,3:4,k,:)
      Q(nx-1:nx,ny-1:ny,k,:) = Q(3:4,3:4,k,:)
    enddo

    do j = 1, ny
      do i = 1, nx
        Q(i,j,1:2,:) = Q(i,j,nz-3:nz-2,:)
        Q(i,j,nz-1:nz,:) = Q(i,j,3:4,:)
      enddo
    enddo
  end subroutine set_bc_init4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_bc_mut(mut)
    real(8), intent(inout), device :: mut(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2) <<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        mut(1,j,k) = mut(nx-1,j,k)
        mut(nx,j,k) = mut(2,j,k)
      enddo
    enddo

    !$cuf kernel do(2) <<<*,*>>>
    do k = 2, nz-1
      do i = 2, nx-1
        mut(i,1,k) = mut(i,ny-1,k)
        mut(i,ny,k) = mut(i,2,k)
      enddo
    enddo

    !$cuf kernel do(1) <<<*,*>>>
    do k = 2, nz-1
      mut(1,1,k) = mut(nx-1,ny-1,k)
      mut(nx,1,k) = mut(2,ny-1,k)
      mut(1,ny,k) = mut(nx-1,2,k)
      mut(nx,ny,k) = mut(2,2,k)
    enddo

    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        mut(i,j,1) = mut(i,j,nz-1)
        mut(i,j,nz) = mut(i,j,2)
      enddo
    enddo
  end subroutine set_bc_mut

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_bc2(id_accuracy,Q,T)
    integer(kind=2), intent(in), value :: id_accuracy
    real(8), intent(inout), device :: Q(nx,ny,nz,5)
    real(8), intent(inout), device :: T(nx,ny,nz)
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do l = 1, 5
      do k = 2, nz-1
        do j = 2, ny-1
          Q(1,j,k,l) = Q(nx-1,j,k,l)
          Q(nx,j,k,l) = Q(2,j,k,l)
        enddo
      enddo
    enddo

    !$cuf kernel do(3) <<<*,*>>>
    do l = 1, 5
      do k = 2, nz-1
        do i = 2, nx-1
          Q(i,1,k,l) = Q(i,ny-1,k,l)
          Q(i,ny,k,l) = Q(i,2,k,l)
        enddo
      enddo
    enddo

    !$cuf kernel do(2) <<<*,*>>>
    do l = 1, 5
      do k = 2, nz-1
        Q(1,1,k,l) = Q(nx-1,ny-1,k,l)
        Q(nx,1,k,l) = Q(2,ny-1,k,l)
        Q(1,ny,k,l) = Q(nx-1,2,k,l)
        Q(nx,ny,k,l) = Q(2,2,k,l)
      enddo
    enddo

    !$cuf kernel do(3) <<<*,*>>>
    do l = 1, 5
      do j = 1, ny
        do i = 1, nx
          Q(i,j,1,l) = Q(i,j,nz-1,l)
          Q(i,j,nz,l) = Q(i,j,2,l)
        enddo
      enddo
    enddo
  end subroutine set_bc2

  subroutine set_bc4(id_accuracy,Q,T)
    integer(kind=4), intent(in), value :: id_accuracy
    real(8), intent(inout), device :: Q(nx,ny,nz,5)
    real(8), intent(inout), device :: T(nx,ny,nz)
    integer i, j, k, l
  
    !$cuf kernel do(3) <<<*,*>>>
    do l = 1, 5
      do k = 3, nz-2
        do j = 3, ny-2
          Q(1:2,j,k,l) = Q(nx-3:nx-2,j,k,l)
          Q(nx-1:nx,j,k,l) = Q(3:4,j,k,l)
        enddo
      enddo
    enddo
  
    !$cuf kernel do(3) <<<*,*>>>
    do l = 1, 5
      do k = 3, nz-2
        do i = 3, nx-2
          Q(i,1:2,k,l) = Q(i,ny-3:ny-2,k,l)
          Q(i,ny-1:ny,k,l) = Q(i,3:4,k,l)
        enddo
      enddo
    enddo
  
    !$cuf kernel do(2) <<<*,*>>>
    do l = 1, 5
      do k = 3, nz-2
        Q(1:2,1:2,k,l) = Q(nx-3:nx-2,ny-3:ny-2,k,l)
        Q(nx-1:nx,1:2,k,l) = Q(3:4,ny-3:ny-2,k,l)
        Q(1:2,ny-1:ny,k,l) = Q(nx-3:nx-2,3:4,k,l)
        Q(nx-1:nx,ny-1:ny,k,l) = Q(3:4,3:4,k,l)
      enddo
    enddo
  
    !$cuf kernel do(3) <<<*,*>>>
    do l = 1, 5
      do j = 1, ny
        do i = 1, nx
          Q(i,j,1:2,l) = Q(i,j,nz-3:nz-2,l)
          Q(i,j,nz-1:nz,l) = Q(i,j,3:4,l)
        enddo
      enddo
    enddo
  end subroutine set_bc4
end module set

