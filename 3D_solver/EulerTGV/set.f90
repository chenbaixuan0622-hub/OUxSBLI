module set
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, gamma, R, RHO0, M0, V0, p0, T, dtn
  implicit none
contains
  function linspace(x1, x2, n) result(x)
    real(8), intent(in) :: x1, x2
    integer, intent(in) :: n
    integer i
    real(8) x(n)
    x = x1 + (x1 + x2) * (/ (dble(i - 1) / dble(n - 1), i = 1, n) /)
  end function linspace

  subroutine set_grid(nx,ny,nz,x,y,z,dx,dy,dz)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx), dy(ny), dz(nz)
    integer i, j, k
    real(8) dx1, dy1, dz1
    dx1 = Lx / dble(nx-4)
    dy1 = Ly / dble(ny-4)
    dz1 = Lz / dble(nz-4)
    x(1) = 0.d0
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo
    y(1) = 0.d0
    do j = 1, ny-1
      dy(j) = dy1
      y(j+1) = y(j) + dy(j)
    enddo
    z(1) = 0.d0
    do k = 1, nz-1
      dz(k) = dz1
      z(k+1) = z(k) + dz(k)
    enddo
  end subroutine set_grid
  
  subroutine set_init(nx,ny,nz,xs,ys,zs,Q)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: xs(nx), ys(ny), zs(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer i, j, k
    integer :: accuracy = 4
    integer :: offset   = 2
    real(8) x(nx-accuracy), y(ny-accuracy), z(nz-accuracy)
    real(8) dx, dy, dz
    dx = Lx / dble(nx-4)
    dy = Ly / dble(ny-4)
    dz = Lz / dble(nz-4)
    x = linspace(0.5d0 * dx, Lx - 0.5d0 * dx, nx-accuracy)
    y = linspace(0.5d0 * dy, Ly - 0.5d0 * dy, ny-accuracy)
    z = linspace(0.5d0 * dz, Lz - 0.5d0 * dx, nz-accuracy)
    ! 2nd-order accuracy : offset = 1
    ! 4th-order accuracy : offset = 2
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          ! rho
          Q(i,j,k,1) =  RHO0
          ! rho u
          Q(i,j,k,2) =  RHO0 * M0 * sin(x(i-offset)) * cos(y(j-offset)) * cos(z(k-offset))
          ! rho v
          Q(i,j,k,3) = -RHO0 * M0 * cos(x(i-offset)) * sin(y(j-offset)) * cos(z(k-offset))
          ! rho w0
          Q(i,j,k,4) = 0.d0
          ! p / (gamma - 1) + 0.5 * (rhou ** 2 + rhov ** 2 ) / rho
          Q(i,j,k,5) = (1.d0/gamma+0.0625d0*RHO0*(M0**2)*(cos(2.d0*x(i-offset))+cos(2.d0*y(j-offset)))*(cos(2.d0*z(k-offset))+2.d0))&
                       / (gamma - 1.d0) + 0.5d0 * (Q(i,j,k,2)**2 + Q(i,j,k,3)**2 + Q(i,j,k,4)**2) / Q(i,j,k,1)
    enddo;enddo;enddo
    !call set_bc_init2(Q)
    call set_bc_init4(Q)
  end subroutine set_init
  
  subroutine set_bc_init2(Q)
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer i, j, k
    do k = 2, nz-1
      do j = 2, ny-1
        Q(1,j,k,:) = Q(nx-1,j,k,:)
        Q(nx,j,k,:) = Q(2,j,k,:)
    enddo;enddo

    do k = 2, nz-1
      do i = 2, nx-1
        Q(i,1,k,:) = Q(i,ny-1,k,:)
        Q(i,ny,k,:) = Q(i,2,k,:)
    enddo;enddo

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
    enddo;enddo
  end subroutine set_bc_init2

  subroutine set_bc_init4(Q)
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer i, j, k

    do k = 3, nz-2
      do j = 3, ny-2
        Q(1:2,j,k,:) = Q(nx-3:nx-2,j,k,:)
        Q(nx-1:nx,j,k,:) = Q(3:4,j,k,:)
    enddo;enddo

    do k = 3, nz-2
      do i = 3, nx-2
        Q(i,1:2,k,:) = Q(i,ny-3:ny-2,k,:)
        Q(i,ny-1:ny,k,:) = Q(i,3:4,k,:)
    enddo;enddo

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
    enddo;enddo
  end subroutine set_bc_init4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_bc2(nx,ny,nz,Jacobian,Q)
    integer, intent(in), value      :: nx, ny, nz
    real(8), intent(in), device     :: Jacobian(nx,ny)
    real(8), intent(inout), device  :: Q(nx,ny,nz,5)
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do l = 1, 5
      do k = 2, nz-1
        do j = 2, ny-1
          Q(1,j,k,l) = Q(nx-1,j,k,l)
          Q(nx,j,k,l) = Q(2,j,k,l)
    enddo;enddo;enddo

    !$cuf kernel do(3) <<<*,*>>>
    do l = 1, 5
      do k = 2, nz-1
        do i = 2, nx-1
          Q(i,1,k,l) = Q(i,ny-1,k,l)
          Q(i,ny,k,l) = Q(i,2,k,l)
    enddo;enddo;enddo

    !$cuf kernel do(2) <<<*,*>>>
    do l = 1, 5
      do k = 2, nz-1
        Q(1,1,k,l) = Q(nx-1,ny-1,k,l)
        Q(nx,1,k,l) = Q(2,ny-1,k,l)
        Q(1,ny,k,l) = Q(nx-1,2,k,l)
        Q(nx,ny,k,l) = Q(2,2,k,l)
    enddo;enddo

    !$cuf kernel do(3) <<<*,*>>>
    do l = 1, 5
      do j = 1, ny
        do i = 1, nx
          Q(i,j,1,l) = Q(i,j,nz-1,l)
          Q(i,j,nz,l) = Q(i,j,2,l)
    enddo;enddo;enddo
  end subroutine set_bc2

  subroutine set_bc(nx,ny,nz,Jacobian,Q)
    integer, intent(in), value      :: nx, ny, nz
    real(8), intent(in), device     :: Jacobian(nx,ny)
    real(8), intent(inout), device  :: Q(nx,ny,nz,5)
    integer i, j, k, l
  
    !$cuf kernel do(3) <<<*,*>>>
    do l = 1, 5
      do k = 3, nz-2
        do j = 3, ny-2
          Q(1,j,k,l) = Q(nx-3,j,k,l)
          Q(2,j,k,l) = Q(nx-2,j,k,l)
          Q(nx-1,j,k,l) = Q(3,j,k,l)
          Q(nx,j,k,l) = Q(4,j,k,l)
    enddo;enddo;enddo
  
    !$cuf kernel do(3) <<<*,*>>>
    do l = 1, 5
      do k = 3, nz-2
        do i = 3, nx-2
          Q(i,1,k,l) = Q(i,ny-3,k,l)
          Q(i,2,k,l) = Q(i,ny-2,k,l)
          Q(i,ny-1,k,l) = Q(i,3,k,l)
          Q(i,ny,k,l) = Q(i,4,k,l)
    enddo;enddo;enddo
  
    !$cuf kernel do(2) <<<*,*>>>
    do l = 1, 5
      do k = 3, nz-2
        Q(1,1,k,l) = Q(nx-3,ny-3,k,l)
        Q(1,2,k,l) = Q(nx-3,ny-2,k,l)
        Q(2,1,k,l) = Q(nx-2,ny-3,k,l)
        Q(2,2,k,l) = Q(nx-2,ny-2,k,l)

        Q(nx-1,1,k,l) = Q(3,ny-3,k,l)
        Q(nx-1,2,k,l) = Q(3,ny-2,k,l)
        Q(nx,1,k,l)   = Q(4,ny-3,k,l)
        Q(nx,2,k,l)   = Q(4,ny-2,k,l)

        Q(1,ny-1,k,l) = Q(nx-3,3,k,l)
        Q(1,ny,k,l)   = Q(nx-3,4,k,l)
        Q(2,ny-1,k,l) = Q(nx-2,3,k,l)
        Q(2,ny,k,l)   = Q(nx-2,4,k,l)

        Q(nx-1,ny-1,k,l) = Q(3,3,k,l)
        Q(nx-1,ny,k,l)   = Q(3,4,k,l)
        Q(nx,ny-1,k,l)   = Q(4,3,k,l)
        Q(nx,ny,k,l)     = Q(4,4,k,l)
    enddo;enddo
  
    !$cuf kernel do(3) <<<*,*>>>
    do l = 1, 5
      do j = 1, ny
        do i = 1, nx
          Q(i,j,1,l) = Q(i,j,nz-3,l)
          Q(i,j,2,l) = Q(i,j,nz-2,l)
          Q(i,j,nz-1,l) = Q(i,j,3,l)
          Q(i,j,nz,l) = Q(i,j,4,l)
    enddo;enddo;enddo
  end subroutine set_bc

  subroutine set_bc_mut(nx,ny,nz,mut)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2) <<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        mut(1,j,k) = mut(nx-1,j,k)
        mut(nx,j,k) = mut(2,j,k)
    enddo;enddo

    !$cuf kernel do(2) <<<*,*>>>
    do k = 2, nz-1
      do i = 2, nx-1
        mut(i,1,k) = mut(i,ny-1,k)
        mut(i,ny,k) = mut(i,2,k)
    enddo;enddo

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
    enddo;enddo
  end subroutine set_bc_mut
end module set

