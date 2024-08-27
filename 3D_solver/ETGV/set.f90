module set
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, gamma, R, RHO0, M0, V0, p0, T, dtn
  implicit none
contains
  subroutine set_grid(nx,ny,nz,xc,yc,zc,dx,dy,dz)
    use mod_globals, only : id_accuracy
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(out) :: xc(nx), yc(ny), zc(nz), dx(nx), dy(ny), dz(nz)
    real(8) dx1, dy1, dz1, x(nx+1), y(ny+1), z(nz+1)
    integer i, j, k
    if (kind(id_accuracy) == 2) then
      dx1 = Lx / dble(nx-2)
      dy1 = Ly / dble(ny-2)
      dz1 = Lz / dble(nz-2)
      dx(:) = dx1
      dy(:) = dy1
      dz(:) = dz1
    
      ! x direction
      do i = 2, nx
        x(i) = dx1 * dble(i-2)
      enddo
      x(1)    = x(2)  - dx1
      x(nx+1) = x(nx) + dx1
    
      ! y direction
      do j = 2, ny
        y(j) = dy1 * dble(j-2)
      enddo
      y(1)    = y(2)  - dy1
      y(ny+1) = y(ny) + dy1

      ! z direction
      do k = 2, nz
        z(k) = dz1 * dble(k-2)
      enddo
      z(1)    = z(2)  - dz1
      z(nz+1) = z(nz) + dz1
    elseif (kind(id_accuracy) == 4) then
      dx1 = Lx / dble(nx-4)
      dy1 = Ly / dble(ny-4)
      dz1 = Lz / dble(nz-4)
      dx(:) = dx1
      dy(:) = dy1
      dz(:) = dz1
    
      ! x direction
      do i = 3, nx-1
        x(i) = dx1 * dble(i-3)
      enddo
      x(1)    = x(3)    - 2.d0 * dx1
      x(2)    = x(3)    - dx1
      x(nx)   = x(nx-1) + dx1
      x(nx+1) = x(nx-1) + 2.d0 * dx1
    
      ! y direction
      do j = 3, ny-1
        y(j) = dy1 * dble(j-3)
      enddo
      y(1)    = y(3)    - 2.d0 * dy1
      y(2)    = y(3)    - dy1
      y(ny)   = y(ny-1) + dy1
      y(ny+1) = y(ny-1) + 2.d0 * dy1

      ! z direction
      do k = 3, nz-1
        z(k) = dz1 * dble(k-3)
      enddo
      z(1)    = z(3)    - 2.d0 * dz1
      z(2)    = z(3)    - dz1
      z(nz)   = z(nz-1) + dz1
      z(nz+1) = z(nz-1) + 2.d0 * dz1
    elseif (kind(id_accuracy) == 8) then
      dx1 = Lx / dble(nx-6)
      dy1 = Ly / dble(ny-6)
      dz1 = Lz / dble(nz-6)
      dx(:) = dx1
      dy(:) = dy1
      dz(:) = dz1
      ! x direction
      do i = 4, nx-2
        x(i) = dx1 * dble(i-4)
      enddo
      x(1)    = x(4)    - 3.d0 * dx1
      x(2)    = x(4)    - 2.d0 * dx1
      x(3)    = x(4)    - dx1
      x(nx-1) = x(nx-2) + dx1
      x(nx)   = x(nx-2) + 2.d0 * dx1
      x(nx+1) = x(nx-2) + 3.d0 * dx1
    
      ! y direction
      do j = 4, ny-2
        y(j) = dy1 * dble(j-4)
      enddo
      y(1)    = y(4)    - 3.d0 * dy1
      y(2)    = y(4)    - 2.d0 * dy1
      y(3)    = y(4)    - dy1
      y(ny-1) = y(ny-2) + dy1
      y(ny)   = y(ny-2) + 2.d0 * dy1
      y(ny+1) = y(ny-2) + 3.d0 * dy1

      ! z direction
      do k = 4, nz-2
        z(k) = dz1 * dble(k-4)
      enddo
      z(1)    = z(4)    - 3.d0 * dz1
      z(2)    = z(4)    - 2.d0 * dz1
      z(3)    = z(4)    - dz1
      z(nz-1) = z(nz-2) + dz1
      z(nz)   = z(nz-2) + 2.d0 * dz1
      z(nz+1) = z(nz-2) + 3.d0 * dz1
    endif

    ! cell centered
    do i = 1, nx
      xc(i) = 0.5d0 * (x(i) + x(i+1))
    enddo
    do j = 1, ny
      yc(j) = 0.5d0 * (y(j) + y(j+1))
    enddo
    do k = 1, nz
      zc(k) = 0.5d0 * (z(k) + z(k+1))
    enddo
  end subroutine set_grid
  
  subroutine set_init(nx,ny,nz,x,y,z,Q)
    use mod_globals, only : id_accuracy
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer i, j, k, offset
    if (kind(id_accuracy) == 2) then
      offset = 1
    elseif (kind(id_accuracy) == 4) then
      offset = 2
    elseif (kind(id_accuracy) == 8) then
      offset = 3
    endif
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          ! rho
          Q(i,j,k,1) =  RHO0
          ! rho u
          Q(i,j,k,2) =  RHO0 * M0 * sin(x(i)) * cos(y(j)) * cos(z(k))
          ! rho v
          Q(i,j,k,3) = -RHO0 * M0 * cos(x(i)) * sin(y(j)) * cos(z(k))
          ! rho w0
          Q(i,j,k,4) = 0.d0
          ! p / (gamma - 1) + 0.5 * (rhou ** 2 + rhov ** 2 ) / rho
          Q(i,j,k,5) = (1.d0/gamma+0.0625d0*RHO0*(M0**2)*(cos(2.d0*x(i))+cos(2.d0*y(j)))*(cos(2.d0*z(k))+2.d0))&
                       / (gamma - 1.d0) + 0.5d0 * (Q(i,j,k,2)**2 + Q(i,j,k,3)**2 + Q(i,j,k,4)**2) / Q(i,j,k,1)
    enddo;enddo;enddo
    if (kind(id_accuracy) == 2) then
      call set_bc_init2(Q)
    elseif (kind(id_accuracy) == 4) then
      call set_bc_init4(Q)
    elseif (kind(id_accuracy) == 8) then
      call set_bc_init6(Q)
    endif
  end subroutine set_init
  
  subroutine set_bc_init2(Q)
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer i, j, k

    do k = 2, nz-1
      do j = 2, ny-1
        Q(1,j,k,:)  = Q(nx-1,j,k,:)
        Q(nx,j,k,:) = Q(2,j,k,:)
    enddo;enddo

    do k = 2, nz-1
      do i = 2, nx-1
        Q(i,1,k,:)  = Q(i,ny-1,k,:)
        Q(i,ny,k,:) = Q(i,2,k,:)
    enddo;enddo

    do k = 2, nz-1
      Q(1,1,k,:)   = Q(nx-1,ny-1,k,:)
      Q(nx,1,k,:)  = Q(2,ny-1,k,:)
      Q(1,ny,k,:)  = Q(nx-1,2,k,:)
      Q(nx,ny,k,:) = Q(2,2,k,:)
    enddo

    do j = 1, ny
      do i = 1, nx
        Q(i,j,1,:)  = Q(i,j,nz-1,:)
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

  subroutine set_bc_init6(Q)
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer i, j, k

    do k = 4, nz-3
      do j = 4, ny-3
        Q(1:3,j,k,:) = Q(nx-5:nx-3,j,k,:)
        Q(nx-2:nx,j,k,:) = Q(4:6,j,k,:)
    enddo;enddo

    do k = 4, nz-3
      do i = 4, nx-3
        Q(i,1:3,k,:) = Q(i,ny-5:ny-3,k,:)
        Q(i,ny-2:ny,k,:) = Q(i,4:6,k,:)
    enddo;enddo

    do k = 4, nz-3
      Q(1:3,1:3,k,:) = Q(nx-5:nx-3,ny-5:ny-3,k,:)
      Q(nx-2:nx,1:3,k,:) = Q(4:6,ny-5:ny-3,k,:)
      Q(1:3,ny-2:ny,k,:) = Q(nx-5:nx-3,4:6,k,:)
      Q(nx-2:nx,ny-2:ny,k,:) = Q(4:6,4:6,k,:)
    enddo

    do j = 1, ny
      do i = 1, nx
        Q(i,j,1:3,:) = Q(i,j,nz-5:nz-3,:)
        Q(i,j,nz-2:nz,:) = Q(i,j,4:6,:)
    enddo;enddo
  end subroutine set_bc_init6

  subroutine set_bc(nx,ny,nz,Jacobian,Q,Qre)
    use mod_globals, only : id_accuracy
    integer, intent(in), value      :: nx, ny, nz
    real(8), intent(in), device     :: Jacobian(nx,ny,nz)
    real(8), intent(inout), device  :: Q(nx,ny,nz,5)
    real(8), intent(in), device     :: Qre(2,ny,nz,5)
    integer i, j, k, l
    if (kind(id_accuracy) == 2) then
      !$cuf kernel do(3)<<<*,*>>>
      do l = 1, 5
        do k = 2, nz-1
          do j = 2, ny-1
            Q(1,j,k,l)  = Q(nx-1,j,k,l)
            Q(nx,j,k,l) = Q(2,j,k,l)
      enddo;enddo;enddo

      !$cuf kernel do(3)<<<*,*>>>
      do l = 1, 5
        do k = 2, nz-1
          do i = 2, nx-1
            Q(i,1,k,l)  = Q(i,ny-1,k,l)
            Q(i,ny,k,l) = Q(i,2,k,l)
      enddo;enddo;enddo

      !$cuf kernel do(2)<<<*,*>>>
      do l = 1, 5
        do k = 2, nz-1
          Q(1,1,k,l)   = Q(nx-1,ny-1,k,l)
          Q(nx,1,k,l)  = Q(2,ny-1,k,l)
          Q(1,ny,k,l)  = Q(nx-1,2,k,l)
          Q(nx,ny,k,l) = Q(2,2,k,l)
      enddo;enddo

      !$cuf kernel do(3)<<<*,*>>>
      do l = 1, 5
        do j = 1, ny
          do i = 1, nx
            Q(i,j,1,l)  = Q(i,j,nz-1,l)
            Q(i,j,nz,l) = Q(i,j,2,l)
      enddo;enddo;enddo
    elseif (kind(id_accuracy) == 4) then
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
    elseif (kind(id_accuracy) == 8) then

      !$cuf kernel do(3)<<<*,*>>>
      do l = 1, 5
        do k = 4, nz-3
          do j = 4, ny-3
            Q(1,j,k,l) = Q(nx-5,j,k,l)
            Q(2,j,k,l) = Q(nx-4,j,k,l)
            Q(3,j,k,l) = Q(nx-3,j,k,l)
            Q(nx-2,j,k,l) = Q(4,j,k,l)
            Q(nx-1,j,k,l) = Q(5,j,k,l)
            Q(nx,j,k,l)   = Q(6,j,k,l)
      enddo;enddo;enddo

      !$cuf kernel do(3)<<<*,*>>>
      do l = 1, 5
        do k = 4, nz-3
          do i = 4, nx-3
            Q(i,1,k,l) = Q(i,ny-5,k,l)
            Q(i,2,k,l) = Q(i,ny-4,k,l)
            Q(i,3,k,l) = Q(i,ny-3,k,l)
            Q(i,ny-2,k,l) = Q(i,4,k,l)
            Q(i,ny-1,k,l) = Q(i,5,k,l)
            Q(i,ny,k,l)   = Q(i,6,k,l)
      enddo;enddo;enddo

      !$cuf kernel do(2)<<<*,*>>>
      do l = 1, 5
        do k = 4, nz-3
          Q(1,1,k,l) = Q(nx-5,ny-5,k,l)
          Q(1,2,k,l) = Q(nx-5,ny-4,k,l)
          Q(1,3,k,l) = Q(nx-5,ny-3,k,l)
          Q(2,1,k,l) = Q(nx-4,ny-5,k,l)
          Q(2,2,k,l) = Q(nx-4,ny-4,k,l)
          Q(2,3,k,l) = Q(nx-4,ny-3,k,l)
          Q(3,1,k,l) = Q(nx-3,ny-5,k,l)
          Q(3,2,k,l) = Q(nx-3,ny-4,k,l)
          Q(3,3,k,l) = Q(nx-3,ny-3,k,l)
          Q(nx-2,1,k,l) = Q(4,ny-5,k,l)
          Q(nx-2,2,k,l) = Q(4,ny-4,k,l)
          Q(nx-2,3,k,l) = Q(4,ny-3,k,l)
          Q(nx-1,1,k,l) = Q(5,ny-5,k,l)
          Q(nx-1,2,k,l) = Q(5,ny-4,k,l)
          Q(nx-1,3,k,l) = Q(5,ny-3,k,l)
          Q(nx,1,k,l)   = Q(6,ny-5,k,l)
          Q(nx,2,k,l)   = Q(6,ny-4,k,l)
          Q(nx,3,k,l)   = Q(6,ny-3,k,l)
          Q(1,ny-2,k,l) = Q(nx-5,4,k,l)
          Q(1,ny-1,k,l) = Q(nx-5,5,k,l)
          Q(1,ny,k,l)   = Q(nx-5,6,k,l)
          Q(2,ny-2,k,l) = Q(nx-4,4,k,l)
          Q(2,ny-1,k,l) = Q(nx-4,5,k,l)
          Q(2,ny,k,l)   = Q(nx-4,6,k,l)
          Q(3,ny-2,k,l) = Q(nx-3,4,k,l)
          Q(3,ny-1,k,l) = Q(nx-3,5,k,l)
          Q(3,ny,k,l)   = Q(nx-3,6,k,l)
          Q(nx-2,ny-2,k,l) = Q(4,4,k,l)
          Q(nx-2,ny-1,k,l) = Q(4,5,k,l)
          Q(nx-2,ny,k,l)   = Q(4,6,k,l)
          Q(nx-1,ny-2,k,l) = Q(5,4,k,l)
          Q(nx-1,ny-1,k,l) = Q(5,5,k,l)
          Q(nx-1,ny,k,l)   = Q(5,6,k,l)
          Q(nx,ny-2,k,l)   = Q(6,4,k,l)
          Q(nx,ny-1,k,l)   = Q(6,5,k,l)
          Q(nx,ny,k,l)     = Q(6,6,k,l)
      enddo;enddo

      !$cuf kernel do(3)<<<*,*>>>
      do l = 1, 5
        do j = 1, ny
          do i = 1, nx
            Q(i,j,1,l) = Q(i,j,nz-5,l)
            Q(i,j,2,l) = Q(i,j,nz-4,l)
            Q(i,j,3,l) = Q(i,j,nz-3,l)
            Q(i,j,nz-2,l) = Q(i,j,4,l)
            Q(i,j,nz-1,l) = Q(i,j,5,l)
            Q(i,j,nz,l)   = Q(i,j,6,l)
      enddo;enddo;enddo
    endif
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

