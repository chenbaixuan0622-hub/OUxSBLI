module set
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, gamma, R, dtn
  implicit none
contains
  subroutine set_grid(nx,ny,nz,xc,yc,zc,dx,dy,dz)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(out) :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    real(8) x(nx+1), y(ny+1), z(nz+1)
    integer i, j, k
    dx(:) = Lx / dble(nx-6)
    dy(:) = Ly / dble(ny-6)
    dz(:) = Lz / dble(nz-6)

    ! x direction
    do i = 4, nx-2
      x(i) = dx(1) * dble(i-4)
    enddo
    x(1)    = x(4)    - 3.d0 * dx(1)
    x(2)    = x(4)    - 2.d0 * dx(1)
    x(3)    = x(4)    - dx(1)
    x(nx-1) = x(nx-2) + dx(1)
    x(nx)   = x(nx-2) + 2.d0 * dx(1)
    x(nx+1) = x(nx-2) + 3.d0 * dx(1)

    ! y direction
    do j = 4, ny-2
      y(j) = dy(1) * dble(j-4)
    enddo
    y(1)    = y(4)    - 3.d0 * dy(1)
    y(2)    = y(4)    - 2.d0 * dy(1)
    y(3)    = y(4)    - dy(1)
    y(ny-1) = y(ny-2) + dy(1)
    y(ny)   = y(ny-2) + 2.d0 * dy(1)
    y(ny+1) = y(ny-2) + 3.d0 * dy(1)

    ! z direction
    do k = 4, nz-2
      z(k) = dz(1) * dble(k-4)
    enddo
    z(1)    = z(4)    - 3.d0 * dz(1)
    z(2)    = z(4)    - 2.d0 * dz(1)
    z(3)    = z(4)    - dz(1)
    z(nz-1) = z(nz-2) + dz(1)
    z(nz)   = z(nz-2) + 2.d0 * dz(1)
    z(nz+1) = z(nz-2) + 3.d0 * dz(1)

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
    use mod_globals, only : M0, rho0, p0, T0, u0, Rc, beta
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer i, j
    real(8) xc, yc, ex, T, rho, u, v, p
    real(8) :: Cp = R * gamma / (gamma - 1.d0)
    xc = x(int(nx/2))
    yc = y(int(ny/2))
    do j = 1, ny
      do i = 1, nx
        ex  = exp(-0.5d0 * ((x(i) - xc)**2 + (y(j) - yc)**2) / (Rc**2))
        T   = T0 - 0.5d0 * (u0 * beta)**2 / Cp * ex**2
        rho = rho0 * (T / T0)**(1.d0  / (gamma - 1.d0))
        u   = u0 * (1.d0 - beta * (y(j) - yc) / Rc * ex)
        v   = u0 *         beta * (x(i) - xc) / Rc * ex
        p   = rho * R * T
        Q(i,j,:,1) = rho
        Q(i,j,:,2) = rho * u
        Q(i,j,:,3) = rho * v
        Q(i,j,:,4) = 0.d0
        ! p / (gamma - 1) + 0.5 * (rhou ** 2 + rhov ** 2 ) / rho
        Q(i,j,:,5) = p / (gamma - 1.d0) + 0.5d0 * rho * (u**2 + v**2)
    enddo;enddo
  end subroutine set_init
  
  subroutine set_bc(nx,ny,nz,Jacobian,Q,Qre)
    integer, intent(in), value      :: nx, ny, nz
    real(8), intent(in), device     :: Jacobian(nx,ny,nz)
    real(8), intent(inout), device  :: Q(nx,ny,nz,5)
    real(8), intent(in), device     :: Qre(2,ny,nz,5)
    integer i, j, k, l
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
          Q(i,j,1,l) = Q(i,j,nz-3,l)
          Q(i,j,2,l) = Q(i,j,nz-3,l)
          Q(i,j,3,l) = Q(i,j,nz-3,l)
          Q(i,j,nz-2,l) = Q(i,j,4,l)
          Q(i,j,nz-1,l) = Q(i,j,4,l)
          Q(i,j,nz,l)   = Q(i,j,4,l)
    enddo;enddo;enddo
  end subroutine set_bc

  subroutine set_bc_mut(nx,ny,nz,mut,qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2) <<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        mut(1,j,k) = mut(nx-1,j,k)
        mut(nx,j,k) = mut(2,j,k)
        qc2(1,j,k) = qc2(nx-1,j,k)
        qc2(nx,j,k) = qc2(2,j,k)
    enddo;enddo

    !$cuf kernel do(2) <<<*,*>>>
    do k = 2, nz-1
      do i = 2, nx-1
        mut(i,1,k) = mut(i,ny-1,k)
        mut(i,ny,k) = mut(i,2,k)
        qc2(i,1,k) = qc2(i,ny-1,k)
        qc2(i,ny,k) = qc2(i,2,k)
    enddo;enddo

    !$cuf kernel do(1) <<<*,*>>>
    do k = 2, nz-1
      mut(1,1,k) = mut(nx-1,ny-1,k)
      mut(nx,1,k) = mut(2,ny-1,k)
      mut(1,ny,k) = mut(nx-1,2,k)
      mut(nx,ny,k) = mut(2,2,k)
      qc2(1,1,k) = qc2(nx-1,ny-1,k)
      qc2(nx,1,k) = qc2(2,ny-1,k)
      qc2(1,ny,k) = qc2(nx-1,2,k)
      qc2(nx,ny,k) = qc2(2,2,k)
    enddo

    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        mut(i,j,1) = mut(i,j,nz-1)
        mut(i,j,nz) = mut(i,j,2)
        qc2(i,j,1) = qc2(i,j,nz-1)
        qc2(i,j,nz) = qc2(i,j,2)
    enddo;enddo
  end subroutine set_bc_mut
end module set

