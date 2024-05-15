module set
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, gamma, u0, &
  & beta, theta, Ms, Ms2, a1, rho0, rho2, p0, p2, u1, u2, v1, v2, u_magnitude, ux, uy
  implicit none
contains
  subroutine set_grid(x,y,z,dx,dy)
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx), dy(ny)
    integer i, j, k
    real(8) :: dx1 = Lx / dble(nx-1)
    real(8) :: dy1 = Ly / dble(ny-1)
    real(8) :: dz1 = Lz / dble(nz-1)
    do i = 1, nx
      x(i) = dble(i-1) * dx1
      dx(i) = dx1
    enddo
    do j = 1, ny
      y(j) = dble(j-1) * dy1
      dy(j) = dy1
    enddo
    do k = 1, nz
      z(k) = dble(k-1) * dz1
    enddo
  end subroutine set_grid

  subroutine set_init(xs,ys,zs,Q,Vin)
    real(8), intent(in) :: xs(nx), ys(ny), zs(nz)
    real(8), intent(out), dimension(nx,ny,nz,5) :: Q
    real(8), intent(in), dimension(ny,2) :: Vin
    integer i, k
    integer :: No = int(0.25 * nx)
    Q(:,:,:,1) = rho0
    Q(:,:,:,2) = rho0 * u0
    Q(:,:,:,3) = 0.d0
    Q(:,:,:,4) = 0.d0
    Q(:,:,:,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2

    ! top
    do k = 2, nz-1
      do i = 1, No
        Q(i,ny,k,1) = rho0
        Q(i,ny,k,2) = rho0 * u0 
        Q(i,ny,k,3) = 0.d0
        Q(i,ny,k,4) = 0.d0
        Q(i,ny,k,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
    enddo;enddo

    do k = 2, nz-1
      do i = No+1, nx
        Q(i,ny,k,1) = rho2
        Q(i,ny,k,2) = rho2 * ux
        Q(i,ny,k,3) = rho2 * uy
        Q(i,ny,k,4) = 0.d0
        Q(i,ny,k,5) = p2 / (gamma - 1.d0) + 0.5d0 * rho2 * (ux**2 + uy**2)
    enddo;enddo

    ! bottom
    Q(:,2,:,1) = Q(:,3,:,1)
    Q(:,2,:,2) = Q(:,3,:,2)
    Q(:,2,:,3) = 0.d0
    Q(:,2,:,4) = Q(:,3,:,4)
    Q(:,2,:,5) = Q(:,3,:,5)
    ! imaginary
    Q(:,1,:,1) = Q(:,3,:,1)
    Q(:,1,:,2) = Q(:,3,:,2)
    Q(:,1,:,3) = -Q(:,3,:,3)
    Q(:,1,:,4) = Q(:,3,:,4)
    Q(:,1,:,5) = Q(:,3,:,5)
  end subroutine set_init
  
  subroutine set_bc(id_accuracy,Q)
    integer(kind=2), intent(in), value :: id_accuracy
    real(8), intent(inout), device :: Q(nx,ny,nz,5)
    integer i, j, k, l
    integer :: No = int(0.25 * nx)
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do j = 3, ny-1
        ! inlet
        Q(1,j,k,1) = rho0
        Q(1,j,k,2) = rho0 * u0
        Q(1,j,k,3) = 0.d0
        Q(1,j,k,4) = 0.d0
        Q(1,j,k,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
        ! outlet
        Q(nx,j,k,1) = Q(nx-1,j,k,1)
        Q(nx,j,k,2) = Q(nx-1,j,k,2)
        Q(nx,j,k,3) = Q(nx-1,j,k,3)
        Q(nx,j,k,4) = Q(nx-1,j,k,4)
        Q(nx,j,k,5) = Q(nx-1,j,k,5)
    enddo;enddo

    ! top
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do i = 1, No
        Q(i,ny,k,1) = rho0
        Q(i,ny,k,2) = rho0 * u0 
        Q(i,ny,k,3) = 0.d0
        Q(i,ny,k,4) = 0.d0
        Q(i,ny,k,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
    enddo;enddo

    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do i = No+1, nx
        Q(i,ny,k,1) = rho2
        Q(i,ny,k,2) = rho2 * ux
        Q(i,ny,k,3) = rho2 * uy
        Q(i,ny,k,4) = 0.d0
        Q(i,ny,k,5) = p2 / (gamma - 1.d0) + 0.5d0 * rho2 * (ux**2 + uy**2)
    enddo;enddo

    ! bottom
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do i = 1, nx
        ! slip wall
        Q(i,2,k,1) = Q(i,3,k,1)
        Q(i,2,k,2) = Q(i,3,k,2)
        Q(i,2,k,3) = 0.d0
        Q(i,2,k,4) = Q(i,3,k,4)
        Q(i,2,k,5) = Q(i,3,k,5)
        ! imaginary
        Q(i,1,k,1) = Q(i,3,k,1)
        Q(i,1,k,2) = Q(i,3,k,2)
        Q(i,1,k,3) = -Q(i,3,k,3)
        Q(i,1,k,4) = Q(i,3,k,4)
        Q(i,1,k,5) = Q(i,3,k,5)
    enddo;enddo

    ! cyclic
    !$cuf kernel do(3)<<<*,*>>>
    do l = 1, 5
      do j = 1, ny
        do i = 1, nx
          Q(i,j,1,l) = Q(i,j,nz-3,l)
          Q(i,j,2,l) = Q(i,j,nz-2,l)
          Q(i,j,nz-1,l) = Q(i,j,3,l)
          Q(i,j,nz,l) = Q(i,j,4,l)
    enddo;enddo;enddo
  end subroutine set_bc

  subroutine set_bc_mut(mut)
    real(8), intent(inout), device :: mut(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2) <<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        ! inlet
        mut(1,j,k) = mut(2,j,k)
        ! outlet
        mut(nx,j,k) = mut(nx-1,j,k)
    enddo;enddo

    !$cuf kernel do(2) <<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        ! wall
        mut(i,1,k) = 0.d0
        ! top
        mut(i,ny,k) = mut(i,ny-1,k)
    enddo;enddo

    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        ! span
        mut(i,j,1) = mut(i,j,2)
        mut(i,j,nz) = mut(i,j,nz-1)
    enddo;enddo
  end subroutine set_bc_mut
end module set

