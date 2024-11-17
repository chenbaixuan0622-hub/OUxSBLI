module set
  use mod_globals, only : id_rescale, nx, ny, nz, Lx, Ly, Lz, gamma, R, Cp, Pr, u0, p0, T0, M0, blt, rho2, p2, ux, uy
  implicit none
contains
  subroutine set_grid(nx,ny,nz,x,y,z,dx,dy,dz)
    integer, intent(in)   :: nx, ny, nz
    real(8), intent(out)  :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    integer i, j, k
    real(8) dx1, dy1, dz1
    dx1 = Lx / dble(nx-1)
    dy1 = 10.d-3 / dble(256)
    dz1 = Lz / dble(nz-1)
    x(1) = 0.d0
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo

    y(1) = 0.d0
    do j = 1, 256
      ! DNS
      dy(j)  = min(1.d0, max(0.05d0, dble(j)/dble(128))) * dy1
      y(j+1) = y(j) + dy(j)
    enddo
    ! buffer
    do j = 257, ny-1
      dy(j)  = 2.d0 * dy1
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
    real(8) :: blt0 = 0.5d0 * blt
    real(8) :: Cp   = gamma * R / (gamma - 1.d0)
    real(8) :: eta, rho, u, v, w, T, Tw, Taw, p_wall

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          eta = 5.d0 * ys(j) / blt0
          u   = min(u0, u0 * (0.0015d0 * eta**4 - 0.0181d0 * eta**3 + 0.029d0 * eta**2 + 0.3192 * eta + 0.0003d0))
          v   = 0.d0
          Taw = T0 + 0.5d0 * u0**2 / Cp
          Tw  = Taw
          T   = Tw + (Taw - Tw) * u / u0 - 0.5d0 * Pr**(1.d0/3.d0) * u**2 / Cp
          rho = p0 / (R * T)
          Q(i,j,k,1) = rho
          Q(i,j,k,2) = Q(i,j,k,1) * u
          Q(i,j,k,3) = Q(i,j,k,1) * v
          Q(i,j,k,4) = 0.d0
          Q(i,j,k,5) = p0 / (gamma - 1.d0) + 0.5d0 * (Q(i,j,k,2)**2 + Q(i,j,k,3)**2 + Q(i,j,k,4)**2) / Q(i,j,k,1)
    enddo;enddo;enddo

    ! bottom
    Q(:,1,:,1) = Q(:,2,:,1)
    Q(:,1,:,2) = 0.d0
    Q(:,1,:,3) = 0.d0
    Q(:,1,:,4) = 0.d0
    p_wall = (gamma - 1.d0) * (Q(2,2,2,5) - 0.5d0 * (Q(2,2,2,2)**2 + Q(2,2,2,3)**2 + Q(2,2,2,4)**2) / Q(2,2,2,1))
    Q(:,1,:,5) = p_wall / (gamma - 1.d0)
  end subroutine set_init
  
  subroutine set_bc(nx,ny,Jacobian,QJ)
    integer, intent(in), value     :: nx, ny
    real(8), intent(in), device    :: Jacobian(ny)
    real(8), intent(inout), device :: QJ(nx,ny,4) ! Q / Jacobian
    integer i, j, k, No
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    real(8) :: p_wall
    ! Riemann invariants
    real(8) :: pin, cin, vin, Rp, Rm, rhob, vb, cb, pb
    real(8) :: Taw, Tw, T, v0 = 0.d0
    No = int(0.4 * nx)

    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, 4
      do j = 2, ny-1
        ! inlet
        QJ(1,j,k) = QJ(No,j,k)
        ! outlet
        QJ(nx,j,k) = QJ(nx-1,j,k)
    enddo;enddo

    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      ! Neumann boundary condition
      QJ(i,ny,1) = QJ(i,ny-1,1)
      QJ(i,ny,2) = QJ(i,ny-1,2)
      QJ(i,ny,3) = QJ(i,ny-1,3)
      QJ(i,ny,4) = QJ(i,ny-1,4)
      ! NoSlip
      QJ(i,1,1) = QJ(i,2,1)
      QJ(i,1,2) = 0.d0
      QJ(i,1,3) = 0.d0
      p_wall = (gamma - 1.d0) * (QJ(i,2,4) - 0.5d0 * (QJ(i,2,2)**2 + QJ(i,2,3)**2) / QJ(i,2,1))
      QJ(i,1,4) = p_wall / (gamma - 1.d0)
    enddo

    !!$cuf kernel do(1)<<<*,*>>>
    !do i = No, nx
    !  QJ(i,ny,1) = rho2 / Jacobian(ny)
    !  QJ(i,ny,2) = rho2 * ux / Jacobian(ny)
    !  QJ(i,ny,3) = rho2 * uy / Jacobian(ny)
    !  QJ(i,ny,4) = (p2 / (gamma - 1.d0) + 0.5d0 * rho2 * (ux**2 + uy**2)) / Jacobian(ny)
    !enddo
  end subroutine set_bc
end module set

