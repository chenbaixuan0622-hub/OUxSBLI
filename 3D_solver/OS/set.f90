module set
  use cudafor
  use mpi
  use mod_globals, only : Lx, Ly, Lz, gamma, R, u0, T0, &
  & beta, theta, M0, Ms, Ms2, a1, rho0, rho2, p0, p2, u1, u2, v1, v2, u_magnitude, ux, uy
  implicit none
contains
  subroutine set_grid(nx,ny,nz,x,y,z,dx,dy,dz)
    integer, intent(in)   :: nx, ny, nz
    real(8), intent(out)  :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    real(8) dx1, dy1, dz1
    integer i, j, k
    dx1 = Lx / dble(nx-1)
    dy1 = Ly / dble(ny-1)
    dz1 = Lz / dble(nz-1)
    x(1) = 0.d0
    do i = 1, nx-1
      dx(i)  = dx1
      x(i+1) = x(i) + dx(i)
    enddo

    y(1) = 0.d0
    do j = 1, ny-1
      dy(j)  = dy1
      y(j+1) = y(j) + dy(j)
    enddo

    z(1) = 0.d0
    do k = 1, nz-1
      dz(k)  = dz1
      z(k+1) = z(k) + dz(k)
    enddo
  end subroutine set_grid

  subroutine set_init(nx,ny,nz,xs,ys,zs,Q)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: xs(nx), ys(ny), zs(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer i, j, k, No
    real(8) rho, u, v, w, T
    rho = p0 / (R * T0)
    Q(:,:,:,1) = rho
    Q(:,:,:,2) = rho * u0
    Q(:,:,:,3) = 0.d0
    Q(:,:,:,4) = 0.d0
    Q(:,:,:,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho * u0**2
    ! oblique shock
    No = int(0.1 * nx)
    do k = 2, nz-1
      do i = No+1, nx
        Q(i,ny,k,1) = rho2
        Q(i,ny,k,2) = rho2 * ux
        Q(i,ny,k,3) = rho2 * uy
        Q(i,ny,k,4) = 0.d0
        Q(i,ny,k,5) = p2 / (gamma - 1.d0) + 0.5d0 * rho2 * (ux**2 + uy**2)
    enddo;enddo
  end subroutine set_init

  subroutine set_bc(nx,ny,nz,Jacobian,QJ,Qre)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny,nz)
    real(8), intent(inout), device :: QJ(nx,ny,nz,5)
    real(8), intent(in), device    :: Qre(2,ny,nz,5)
    ! Riemann boundary condition
    real(8) :: pin, cin, vin, Rp, Rm, vb, rhob, cb, pb, v0 = 0.d0, c0 = sqrt(gamma * p0 / rho0)
    integer i, j, k, l, No
    No = int(0.1 * nx)
    !$cuf kernel do(2)<<<*,*>>>
    do k = 3, nz-2
      do j = 2, ny-1
        ! inlet
        QJ(1,j,k,1)  = rho0 / Jacobian(1,j,k)
        QJ(1,j,k,2)  = rho0 * u0 / Jacobian(1,j,k)
        QJ(1,j,k,3)  = 0.d0
        QJ(1,j,k,4)  = 0.d0
        QJ(1,j,k,5)  = (p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2) / Jacobian(1,j,k)
        ! outlet
        QJ(nx,j,k,1) = QJ(nx-1,j,k,1)
        QJ(nx,j,k,2) = QJ(nx-1,j,k,2)
        QJ(nx,j,k,3) = QJ(nx-1,j,k,3)
        QJ(nx,j,k,4) = QJ(nx-1,j,k,4)
        QJ(nx,j,k,5) = QJ(nx-1,j,k,5)
    enddo;enddo

    ! Riemann boundary condition
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do i = 1, No
        pin = (gamma - 1.d0) * (QJ(i,ny-1,k,5) - 0.5d0 * (QJ(i,ny-1,k,2)**2 + QJ(i,ny-1,k,3)**2 + QJ(i,ny-1,k,4)**2) / QJ(i,ny-1,k,1)) &
        & * Jacobian(i,ny-1,k)
        cin = sqrt(gamma * pin / (QJ(i,ny-1,k,1) * Jacobian(i,ny-1,k)))
        vin = QJ(i,ny-1,k,3) / QJ(i,ny-1,k,1)
        Rp = vin + 2.d0 * cin / (gamma - 1.d0)
        Rm = v0  - 2.d0 * c0  / (gamma - 1.d0)
        vb = v0 + (0.5d0 * (Rp + Rm) - v0)

        rhob = QJ(i,ny-1,k,1) * Jacobian(i,ny-1,k)
        QJ(i,ny,k,1) = rhob / Jacobian(i,ny,k)
        QJ(i,ny,k,2) = QJ(i,ny-1,k,1) * u0 
        QJ(i,ny,k,3) = QJ(i,ny-1,k,1) * vb
        QJ(i,ny,k,4) = 0.d0
        cb = 0.25d0 * (gamma - 1.d0) * (Rp - Rm)
        pb = (rhob * cb**2) / gamma
        QJ(i,ny,k,5) = (pb / (gamma - 1.d0)) / Jacobian(i,ny,k)  + 0.5d0 * (QJ(i,ny,k,2)**2 + QJ(i,ny,k,3)**2 + QJ(i,ny,k,4)**2) / QJ(i,ny,k,1)
    enddo;enddo

    ! oblique shock
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do i = No+1, nx
        QJ(i,ny,k,1) = rho2 / Jacobian(i,ny,k)
        QJ(i,ny,k,2) = rho2 * ux / Jacobian(i,ny,k)
        QJ(i,ny,k,3) = rho2 * uy / Jacobian(i,ny,k)
        QJ(i,ny,k,4) = 0.d0
        QJ(i,ny,k,5) = (p2 / (gamma - 1.d0) + 0.5d0 * rho2 * (ux**2 + uy**2)) / Jacobian(i,ny,k)
    enddo;enddo

    ! bottom slip-wall
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        QJ(i,1,k,1) = QJ(i,2,k,1)
        QJ(i,1,k,2) = QJ(i,2,k,2)
        QJ(i,1,k,3) = 0.d0
        QJ(i,1,k,4) = 0.d0
        QJ(i,1,k,5) = QJ(i,2,k,5)
    enddo;enddo

    ! cyclic
    !$cuf kernel do(3)<<<*,*>>>
    do l = 1, 5
      do j = 1, ny
        do i = 1, nx
          QJ(i,j,1,l)    = QJ(i,j,nz-5,l)
          QJ(i,j,2,l)    = QJ(i,j,nz-4,l)
          QJ(i,j,3,l)    = QJ(i,j,nz-3,l)
          QJ(i,j,nz-2,l) = QJ(i,j,4,l)
          QJ(i,j,nz-1,l) = QJ(i,j,5,l)
          QJ(i,j,nz,l)   = QJ(i,j,6,l)
    enddo;enddo;enddo
  end subroutine set_bc

  subroutine set_bc_mut(nx,ny,nz,mut,qc2)
    integer, intent(in), value      :: nx, ny, nz
    real(8), intent(inout), device  :: mut(nx,ny,nz), qc2(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2) <<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        ! inlet
        mut(1,j,k) = mut(2,j,k)
        qc2(1,j,k) = qc2(2,j,k)
        ! outlet
        mut(nx,j,k) = mut(nx-1,j,k)
        qc2(nx,j,k) = qc2(nx-1,j,k)
    enddo;enddo

    !$cuf kernel do(2) <<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        ! wall
        mut(i,1,k) = 0.d0
        qc2(i,1,k) = 0.d0
        ! top
        mut(i,ny,k) = mut(i,ny-1,k)
        qc2(i,ny,k) = qc2(i,ny-1,k)
    enddo;enddo

    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        ! span
        mut(i,j,1)    = mut(i,j,nz-5)
        mut(i,j,2)    = mut(i,j,nz-4)
        mut(i,j,3)    = mut(i,j,nz-3)
        mut(i,j,nz-2) = mut(i,j,4)
        mut(i,j,nz-1) = mut(i,j,5)
        mut(i,j,nz)   = mut(i,j,6)
        qc2(i,j,1)    = qc2(i,j,nz-5)
        qc2(i,j,2)    = qc2(i,j,nz-4)
        qc2(i,j,3)    = qc2(i,j,nz-3)
        qc2(i,j,nz-2) = qc2(i,j,4)
        qc2(i,j,nz-1) = qc2(i,j,5)
        qc2(i,j,nz)   = qc2(i,j,6)
    enddo;enddo
  end subroutine set_bc_mut
end module set

