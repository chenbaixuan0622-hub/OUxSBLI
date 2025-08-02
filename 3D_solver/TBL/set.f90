module set
  use mod_globals, only : id_rescale, nx, ny, nz, nre2, gamma, R, Cp, Pr, u0, p0, T0, M0, blt, rho2, p2, ux, uy
  use set_init_common
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    integer i, j, k
    real(8) dx1, dy1, dz1
    dx1 = Lx / dble(nx-1)
    dy1 = dx1
    dz1 = Lz / dble(nz-1)
    x(1) = 0.d0
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo
    
    y(1) = 0.d0
    do j = 1, ny
      if (y(j) <= 2.d0 * blt) then
        dy(j) = min(1.d0, max(0.07d0, dble(j)/dble(128))) * dy1
      elseif (2.d0 * blt <= y(j) .and. y(j) <= 8.d0 * blt) then
        dy(j) = 1.5d0 * dy1
      else
        dy(j) = 1.75d0 * dy1
      endif
      y(j+1) = y(j) + dy(j)
    enddo

    z(1) = 0.d0
    do k = 1, nz-1
      dz(k) = dz1
      z(k+1) = z(k) + dz(k)
    enddo
  end subroutine set_grid

  subroutine set_init(myrank, nx, ny, nz, xs, ys, zs, Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: xs(nx), ys(ny), zs(nz)
    real(8), intent(out) :: Q(5,nx,ny,nz)
    real(8) :: rf = 0.89d0
    call set_init_tbl(nx, ny, nz, xs, ys, zs, 0.75d0*blt, blt, rf, u0, p0, T0, M0, Q)
  end subroutine set_init
  
  subroutine set_bc(myrank, nx, ny, nz, Jacobian, QJ, Qre)
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(ny)
    real(8), intent(inout), device :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(in), device    :: Qre(ny*(nz-6)*5)
    integer i, j, k, l, No
    real(8) :: Cp = gamma * R / (gamma - 1.d0), rf = 0.89d0
    real(8) :: p_wall
    ! Riemann invariants
    real(8) :: rhoin, pin, cin, vin, Rp, Rm, rhob, ub, vb, cb, pb
    real(8) :: rho0, c0, v0 = 0.d0, Taw, T
    No = int(0.4 * nx)

    if (kind(id_rescale) == 4) then
      !!$cuf kernel do(3)<<<*,*>>>
      !do l = 1, 5
      !  do k = 1, nz-6
      !    do j = 2, ny-1
      !      ! inlet
      !      QJ(1,j,k+3,l)  = Qre(ny*(nz-6)*(l-1)+ny*(k-1)+j)
      !      ! outlet
      !      QJ(nx,j,k+3,l) = QJ(nx-1,j,k+3,l)
      !enddo;enddo;enddo
    else
      !$cuf kernel do(3)<<<*,*>>>
      do k = 4, nz-3
        do j = 2, ny-1
          do l = 1, 5
            ! inlet
            QJ(l,1,j,k) = QJ(l,nx-5,j,k)
            QJ(l,2,j,k) = QJ(l,nx-4,j,k)
            QJ(l,3,j,k) = QJ(l,nx-3,j,k)
            ! outlet
            QJ(l,nx-2,j,k) = QJ(l,4,j,k)
            QJ(l,nx-1,j,k) = QJ(l,5,j,k)
            QJ(l,nx,j,k)   = QJ(l,6,j,k)
      enddo;enddo;enddo
    endif

    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz-3
      do i = 1, nx
        ! top
        ! Riemann invariants
        pin   = (gamma - 1.d0) * (QJ(5,i,ny-1,k) - 0.5d0 * (QJ(2,i,ny-1,k)**2 + QJ(3,i,ny-1,k)**2 + QJ(4,i,ny-1,k)**2) &
                / QJ(1,i,ny-1,k)) * Jacobian(ny-1)
        rhoin = QJ(1,i,ny-1,k) * Jacobian(ny-1)
        cin   = sqrt(gamma * pin / rhoin)
        vin   = QJ(3,i,ny-1,k) / QJ(1,i,ny-1,k)
        ! temperature and density at top
        Taw   = T0 * (1.d0 + rf * 0.5d0 * (gamma - 1.d0) * M0**2)
        T     = Taw - rf * u0**2 / (2.d0 * (gamma * R / (gamma - 1.d0)))
        rho0  = p0 / (R * T)
        c0    = sqrt(gamma * p0 / rho0)

        Rp    = vin + 2.d0 * cin / (gamma - 1.d0)
        Rm    = v0  - 2.d0 * c0  / (gamma - 1.d0)
        vb    = 0.5d0 * (Rp + Rm)
        cb    = 0.25d0 * (gamma - 1.d0) * (Rp - Rm)
        rhob  = (cb / c0)**(2.d0 / (gamma - 1.d0)) * rho0
        pb    = (rhob * cb**2) / gamma

        QJ(1,i,ny,k) = rhob / Jacobian(ny)
        QJ(2,i,ny,k) = rhob * u0 / Jacobian(ny)
        QJ(3,i,ny,k) = rhob * vb / Jacobian(ny)
        QJ(4,i,ny,k) = 0.d0
        QJ(5,i,ny,k) = (pb / (gamma - 1.d0) + 0.5d0 * rhob * (u0**2 + vb**2)) / Jacobian(ny)

        ! NoSlip
        QJ(1,i,1,k) = QJ(1,i,2,k)
        QJ(2,i,1,k) = 0.d0
        QJ(3,i,1,k) = 0.d0
        QJ(4,i,1,k) = 0.d0
        p_wall = (gamma - 1.d0) * (QJ(5,i,2,k) - 0.5d0 * (QJ(2,i,2,k)**2 + QJ(3,i,2,k)**2 + QJ(4,i,2,k)**2) / QJ(1,i,2,k))
        QJ(5,i,1,k) = p_wall / (gamma - 1.d0)
    enddo;enddo

!    !$cuf kernel do(2)<<<*,*>>>
!    do k = 1, nz
!      do i = No, nx/2
!        QJ(1,i,ny,k) = rho2 / Jacobian(ny)
!        QJ(2,i,ny,k) = rho2 * ux / Jacobian(ny)
!        QJ(3,i,ny,k) = rho2 * uy / Jacobian(ny)
!        QJ(4,i,ny,k) = 0.d0
!        QJ(5,i,ny,k) = (p2 / (gamma - 1.d0) + 0.5d0 * rho2 * (ux**2 + uy**2)) / Jacobian(ny)
!    enddo;enddo
!
!    !$cuf kernel do(2)<<<*,*>>>
!    do k = 1, nz
!      do i = nx/2, nx
!        pin   = (gamma - 1.d0) * (QJ(5,i,ny-1,k) - 0.5d0 * (QJ(2,i,ny-1,k)**2 + QJ(3,i,ny-1,k)**2 + QJ(4,i,ny-1,k)**2) &
!                / QJ(1,i,ny-1,k)) * Jacobian(ny-1)
!        rhoin = QJ(1,i,ny-1,k) * Jacobian(ny-1)
!        cin   = sqrt(gamma * pin / rhoin)
!        vin   = QJ(3,i,ny-1,k) / QJ(1,i,ny-1,k)
!        c0    = sqrt(gamma * p2 / rho2)
!        Rp    = vin + 2.d0 * cin / (gamma - 1.d0)
!        Rm    = uy  - 2.d0 * c0  / (gamma - 1.d0)
!        vb    = 0.5d0 * (Rp + Rm)
!        cb    = 0.25d0 * (gamma - 1.d0) * (Rp - Rm)
!        rhob  = cin * rhoin / cb
!        pb    = (rhob * cb**2) / gamma
!        ub    = sqrt(2.d0 * gamma * (p2 / rho2 - pb / rhob) / (gamma - 1.d0) + ux**2 + uy**2 - vb**2)
!        QJ(1,i,ny,k) = rhob / Jacobian(ny)
!        QJ(2,i,ny,k) = rhob * ub / Jacobian(ny)
!        QJ(3,i,ny,k) = rhob * vb / Jacobian(ny)
!        QJ(4,i,ny,k) = 0.d0
!        QJ(5,i,ny,k) = (pb / (gamma - 1.d0) + 0.5d0 * rhob * (ub**2 + vb**2)) / Jacobian(ny)
!    enddo;enddo
!
    ! cyclic
    !$cuf kernel do(3)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        do l = 1, 5
          QJ(l,i,j,1) = QJ(l,i,j,nz-5)
          QJ(l,i,j,2) = QJ(l,i,j,nz-4)
          QJ(l,i,j,3) = QJ(l,i,j,nz-3)
          QJ(l,i,j,nz-2) = QJ(l,i,j,4)
          QJ(l,i,j,nz-1) = QJ(l,i,j,5)
          QJ(l,i,j,nz)   = QJ(l,i,j,6)
    enddo;enddo;enddo
  end subroutine set_bc

  subroutine set_bc_mut(nx, ny, nz, mut, qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2) <<<*,*>>>
    do k = 4, nz-3
      do j = 2, ny-1
        ! inlet
        mut(1,j,k)  = mut(nre2,j,k)
        qc2(1,j,k)  = qc2(nre2,j,k)
        ! outlet
        mut(nx,j,k) = mut(nx-1,j,k)
        qc2(nx,j,k) = qc2(nx-1,j,k)
    enddo;enddo

    !$cuf kernel do(2) <<<*,*>>>
    do k = 4, nz-3
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

  subroutine calc_forcing(nx, ny, nz, dx, dy, dz, rho, u, v, w, p, fx, fy, fz)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: rho(nx,ny,nz), u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), p(nx,ny,nz)
    real(8), intent(out), device :: fx(nx-2,ny-2,nz-2), fy(nx-2,ny-2,nz-2), fz(nx-2,ny-2,nz-2)
  end subroutine calc_forcing
end module set

