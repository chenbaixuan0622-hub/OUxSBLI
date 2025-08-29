module set
  use cudafor
  use mpi
  use mod_globals, only : id_rescale, ny1, nre2, gamma, R, Cp, Pr, u0, p0, T0, M0, blt, rho2, p2, ux, uy
  use mod_constant, only : Cp, gamma_1, over_gamma_1
  use set_bc_common
  use set_init_common
  use calc_para
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
    do j = 1, ny-1
      if (y(j) <= 3.d0 * blt) then
        dy(j) = min(1.d0, max(0.07d0, dble(j)/dble(128))) * dy1
      elseif (3.d0 * blt <= y(j) .and. y(j) <= 8.d0 * blt) then
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
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(in), device, optional :: Qre(ny*(nz-6)*5)
    integer i, j, k, l, No, ireq, ierr, istat(MPI_STATUS_SIZE)
    real(8) :: p_wall, rf = 0.89d0
    ! Riemann invariants
    real(8) :: rhoin, pin, cin, vin, Rp, Rm, rhob, ub, vb, cb, pb
    real(8) :: rho0, c0, v0 = 0.d0, Taw, T
    ! cache
    real(8) Jacobian_tmp
    ! temperature and density at top
    Taw  = T0 * (1.d0 + rf * 0.5d0 * gamma_1 * M0**2)
    T    = Taw - rf * u0**2 / (2.d0 * (gamma * R * over_gamma_1))
    rho0 = p0 / (R * T)
    c0   = sqrt(gamma * p0 / rho0)
    if (kind(id_rescale) == 4 .and. present(Qre)) then
      !block
      !  real(8) Qre_cpu(ny*(nz-6)*5)
      !  Qre_cpu = Qre
      !  do k = 1, nz-6
      !    do j = 2, ny-1
      !      do l = 1, 5
      !        print *, Qre_cpu(ny*5*(k-1)+5*(j-1)+l)
      !  enddo;enddo;enddo
      !end block
      !$cuf kernel do(2)<<<*,*>>>
      do k = 1, nz-6
        do j = 2, ny-1
          do l = 1, 5
            ! inlet
            QJ(l,1,j,k+3)  = Qre(ny*5*(k-1)+5*(j-1)+l)
            ! outlet
            QJ(l,nx,j,k+3) = QJ(l,nx-1,j,k+3)
      enddo;enddo;enddo
    else
      !$cuf kernel do(2)<<<*,*>>>
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
        Jacobian_tmp = 1.d0 / Jacobian(i,ny)
        pin   = gamma_1 * (QJ(5,i,ny-1,k) - 0.5d0 * (QJ(2,i,ny-1,k)**2 + QJ(3,i,ny-1,k)**2 + QJ(4,i,ny-1,k)**2) &
                / QJ(1,i,ny-1,k)) * Jacobian(i,ny-1)
        rhoin = QJ(1,i,ny-1,k) * Jacobian(i,ny-1)
        cin   = sqrt(gamma * pin / rhoin)
        vin   = QJ(3,i,ny-1,k) / QJ(1,i,ny-1,k)
        Rp   = vin + 2.d0 * cin * over_gamma_1
        Rm   = v0  - 2.d0 * c0  * over_gamma_1
        vb   = 0.5d0 * (Rp + Rm)
        cb   = 0.25d0 * gamma_1 * (Rp - Rm)
        rhob = (cb / c0)**(2.d0 * over_gamma_1) * rho0
        pb   = (rhob * cb**2) / gamma
        QJ(1,i,ny,k) = rhob * Jacobian_tmp
        QJ(2,i,ny,k) = rhob * u0 * Jacobian_tmp
        QJ(3,i,ny,k) = rhob * vb * Jacobian_tmp
        QJ(4,i,ny,k) = 0.d0
        QJ(5,i,ny,k) = (pb * over_gamma_1 + 0.5d0 * rhob * (u0**2 + vb**2)) * Jacobian_tmp
        ! NoSlip
        QJ(1,i,1,k) = QJ(1,i,2,k)
        QJ(2,i,1,k) = 0.d0
        QJ(3,i,1,k) = 0.d0
        QJ(4,i,1,k) = 0.d0
        p_wall = gamma_1 * (QJ(5,i,2,k) - 0.5d0 * (QJ(2,i,2,k)**2 + QJ(3,i,2,k)**2 + QJ(4,i,2,k)**2) / QJ(1,i,2,k))
        QJ(5,i,1,k) = p_wall * over_gamma_1
    enddo;enddo

    ! cyclic
    !$cuf kernel do(2)<<<*,*>>>
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

  subroutine set_bc_mut(nx,ny,nz,mut,qc2)
    integer, intent(in), value      :: nx, ny, nz
    real(8), intent(inout), device  :: mut(nx,ny,nz), qc2(nx,ny,nz)
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
end module set

