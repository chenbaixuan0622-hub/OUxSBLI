module set
  use cudafor
  use mpi
  use mod_globals, only : id_rescale, ny1, nre2, gamma, R, Cp, Pr, u0, p0, T0, M0, blt, rho2, p2, ux, uy, rf, Taw
  use mod_constant, only : Cp, gamma_1, over_gamma, over_gamma_1
  use set_bc_common
  use set_init_common
  use calc_para
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, Lx, Ly, x, y, dx, dy)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: Lx, Ly
    real(8), intent(out) :: x(nx), y(ny), dx(nx-1), dy(ny-1)
    integer i, j, ny_b
    real(8) dx1, dy1 !dz1
    dx1 = Lx / dble(nx-1)
    dy1 = dx1
    !dz1 = Lz / dble(nz-1)

    x(1) = 0.d0
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo

    y(1) = 0.d0
    do j = 1, ny-1
      if (y(j) <= 3.d0 * blt) then
        dy(j) = min(1.d0, max(0.07d0, dble(j)/dble(128))) * dy1
        ny_b  = j
      else
        dy(j) = dy1 * (1.d0 + 0.75d0 * dble(j-ny_b) / dble(ny-ny_b))
      endif
      y(j+1) = y(j) + dy(j)
    enddo

    ! z(1) = 0.d0
    ! do k = 1, nz-1
    !   dz(k) = dz1
    !   z(k+1) = z(k) + dz(k)
    ! enddo
    ! z(:) = z(:) - 0.5d0 * Lz
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, xs, ys, Q)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: xs(nx), ys(ny)
    real(8), intent(out) :: Q(nx,4,ny)
    call set_init_tbl(nx, ny, xs, ys, 0.1d0, 0.75d0*blt, blt, u0, p0, T0, M0, Q)
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, Jacobian, QJ, Qre)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ(nx,4,ny) ! Q / Jacobian
    real(8), intent(in), device, optional :: Qre(ny*5)   !Qre(ny*(nz-6)*5)
    integer i, j, l, ireq, ierr, istat(MPI_STATUS_SIZE)
    real(8) :: p_wall
    ! Riemann invariants
    real(8) :: rhoin, pin, cin, vin, Rp, Rm, rhob, ub, vb, cb, pb, v0 = 0.d0
    ! cache
    real(8) Jacobian_tmp
    ! temperature and density at top
    real(8), parameter :: T       = Taw - rf * u0**2 / (2.d0 * Cp)
    real(8), parameter :: rho0    = p0 / (R * T)
    real(8), parameter :: c0      = sqrt(gamma * p0 / rho0)
    real(8), parameter :: over_c0 = 1.d0 / c0
    if (kind(id_rescale) == 4 .and. present(Qre)) then
      !$cuf kernel do(1)<<<*,*>>>
      !do k = 1, nz-6
        do j = 2, ny-1
          do l = 1, 4
            ! inlet
            QJ(1,l,j)  = Qre(ny*5+5*(j-1)+l)
            ! outlet
            QJ(nx,l,j) = QJ(nx-1,l,j)
        enddo;enddo
    else
      !$cuf kernel do(1)<<<*,*>>>
      !do k = 4, nz-3
        do j = 2, ny-1
          do l = 1, 4
            ! inlet
            QJ(1,l,j) = QJ(nx-5,l,j)
            QJ(2,l,j) = QJ(nx-4,l,j)
            QJ(3,l,j) = QJ(nx-3,l,j)
            ! outlet
            QJ(nx-2,l,j) = QJ(4,l,j)
            QJ(nx-1,l,j) = QJ(5,l,j)
            QJ(nx,l,j)   = QJ(6,l,j)
        enddo;enddo
    endif

    !$cuf kernel do(1)<<<*,*>>>
    !do k = 4, nz-3
      do i = 1, nx
        ! top
        ! Riemann invariants
        Jacobian_tmp = 1.d0 / Jacobian(i,ny)
        pin   = gamma_1 * (QJ(i,4,ny-1) - 0.5d0 * (QJ(i,2,ny-1)**2 + QJ(i,3,ny-1)**2) &
                / QJ(i,1,ny-1)) * Jacobian(i,ny-1)
        rhoin = QJ(i,1,ny-1) * Jacobian(i,ny-1)
        cin   = sqrt(gamma * pin / rhoin)
        vin   = QJ(i,3,ny-1) / QJ(i,1,ny-1)
        Rp   = vin + 2.d0 * cin * over_gamma_1
        Rm   = v0  - 2.d0 * c0  * over_gamma_1
        vb   = 0.5d0 * (Rp + Rm)
        cb   = 0.25d0 * gamma_1 * (Rp - Rm)
        rhob = (cb * over_c0)**(2.d0 * over_gamma_1) * rho0
        pb   = (rhob * cb**2) * over_gamma
        QJ(i,1,ny) = rhob * Jacobian_tmp
        QJ(i,2,ny) = rhob * u0 * Jacobian_tmp
        QJ(i,3,ny) = rhob * vb * Jacobian_tmp
        !QJ(i,4,ny) = 0.d0
        QJ(i,4,ny) = (pb * over_gamma_1 + 0.5d0 * rhob * (u0**2 + vb**2)) * Jacobian_tmp
        ! NoSlip
        QJ(i,1,1) = QJ(i,1,2)
        QJ(i,2,1) = 0.d0
        QJ(i,3,1) = 0.d0
        !QJ(i,4,1) = 0.d0
        p_wall = gamma_1 * (QJ(i,4,2) - 0.5d0 * (QJ(i,2,2)**2 + QJ(i,3,2)**2) / QJ(i,1,2))
        QJ(i,4,1) = p_wall * over_gamma_1
      enddo

    ! cyclic
    !cuf kernel do(2)<<<*,*>>>
    ! do j = 1, ny
    !   do i = 1, nx
    !     do l = 1, 5
    !       QJ(i,l,j,1) = QJ(i,l,j,nz-5)
    !       QJ(i,l,j,2) = QJ(i,l,j,nz-4)
    !       QJ(i,l,j,3) = QJ(i,l,j,nz-3)
    !       QJ(i,l,j,nz-2) = QJ(i,l,j,4)
    !       QJ(i,l,j,nz-1) = QJ(i,l,j,5)
    !       QJ(i,l,j,nz)   = QJ(i,l,j,6)
    ! enddo;enddo;enddo
  end subroutine set_bc

  subroutine set_bc_mut(nx,ny,mut,qc2)
    integer, intent(in), value      :: nx, ny
    real(8), intent(inout), device  :: mut(nx,ny), qc2(nx,ny)
    integer i, j
    !$cuf kernel do(1) <<<*,*>>>
    !do k = 4, nz-3
      do j = 2, ny-1
        ! inlet
        mut(1,j)  = mut(nre2,j)
        qc2(1,j)  = qc2(nre2,j)
        ! outlet
        mut(nx,j) = mut(nx-1,j)
        qc2(nx,j) = qc2(nx-1,j)
      enddo

    !$cuf kernel do(1) <<<*,*>>>
    !do k = 4, nz-3
      do i = 1, nx
        ! wall
        mut(i,1) = 0.d0
        qc2(i,1) = 0.d0
        ! top
        mut(i,ny) = mut(i,ny-1)
        qc2(i,ny) = qc2(i,ny-1)
      enddo

    !cuf kernel do(2) <<<*,*>>>
    ! do j = 1, ny
    !   do i = 1, nx
    !     ! span
    !     mut(i,j,1)    = mut(i,j,nz-5)
    !     mut(i,j,2)    = mut(i,j,nz-4)
    !     mut(i,j,3)    = mut(i,j,nz-3)
    !     mut(i,j,nz-2) = mut(i,j,4)
    !     mut(i,j,nz-1) = mut(i,j,5)
    !     mut(i,j,nz)   = mut(i,j,6)
    !     qc2(i,j,1)    = qc2(i,j,nz-5)
    !     qc2(i,j,2)    = qc2(i,j,nz-4)
    !     qc2(i,j,3)    = qc2(i,j,nz-3)
    !     qc2(i,j,nz-2) = qc2(i,j,4)
    !     qc2(i,j,nz-1) = qc2(i,j,5)
    !     qc2(i,j,nz)   = qc2(i,j,6)
    ! enddo;enddo
  end subroutine set_bc_mut
end module set
