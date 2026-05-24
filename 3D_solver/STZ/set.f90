module set
  use mod_globals, only : nx, ny, nz, gamma, rho_L, p_L, rho_R, p_R, Lx, Ly, Lz
  use set_bc_common
  use set_coordinate
  implicit none
contains

  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    use mpi
    use mod_globals, only : id_accuracy
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    real(8) dx1, dy1, dz1, xf(nx+1), yf(ny+1), zf(nz+1)
    integer nranks, ierr, nz_int, iz_offset, i, j, k
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! interior per rank = nz - 2*overlap_fb; z-rank index = myrank/2
    nz_int    = nz - 2*(kind(id_accuracy)/3 + 1)
    iz_offset = (myrank/2) * nz_int
    dx1 = Lx / dble(nx - 2)
    dy1 = Ly / dble(ny - 2)
    dz1 = Lz / dble((nranks/2) * nz_int)
    dx(:) = dx1;  dy(:) = dy1;  dz(:) = dz1
    ! x: cell-face positions then cell centres
    do i = 2, nx;  xf(i) = dx1 * dble(i-2);  enddo
    xf(1) = xf(2) - dx1;  xf(nx+1) = xf(nx) + dx1
    do i = 1, nx;  x(i) = 0.5d0*(xf(i)+xf(i+1));  enddo
    ! y
    do j = 2, ny;  yf(j) = dy1 * dble(j-2);  enddo
    yf(1) = yf(2) - dy1;  yf(ny+1) = yf(ny) + dy1
    do j = 1, ny;  y(j) = 0.5d0*(yf(j)+yf(j+1));  enddo
    ! z: rank-local slab, with ghost face positions
    do k = 2, nz;  zf(k) = dz1 * dble(iz_offset + k - 2);  enddo
    zf(1) = zf(2) - dz1;  zf(nz+1) = zf(nz) + dz1
    do k = 1, nz;  z(k) = 0.5d0*(zf(k)+zf(k+1));  enddo
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,5,ny,nz)
    integer i, j, k
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (myrank < 2) then
            Q(i,1,j,k) = rho_L
            Q(i,5,j,k) = p_L / (gamma - 1.d0)
          else
            Q(i,1,j,k) = rho_R
            Q(i,5,j,k) = p_R / (gamma - 1.d0)
          endif
          Q(i,2,j,k) = 0.d0   ! rho*u
          Q(i,3,j,k) = 0.d0   ! rho*v
          Q(i,4,j,k) = 0.d0   ! rho*w
        enddo
      enddo
    enddo
  end subroutine set_init


  ! Custom BC for z-decomposition:
  !   - periodic in x and y (for all z including ghost layers)
  !   - zero-gradient in z only for outermost ranks
  !   - does NOT call set_bc_cyclic (which would corrupt z ghost cells)
  subroutine set_bc(myrank, nx, ny, nz, Jacobian, QJ, Qre)
    use mpi
    use mod_globals, only : id_accuracy
    integer, intent(in), value               :: myrank, nx, ny, nz
    real(8), intent(in), device              :: Jacobian(nx,ny)
    real(8), intent(inout), device           :: QJ(nx,5,ny,nz)
    real(8), intent(in), device, optional    :: Qre(ny*(nz-6)*5)
    integer nranks, ierr, i, j, k, l, m, ovlp
    ovlp = kind(id_accuracy)/3 + 1
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! x: Dirichlet slip wall — reflect ρu (l=2) to enforce u=0 at face; extrapolate rest
    !$cuf kernel do(3)<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do m = 0, ovlp-1
          QJ(ovlp-m,      1,j,k) =  QJ(ovlp+1+m,  1,j,k)
          QJ(ovlp-m,      2,j,k) = -QJ(ovlp+1+m,  2,j,k)
          QJ(ovlp-m,      3,j,k) =  QJ(ovlp+1+m,  3,j,k)
          QJ(ovlp-m,      4,j,k) =  QJ(ovlp+1+m,  4,j,k)
          QJ(ovlp-m,      5,j,k) =  QJ(ovlp+1+m,  5,j,k)
          QJ(nx-ovlp+1+m, 1,j,k) =  QJ(nx-ovlp-m, 1,j,k)
          QJ(nx-ovlp+1+m, 2,j,k) = -QJ(nx-ovlp-m, 2,j,k)
          QJ(nx-ovlp+1+m, 3,j,k) =  QJ(nx-ovlp-m, 3,j,k)
          QJ(nx-ovlp+1+m, 4,j,k) =  QJ(nx-ovlp-m, 4,j,k)
          QJ(nx-ovlp+1+m, 5,j,k) =  QJ(nx-ovlp-m, 5,j,k)
        enddo
      enddo
    enddo
    ! y: Dirichlet slip wall — reflect ρv (l=3) to enforce v=0 at face; extrapolate rest
    !$cuf kernel do(3)<<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        do m = 0, ovlp-1
          QJ(i,1,ovlp-m,     k) =  QJ(i,1,ovlp+1+m,  k)
          QJ(i,2,ovlp-m,     k) =  QJ(i,2,ovlp+1+m,  k)
          QJ(i,3,ovlp-m,     k) = -QJ(i,3,ovlp+1+m,  k)
          QJ(i,4,ovlp-m,     k) =  QJ(i,4,ovlp+1+m,  k)
          QJ(i,5,ovlp-m,     k) =  QJ(i,5,ovlp+1+m,  k)
          QJ(i,1,ny-ovlp+1+m,k) =  QJ(i,1,ny-ovlp-m, k)
          QJ(i,2,ny-ovlp+1+m,k) =  QJ(i,2,ny-ovlp-m, k)
          QJ(i,3,ny-ovlp+1+m,k) = -QJ(i,3,ny-ovlp-m, k)
          QJ(i,4,ny-ovlp+1+m,k) =  QJ(i,4,ny-ovlp-m, k)
          QJ(i,5,ny-ovlp+1+m,k) =  QJ(i,5,ny-ovlp-m, k)
        enddo
      enddo
    enddo
    ! z lo boundary: zero-gradient on rank 0 — set all ovlp ghost cells
    if (myrank == 0) then
      !$cuf kernel do(3)<<<*,*>>>
      do k = 1, ovlp
        do j = 1, ny
          do i = 1, nx
            do l = 1, 5
              QJ(i,l,j,k) = QJ(i,l,j,ovlp+1)
            enddo
          enddo
        enddo
      enddo
    endif
    ! z hi boundary: zero-gradient on last compute rank (nranks-2 in even/odd pattern)
    if (myrank == nranks - 2) then
      !$cuf kernel do(3)<<<*,*>>>
      do k = 1, ovlp
        do j = 1, ny
          do i = 1, nx
            do l = 1, 5
              QJ(i,l,j,nz-ovlp+k) = QJ(i,l,j,nz-ovlp)
            enddo
          enddo
        enddo
      enddo
    endif
  end subroutine set_bc


  subroutine set_bc_mut(nx, ny, nz, mut, qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    call set_bc_mut_common(nx, ny, nz, mut, qc2)
  end subroutine set_bc_mut
end module set
