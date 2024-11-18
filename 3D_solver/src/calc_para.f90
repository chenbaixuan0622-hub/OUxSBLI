module calc_para
  use mpi
  implicit none
contains
  subroutine flatten(nx, ny, nz, overlap, Q, Q1d_left, Q1d_right)
    integer, intent(in), value   :: nx, ny, nz, overlap
    real(8), intent(in), device  :: Q(nx,ny,nz,5)
    real(8), intent(out), device :: Q1d_left(overlap*(ny-2)*(nz-6)*5)
    real(8), intent(out), device :: Q1d_right(overlap*(ny-2)*(nz-6)*5)
    integer i, j, k, l, ni, nj, nk
    ni = overlap
    nj = ny-2
    nk = nz-6
    !$cuf kernel do(4)<<<*,*>>>
    do l = 1, 5
      do k = 1, nk
        do j = 1, nj
          do i = 1, ni
            Q1d_left(ni*nj*nk*(l-1)+ni*nj*(k-1)+ni*(j-1)+i)  = Q(overlap+i,j+1,k+3,l)
            Q1d_right(ni*nj*nk*(l-1)+ni*nj*(k-1)+ni*(j-1)+i) = Q(nx-2*overlap+i,j+1,k+3,l)
    enddo;enddo;enddo;enddo
  end subroutine flatten

  subroutine reconstruct(nx, ny, nz, overlap, Q1d_left, Q1d_right, Q)
    integer, intent(in), value   :: nx, ny, nz, overlap
    real(8), intent(in), device  :: Q1d_left(overlap*(ny-2)*(nz-6)*5)
    real(8), intent(in), device  :: Q1d_right(overlap*(ny-2)*(nz-6)*5)
    real(8), intent(out), device :: Q(nx,ny,nz,5)
    integer i, j, k, l, ni, nj, nk
    ni = overlap
    nj = ny-2
    nk = nz-6
    !$cuf kernel do(4)<<<*,*>>>
    do l = 1, 5
      do k = 1, nk
        do j = 1, nj
          do i = 1, ni
            Q(i,j+1,k+3,l)            = Q1d_left(ni*nj*nk*(l-1)+ni*nj*(k-1)+ni*(j-1)+i)
            Q(nx-overlap+i,j+1,k+3,l) = Q1d_right(ni*nj*nk*(l-1)+ni*nj*(k-1)+ni*(j-1)+i)
    enddo;enddo;enddo;enddo
  end subroutine reconstruct

  subroutine exchange(myrank, nranks, overlap, nx, ny, nz, QJ, Jacobian)
    integer, intent(in), value     :: myrank, nranks, overlap, nx, ny, nz
    real(8), intent(inout), device :: QJ(nx,ny,nz,5) ! Q / Jacobian
    real(8), intent(in)            :: Jacobian(ny)
    integer rank1, rank2, ierr, ireq4(4), istat(MPI_STATUS_SIZE), istat4(MPI_STATUS_SIZE,4)
    real(8), dimension(overlap*(ny-2)*(nz-6)*5)         :: Qs_left,   Qs_right,   Qr_left,   Qr_right
    real(8), dimension(overlap*(ny-2)*(nz-6)*5), device :: Qs1d_left, Qs1d_right, Qr1d_left, Qr1d_right
    integer j, k, ni, nj, nk

    if (2 <= myrank .and. myrank <= nranks-4) then
      rank1 = myrank-2
      rank2 = myrank+2
    elseif (myrank == 0 .and. 4 <= nranks) then
      rank1 = nranks-2
      rank2 = myrank+2
    elseif (myrank == nranks-2 .and. 4 <= nranks) then
      rank1 = myrank-2
      rank2 = 0
    endif

    call flatten(nx, ny, nz, overlap, QJ, Qs1d_left, Qs1d_right)
    !call MPI_ISEND(Qs1d_left, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, MPI_COMM_WORLD, ireq4(1), ierr)
    !call MPI_ISEND(Qs1d_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, MPI_COMM_WORLD, ireq4(2), ierr)

    !call MPI_IRECV(Qr1d_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, MPI_COMM_WORLD, ireq4(4), ierr)
    !call MPI_IRECV(Qr1d_left, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, MPI_COMM_WORLD, ireq4(3), ierr)
    !call MPI_WAITALL(4, ireq4, istat4, ierr)

    Qs_left    = Qs1d_left
    call MPI_SENDRECV(Qs_left,  5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, &
                      Qr_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, MPI_COMM_WORLD, istat, ierr)
    Qr1d_right = Qr_right

    Qs_right   = Qs1d_right
    call MPI_SENDRECV(Qs_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, &
                      Qr_left,  5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, MPI_COMM_WORLD, istat, ierr)
    Qr1d_left  = Qr_left

    !Qr_left = Qr1d_left

    call reconstruct(nx, ny, nz, overlap, Qr1d_left, Qr1d_right, QJ)
  end subroutine exchange
end module calc_para

