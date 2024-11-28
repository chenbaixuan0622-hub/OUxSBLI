module calc_para
  use mpi
  implicit none
  interface exchange
    module procedure exchange_cyclic, exchange_rescale
  end interface exchange
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

  subroutine flatten_left(nx, ny, nz, overlap, Q, Q1d_left)
    integer, intent(in), value   :: nx, ny, nz, overlap
    real(8), intent(in), device  :: Q(nx,ny,nz,5)
    real(8), intent(out), device :: Q1d_left(overlap*(ny-2)*(nz-6)*5)
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
    enddo;enddo;enddo;enddo
  end subroutine flatten_left

  subroutine flatten_right(nx, ny, nz, overlap, Q, Q1d_right)
    integer, intent(in), value   :: nx, ny, nz, overlap
    real(8), intent(in), device  :: Q(nx,ny,nz,5)
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
            Q1d_right(ni*nj*nk*(l-1)+ni*nj*(k-1)+ni*(j-1)+i) = Q(nx-2*overlap+i,j+1,k+3,l)
    enddo;enddo;enddo;enddo
  end subroutine flatten_right
  
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

  subroutine reconstruct_left(nx, ny, nz, overlap, Q1d_left, Q)
    integer, intent(in), value   :: nx, ny, nz, overlap
    real(8), intent(in), device  :: Q1d_left(overlap*(ny-2)*(nz-6)*5)
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
            Q(i,j+1,k+3,l) = Q1d_left(ni*nj*nk*(l-1)+ni*nj*(k-1)+ni*(j-1)+i)
    enddo;enddo;enddo;enddo
  end subroutine reconstruct_left
  
  subroutine reconstruct_right(nx, ny, nz, overlap, Q1d_right, Q)
    integer, intent(in), value   :: nx, ny, nz, overlap
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
            Q(nx-overlap+i,j+1,k+3,l) = Q1d_right(ni*nj*nk*(l-1)+ni*nj*(k-1)+ni*(j-1)+i)
    enddo;enddo;enddo;enddo
  end subroutine reconstruct_right

  subroutine exchange_cyclic(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ)
    integer(kind=2), intent(in), value :: id_rescale
    integer, intent(in), value         :: myrank, nranks, overlap, nx, ny, nz
    real(8), intent(inout), device     :: QJ(nx,ny,nz,5) ! Q / Jacobian
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

    Qs_left    = Qs1d_left
    call MPI_SENDRECV(Qs_left,  5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, &
                      Qr_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, MPI_COMM_WORLD, istat, ierr)
    Qr1d_right = Qr_right

    Qs_right   = Qs1d_right
    call MPI_SENDRECV(Qs_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, &
                      Qr_left,  5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, MPI_COMM_WORLD, istat, ierr)
    Qr1d_left  = Qr_left

    call reconstruct(nx, ny, nz, overlap, Qr1d_left, Qr1d_right, QJ)
  end subroutine exchange_cyclic

  subroutine exchange_rescale(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ)
    integer(kind=4), intent(in), value :: id_rescale
    integer, intent(in), value         :: myrank, nranks, overlap, nx, ny, nz
    real(8), intent(inout), device     :: QJ(nx,ny,nz,5) ! Q / Jacobian
    integer rank1, rank2, ierr, ireq4(4), istat(MPI_STATUS_SIZE), istat4(MPI_STATUS_SIZE,4)
    real(8), dimension(overlap*(ny-2)*(nz-6)*5)         :: Qs_left,   Qs_right,   Qr_left,   Qr_right
    real(8), dimension(overlap*(ny-2)*(nz-6)*5), device :: Qs1d_left, Qs1d_right, Qr1d_left, Qr1d_right
    integer j, k, ni, nj, nk

    if (2 <= myrank .and. myrank <= nranks-4) then
      rank1 = myrank-2
      rank2 = myrank+2

      call flatten(nx, ny, nz, overlap, QJ, Qs1d_left, Qs1d_right)

      Qs_left    = Qs1d_left
      call MPI_SENDRECV(Qs_left,  5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, &
                        Qr_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, MPI_COMM_WORLD, istat, ierr)
      Qr1d_right = Qr_right

      Qs_right   = Qs1d_right
      call MPI_SENDRECV(Qs_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, &
                        Qr_left,  5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, MPI_COMM_WORLD, istat, ierr)
      Qr1d_left  = Qr_left

      call reconstruct(nx, ny, nz, overlap, Qr1d_left, Qr1d_right, QJ)
    elseif (myrank == 0 .and. 4 <= nranks) then
      rank2 = myrank+2

      call flatten_right(nx, ny, nz, overlap, QJ, Qs1d_right)

      call MPI_RECV(Qr_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, MPI_COMM_WORLD, istat, ierr)
      Qr1d_right = Qr_right

      Qs_right   = Qs1d_right
      call MPI_SEND(Qs_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, MPI_COMM_WORLD, ierr)

      call reconstruct_right(nx, ny, nz, overlap, Qr1d_right, QJ)
    elseif (myrank == nranks-2 .and. 4 <= nranks) then
      rank1 = myrank-2

      call flatten_left(nx, ny, nz, overlap, QJ, Qs1d_left)

      Qs_left    = Qs1d_left
      call MPI_SEND(Qs_left, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, MPI_COMM_WORLD, ierr)

      call MPI_RECV(Qr_left, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, MPI_COMM_WORLD, istat, ierr)
      Qr1d_left  = Qr_left

      call reconstruct_left(nx, ny, nz, overlap, Qr1d_left, QJ)
    endif
  end subroutine exchange_rescale
end module calc_para

