!> Z-direction MPI halo exchange for curvilinear solver.
!> QJ(nx,5,ny,nz) has the same layout as the Cartesian solver, so pack/unpack
!> kernels are identical to those in 3D_solver/src/calc_para.f90.
module calc_para_curv
  use mpi
  use cudafor
  implicit none
contains

  !> Pack interior slab z=overlap+1..2*overlap into Q1d_lo for sending to rank_lo.
  subroutine flatten_z_lo(nx, ny, nz, overlap, Q, Q1d_lo)
    integer, intent(in), value   :: nx, ny, nz, overlap
    real(8), intent(in), device  :: Q(nx,5,ny,nz)
    real(8), intent(out), device :: Q1d_lo(nx*ny*overlap*5)
    integer i, j, k, l
    !$cuf kernel do(4)<<<*,*>>>
    do k = 1, overlap
      do j = 1, ny
        do l = 1, 5
          do i = 1, nx
            Q1d_lo(ny*nx*5*(k-1)+nx*5*(j-1)+nx*(l-1)+i) = Q(i,l,j,overlap+k)
    enddo;enddo;enddo;enddo
  end subroutine flatten_z_lo


  !> Pack interior slab z=nz-2*overlap+1..nz-overlap into Q1d_hi for sending to rank_hi.
  subroutine flatten_z_hi(nx, ny, nz, overlap, Q, Q1d_hi)
    integer, intent(in), value   :: nx, ny, nz, overlap
    real(8), intent(in), device  :: Q(nx,5,ny,nz)
    real(8), intent(out), device :: Q1d_hi(nx*ny*overlap*5)
    integer i, j, k, l
    !$cuf kernel do(4)<<<*,*>>>
    do k = 1, overlap
      do j = 1, ny
        do l = 1, 5
          do i = 1, nx
            Q1d_hi(ny*nx*5*(k-1)+nx*5*(j-1)+nx*(l-1)+i) = Q(i,l,j,nz-2*overlap+k)
    enddo;enddo;enddo;enddo
  end subroutine flatten_z_hi


  !> Fill ghost cells z=1..overlap from data received from rank_lo.
  subroutine reconstruct_z_lo(nx, ny, nz, overlap, Q1d_lo, Q)
    integer, intent(in), value     :: nx, ny, nz, overlap
    real(8), intent(in), device    :: Q1d_lo(nx*ny*overlap*5)
    real(8), intent(inout), device :: Q(nx,5,ny,nz)
    integer i, j, k, l
    !$cuf kernel do(4)<<<*,*>>>
    do k = 1, overlap
      do j = 1, ny
        do l = 1, 5
          do i = 1, nx
            Q(i,l,j,k) = Q1d_lo(ny*nx*5*(k-1)+nx*5*(j-1)+nx*(l-1)+i)
    enddo;enddo;enddo;enddo
  end subroutine reconstruct_z_lo


  !> Fill ghost cells z=nz-overlap+1..nz from data received from rank_hi.
  subroutine reconstruct_z_hi(nx, ny, nz, overlap, Q1d_hi, Q)
    integer, intent(in), value     :: nx, ny, nz, overlap
    real(8), intent(in), device    :: Q1d_hi(nx*ny*overlap*5)
    real(8), intent(inout), device :: Q(nx,5,ny,nz)
    integer i, j, k, l
    !$cuf kernel do(4)<<<*,*>>>
    do k = 1, overlap
      do j = 1, ny
        do l = 1, 5
          do i = 1, nx
            Q(i,l,j,nz-overlap+k) = Q1d_hi(ny*nx*5*(k-1)+nx*5*(j-1)+nx*(l-1)+i)
    enddo;enddo;enddo;enddo
  end subroutine reconstruct_z_hi


  !> Blocking z ghost-cell exchange using MPI_SENDRECV.
  !> Only even compute ranks call this; step-2 neighbor pattern:
  !>   rank_lo = mod(myrank - 2 + nranks, nranks)
  !>   rank_hi = mod(myrank + 2, nranks)
  subroutine exchange_z_curv(myrank, nranks, overlap, nx, ny, nz, QJ)
    integer, intent(in)              :: myrank, nranks, overlap, nx, ny, nz
    real(8), intent(inout), device   :: QJ(nx,5,ny,nz)
    integer  rank_lo, rank_hi, msglen, stat, ierr
    integer  istat(MPI_STATUS_SIZE)
    real(8), allocatable         :: send_lo(:), send_hi(:), recv_lo(:), recv_hi(:)
    real(8), allocatable, device :: Qs1d_lo(:), Qs1d_hi(:), Qr1d_lo(:), Qr1d_hi(:)
    msglen  = nx * ny * overlap * 5
    rank_lo = mod(myrank - 2 + nranks, nranks)
    rank_hi = mod(myrank + 2, nranks)
    allocate(send_lo(msglen), send_hi(msglen), recv_lo(msglen), recv_hi(msglen))
    allocate(Qs1d_lo(msglen), Qs1d_hi(msglen), Qr1d_lo(msglen), Qr1d_hi(msglen))
    call flatten_z_lo(nx, ny, nz, overlap, QJ, Qs1d_lo)
    call flatten_z_hi(nx, ny, nz, overlap, QJ, Qs1d_hi)
    stat = cudaMemcpy(send_lo, Qs1d_lo, msglen, cudaMemcpyDeviceToHost)
    stat = cudaMemcpy(send_hi, Qs1d_hi, msglen, cudaMemcpyDeviceToHost)
    call MPI_SENDRECV(send_lo, msglen, MPI_REAL8, rank_lo, 30, &
                      recv_hi, msglen, MPI_REAL8, rank_hi, 30, MPI_COMM_WORLD, istat, ierr)
    call MPI_SENDRECV(send_hi, msglen, MPI_REAL8, rank_hi, 31, &
                      recv_lo, msglen, MPI_REAL8, rank_lo, 31, MPI_COMM_WORLD, istat, ierr)
    stat = cudaMemcpy(Qr1d_lo, recv_lo, msglen, cudaMemcpyHostToDevice)
    stat = cudaMemcpy(Qr1d_hi, recv_hi, msglen, cudaMemcpyHostToDevice)
    call reconstruct_z_lo(nx, ny, nz, overlap, Qr1d_lo, QJ)
    call reconstruct_z_hi(nx, ny, nz, overlap, Qr1d_hi, QJ)
    deallocate(send_lo, send_hi, recv_lo, recv_hi, Qs1d_lo, Qs1d_hi, Qr1d_lo, Qr1d_hi)
  end subroutine exchange_z_curv
end module calc_para_curv
