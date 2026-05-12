!> Module for preprocessing and GPU memory management
!> Handles device setup, memory allocation, and initial data transfers
module preprocess
  use cudafor
  use mpi
  use print
  implicit none
contains

  !> Check GPU device properties and availability
  !> Prints device name and capability information
  subroutine check_gpu(mygpu)
    integer, intent(in) :: mygpu                         !< GPU device ID to check
    integer ilen, stat
    type(cudaDeviceProp) prop
    stat = cudaSetDevice(mygpu)
    stat = cudaGetDeviceProperties(prop, mygpu)
    ilen = verify(prop%name, ' ', .true.)
    print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", mygpu, ") is available"
  end subroutine check_gpu


  !> Allocate GPU device memory for simulation variables
  !> Size and allocation depends on viscosity model selection
  subroutine allocate_device_mem(myrank, nx, ny, dtdxdy,  xix, etay,  Jacobian, ruvwp, T, mu, mut, qc2, E, F)
    use mod_globals, only : id_visc
    use calc_flux_base, only : init_sensor
    integer, intent(in)                       :: myrank    !< MPI rank
    integer, intent(in)                       :: nx        !< x grid dimension
    integer, intent(in)                       :: ny        !< y grid dimension
   ! integer, intent(in)                       :: nz        !< z grid dimension
    real(8), intent(out), allocatable, device :: dtdxdy(:,:) !< dt * Sxy
    !real(8), intent(out), allocatable, device :: dtdydz(:,:) !< dt * Syz
    !real(8), intent(out), allocatable, device :: dtdzdx(:,:) !< dt * Szx
    real(8), intent(out), allocatable, device :: xix(:)    !< coordinate transform metric in x
    real(8), intent(out), allocatable, device :: etay(:)   !< coordinate transform metric in y
   ! real(8), intent(out), allocatable, device :: zetaz(:)  !< coordinate transform metric in z
    real(8), intent(out), allocatable, device :: Jacobian(:,:)  !< Jacobian determinant for coordinate transform
    real(8), intent(out), allocatable, device :: ruvwp(:,:,:) !< work array for momentum/velocities
    real(8), intent(out), allocatable, device :: T(:,:)   !< temperature field
    real(8), intent(out), allocatable, device :: mu(:,:)  !< molecular viscosity
    real(8), intent(out), allocatable, device :: mut(:,:) !< turbulent viscosity (LES)
    real(8), intent(out), allocatable, device :: qc2(:,:) !< quadratic constitutive terms
    real(8), intent(out), allocatable, device :: E(:,:,:) !< flux in x direction
    real(8), intent(out), allocatable, device :: F(:,:,:) !< flux in y direction
    !real(8), intent(out), allocatable, device :: G(:,:,:,:) !< flux in z direction
    integer ierr
    allocate(ruvwp(4,nx,ny), E(4,nx-1,ny-2), F(4,nx-2,ny-1),  stat=ierr)
    allocate(dtdxdy(nx-2,ny-2),  xix(nx-1), etay(ny-1), Jacobian(nx,ny), stat=ierr)
    if (kind(id_visc) == 2) then
      allocate(T(nx,ny), mu(1,1), mut(1,1), qc2(1,1), stat=ierr)
    elseif (kind(id_visc) == 4) then
      allocate(T(nx,ny), mu(nx,ny), mut(1,1), qc2(1,1), stat=ierr)
    elseif (kind(id_visc) == 8) then
      allocate(T(nx,ny), mu(nx,ny), mut(nx,ny), qc2(nx,ny), stat=ierr)
    endif
    call init_sensor(nx, ny)
    if (ierr /= 0) then
      print *, "myrank is ", myrank, " memory allocation failed", ierr
    else
      print *, "myrank is ", myrank, " memory allocation has completed"
    endif
  end subroutine allocate_device_mem


  !> Preprocessing: compute metrics, initialize Q, transfer to device
  !> Divides computational domain across MPI ranks
  subroutine pre_calc(nx, ny, myrank, nranks, x, dx_cpu, y, dy_cpu,  Jacobian_cpu, Q, overlap, &
                      dtdxdy, xix, etay, Jacobian, QJ, ke0, entropy0)
    use mod_globals, only : dt
    integer, intent(in)                      :: nx                  !< x grid dimension
    integer, intent(in)                      :: ny                  !< y grid dimension
    !integer, intent(in)                      :: nz                  !< z grid dimension
    integer, intent(in)                      :: myrank              !< MPI rank of this process
    integer, intent(in)                      :: nranks              !< total number of MPI ranks
    real(8), intent(in)                      :: x(nx)               !< x coordinate array (host)
    real(8), intent(in)                      :: dx_cpu(nx-1)        !< inverse x spacing (host)
    real(8), intent(in)                      :: y(ny)               !< y coordinate array (host)
    real(8), intent(in)                      :: dy_cpu(ny-1)        !< inverse y spacing (host)
    !real(8), intent(in)                      :: z(nz)               !< z coordinate array (host)
   ! real(8), intent(in)                      :: dz_cpu(nz-1)        !< inverse z spacing (host)
    real(8), intent(in)                      :: Jacobian_cpu(nx,ny) !< Jacobian determinant (host)
    real(8), intent(inout)                   :: Q(nx,4,ny)       !< conservative variables on host
    integer, intent(out)                     :: overlap             !< ghost cell width for MPI halo exchange
    real(8), intent(out), device, contiguous :: dtdxdy(nx-2,ny-2)   !< dt * Sxy (device)
    !real(8), intent(out), device, contiguous :: dtdydz(ny-2,nz-2)   !< dt * Syz (device)
    !real(8), intent(out), device, contiguous :: dtdzdx(nx-2,nz-2)   !< dt * Szx (device)
    real(8), intent(out), device, contiguous :: xix(nx-1)           !< x coordinate metric (device)
    real(8), intent(out), device, contiguous :: etay(ny-1)          !< y coordinate metric (device)
    !real(8), intent(out), device, contiguous :: zetaz(nz-1)         !< z coordinate metric (device)
    real(8), intent(out), device, contiguous :: Jacobian(nx,ny)     !< Jacobian determinant (device)
    real(8), intent(out), device, contiguous :: QJ(nx,4,ny)      !< Q divided by Jacobian (device)
    real(4), intent(inout)                   :: ke0                 !< reference kinetic energy
    real(4), intent(inout)                   :: entropy0            !< reference entropy
    real(8) xix_cpu(nx-1), etay_cpu(ny-1)
    real(8) dtdxdy_cpu(nx-2,ny-2)
    real(4) rho1d(nx*ny), p1d(nx*ny), v1d(nx*ny*3)
    integer i, j, k, l, ierr
    ! set Q / Jacobian
    !do k = 1, nz
      do j = 1, ny
        do l = 1, 5
          do i = 1, nx
            Q(i,l,j) = Q(i,l,j) / Jacobian_cpu(i,j)
          enddo
        enddo
      enddo
    !enddo
    ! copy on GPU
    xix_cpu   = 1.d0 / dx_cpu
    etay_cpu  = 1.d0 / dy_cpu
    !zetaz_cpu = 1.d0 / dz_cpu
    do j = 1, ny-2
      do i = 1, nx-2
        dtdxdy_cpu(i,j) = 0.25d0 * dt * (dx_cpu(i) + dx_cpu(i+1)) * (dy_cpu(j) + dy_cpu(j+1))
    enddo;enddo
    !do k = 1, nz-2
     ! do j = 1, ny-2
       ! dtdydz_cpu(j,k) = 0.25d0 * dt * (dy_cpu(j) + dy_cpu(j+1)) * (dz_cpu(k) + dz_cpu(k+1))
    !enddo;enddo
    !do k = 1, nz-2
      !do i = 1, nx-2
      !  dtdzdx_cpu(i,k) = 0.25d0 * dt * (dz_cpu(k) + dz_cpu(k+1)) * (dx_cpu(i) + dx_cpu(i+1))
    !enddo;enddo
    dtdxdy = dtdxdy_cpu
    !dtdydz = dtdydz_cpu
    !dtdzdx = dtdzdx_cpu
    xix      = xix_cpu
    etay     = etay_cpu
    !zetaz    = zetaz_cpu
    Jacobian = Jacobian_cpu
    QJ = Q
    ! for multi GPU
    if (kind(id_accuracy) == 8) then
      overlap = 3
    elseif (kind(id_accuracy) == 4) then
      overlap = 2
    else
      overlap = 1
    endif
    call make_1d_for_print(nx, ny, Jacobian_cpu, Q, rho1d, p1d, v1d)
    call print_vtk(0, nx, ny, x, y, rho1d, p1d, v1d)     !(0, nx, ny,  myrank+1, nranks, x, y, rho1d, p1d, v1d, ke0, entropy0)
    call MPI_SEND(ke0,      1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
    call MPI_SEND(entropy0, 1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
  end subroutine pre_calc


  subroutine pre_rescale(myrank, flag_re, flag_req, ny, Qre, Qm, Qm_cpu)
    use mod_globals, only : rerank, id_recal, id_rescale
    integer, intent(in)                         :: myrank, ny
    integer, intent(inout)                      :: flag_re, flag_req
    real(8), intent(inout), allocatable, device :: Qre(:), Qm(:)
    real(8), intent(inout), allocatable         :: Qm_cpu(:)
    character(len=40) filename
    logical exists
    integer j, stat, ilen, ireq, ierr, errorcode
    type(cudaDeviceProp) prop
    if (myrank == rerank) then
      call MPI_IRECV(flag_re, 1, MPI_INTEGER, rerank+1, 1001, MPI_COMM_WORLD, flag_req, ierr)
    endif
    if (mod(myrank,2) == 0) then
      allocate(Qre(ny*5), Qm(ny*5), stat=ierr)
      Qm(:) = 0.d0
      if (ierr /= 0) then
        print *, "myrank is ", myrank, " memory allocation failed (Qm)", ierr
      else
        print *, "myrank is ", myrank, " memory allocation has completed (Qm)"
      endif
    elseif (myrank == rerank+1) then
      stat = cudaSetDevice(0)
      stat = cudaGetDeviceProperties(prop, 0)
      ilen = verify(prop%name, ' ', .true.)
      print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", 0, ") calculates rescaling"
      allocate(Qm_cpu(ny*5))
      if (kind(id_recal) == 4) then
        write(filename, "(a)") "recal/Qm.dat"
        inquire(file=filename, exist=exists)
        if (exists) then
          if (id_rescale == 0) then
            print *, "Qm exists (not fixed value)"
          else
            print *, "Qm exists (fixed value)"
          endif
          !open(10, file=filename, action="read", form="unformatted", access="stream", status="old")
          !read(10) Qm_cpu
          !close(10)
          do j = 1, 5*ny
            if (Qm_cpu(j) /= Qm_cpu(j)) then
              print *, "Qm is NaN"
              call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
            endif
          enddo
        else
          if (id_rescale /= 0) then
            print *, "Qm is required but doesn't exist"
            call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
          endif
          Qm_cpu(:) = 0.d0
          print *, "Qm doesn't exist"
        endif
      endif
    endif
  end subroutine pre_rescale
end module preprocess

