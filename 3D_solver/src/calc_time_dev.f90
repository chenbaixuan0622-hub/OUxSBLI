module calc_time_dev
  use cudafor
  use mpi
  use nvtx
  use mod_globals, only : id_rescale, id_exchange, id_visc, id_LL, id_igr, nt, np, nre2, rerank, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use calc_flux_base
  use calc_steps
  use calc_rescale
  use calc_para
  use calc_rand
  use set
  use preprocess
  use print
  implicit none
  interface RungeKutta
    module procedure RungeKutta_3rd, RungeKutta_4th, Gauss_RungeKutta
  end interface
contains
  subroutine RungeKutta_3rd(id_RungeKutta, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    integer(kind=2), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(5,nx,ny,nz)
    integer i, j, k, l, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireq2(2)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    ! rescal_cpu!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: step = 1, flag_re = 0, flag_req
    real(8), allocatable, device :: Qre(:), Qm(:)
    real(8), allocatable, pinned :: Qm_cpu(:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJ2(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)
    real(8), allocatable, device :: dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)
    ! Landau !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(8), allocatable, device :: seed(:,:,:)
    ! IGR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(4), allocatable, device :: sigma(:,:,:)
    ! for plot
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! count GPU
    stat = cudaGetDeviceCount(ndevices)
    print *, "rank", myrank, " has found ", ndevices, " GPU devices"
    if (mod(myrank,2) == 0) then
      call check_gpu(mygpu)
      call allocate_device_mem(myrank, nx, ny, nz, dx, dy, dz, xix, etay, zetaz, Jacobian, ruvwp, T, mu, mut, qc2, E, F, G)
      allocate(QJ(5,nx,ny,nz), QJ2(5,nx,ny,nz), stat=ierr)
      if (ierr /= 0) then
        print *, "myrank is ", myrank, " memory allocation failed", ierr
      else
        print *, "myrank is ", myrank, " memory allocation has completed"
      endif
      call pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, overlap, &
                    dx, dy, dz, xix, etay, zetaz, Jacobian, QJ, ke0, entropy0)
      if (kind(id_LL) == 4) then
        allocate(seed(nx,ny,nz))
        call init_seed(nx, ny, nz, seed)
      endif
      if (kind(id_igr) == 4) then
        allocate(sigma(nx,ny,nz))
        call set_init_sigma(nx, ny, nz, dx_cpu, dy_cpu, dz_cpu, Jacobian_cpu, Q, sigma)
      endif
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
    endif
    if (kind(id_rescale) == 4) then
      call pre_rescale(myrank, flag_re, flag_req, ny, nz, Qre, Qm, Qm_cpu)
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    print *, "myrank is ", myrank, " start Runge-Kutta"
    do t2 = 1, np
      do t1 = 1, nt
        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(1, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJ, Qm, Qre)
          endif
          call nvtxStartRange("calc flux", 1)
          if (kind(id_LL) == 4) then
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, sigma, seed)
          else
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, sigma)
          endif
          !print *, "myrank is ", myrank, " calc EFG"
          call nvtxEndRange
          call nvtxStartRange("calc step", 2)
          call calc_step1(nx, ny, nz, 1.d0, dx, dy, dz, E, F, G, QJ, QJ2)
          !print *, "myrank is ", myrank, " calc step"
          call nvtxEndRange
          if (ndevices >= 2 .and. kind(id_exchange) == 4) then
            call nvtxStartRange("exchange", 3)
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ2)
            call nvtxEndRange
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ2, Qre)
            call nvtxStartRange("set bc", 4)
          else
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ2)
          endif
          !print *, "myrank is ", myrank, " set bc"
          call nvtxEndRange
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call nvtxStartRange("calc rescale", 5)
          call rescale_recv_send(1, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
          call nvtxEndRange
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(2, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJ2, Qm, Qre)
          endif
          if (kind(id_LL) == 4) then
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G, sigma, seed)
          else
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G, sigma)
          endif
          call calc_step2_3(nx, ny, nz, 0.75d0, 0.25d0, 0.25d0, 1.d0, dx, dy, dz, E, F, G, QJ, QJ2)
          if (ndevices >= 2 .and. kind(id_exchange) == 4) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ2)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ2, Qre)
          else
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ2)
          endif
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(2, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(3, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJ2, Qm, Qre)
          endif
          if (kind(id_LL) == 4) then
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G, sigma, seed)
          else
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G, sigma)
          endif
          call calc_step2_3(nx, ny, nz, 2.d0, 1.d0, 2.d0, 3.d0, dx, dy, dz, E, F, G, QJ2, QJ)
          if (ndevices >= 2 .and. kind(id_exchange) == 4) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ, Qre)
          else
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ)
          endif
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(3, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif
      enddo
      if (mod(myrank, 2) == 0) then
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, QJ, Q, ke0, entropy0)
      else
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
        if (myrank == rerank + 1 .and. kind(id_rescale) == 4) then
          call write_Qm(ny, t2, y, Qm_cpu)
        endif
      endif
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJ2, E, F, G, dx, dy, dz, xix, etay, zetaz, Jacobian)
      if (kind(id_LL) == 4) then
        deallocate(seed)
      endif
      if (kind(id_igr) == 4) then
        deallocate(sigma)
      endif
    endif
    if (kind(id_rescale) == 4) then
      if (myrank == 0 .or. myrank == rerank) then
        deallocate(Qre, Qm)
      elseif (myrank == rerank+1) then
        call write_Qm(ny, np, y, Qm_cpu)
        deallocate(Qm_cpu)
      endif
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_3rd


  subroutine RungeKutta_4th(id_RungeKutta, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    integer(kind=4), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(5,nx,ny,nz)
    integer i, j, k, l, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireq2(2)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    ! rescale !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: step = 1, flag_re = 0, flag_req
    real(8), allocatable, device :: Qre(:), Qm(:)
    real(8), allocatable, pinned :: Qm_cpu(:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)
    real(8), allocatable, device :: dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)
    ! Landau !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(8), allocatable, device :: seed(:,:,:)
    ! IGR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(4), allocatable, device :: sigma(:,:,:)
    ! for plot
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! count GPU
    stat = cudaGetDeviceCount(ndevices)
    if (myrank == 0) then
      print '(2x, i2, a)', ndevices, " GPU devices are found"
    endif
    if (mod(myrank,2) == 0) then
      call check_gpu(mygpu)
      call allocate_device_mem(myrank, nx, ny, nz, dx, dy, dz, xix, etay, zetaz, Jacobian, ruvwp, T, mu, mut, qc2, E, F, G)
      allocate(QJ(5,nx,ny,nz), QJs(5,nx,ny,nz), Rs(5,nx-2,ny-2,nz-2))
      print *, "myrank is ", myrank, " memory allocation has completed"
      call pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, overlap, &
                    dx, dy, dz, xix, etay, zetaz, Jacobian, QJ, ke0, entropy0)
      Rs = 0.d0
      if (kind(id_LL) == 4) then
        allocate(seed(nx,ny,nz))
        call init_seed(nx, ny, nz, seed)
      endif
      if (kind(id_igr) == 4) then
        allocate(sigma(nx,ny,nz))
        call set_init_sigma(nx, ny, nz, dx_cpu, dy_cpu, dz_cpu, Jacobian_cpu, Q, sigma)
      endif
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
    endif
    if (kind(id_rescale) == 4) then
      call pre_rescale(myrank, flag_re, flag_req, ny, nz, Qre, Qm, Qm_cpu)
    endif
    
    do t2 = 1, np
      do t1 = 1, nt
        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(1, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJ, Qm, Qre)
          endif
          if (kind(id_LL) == 4) then
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, sigma, seed)
          else
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, sigma)
          endif
          call calc_step(nx, ny, nz, 0.5d0, 1.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q2
          if (ndevices >= 2) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJs)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs, Qre)
          else
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs)
          endif
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(1, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(2, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJs, Qm, Qre)
          endif
          if (kind(id_LL) == 4) then
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, sigma, seed)
          else
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, sigma)
          endif
          call calc_step(nx, ny, nz, 0.5d0, 2.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q3
          if (ndevices >= 2) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJs)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs, Qre)
          else
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs)
          endif
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(2, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(3, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJs, Qm, Qre)
          endif
          if (kind(id_LL) == 4) then
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, sigma, seed)
          else
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, sigma)
          endif
          call calc_step(nx, ny, nz, 1.0d0, 2.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q4
          if (ndevices >= 2) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJs)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs, Qre)
          else
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs)
          endif
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(3, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(4, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJs, Qm, Qre)
          endif
          if (kind(id_LL) == 4) then
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, sigma, seed)
          else
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, sigma)
          endif
          call calc_step4(nx, ny, nz, dx, dy, dz, E, F, G, Rs, QJ)
          if (ndevices >= 2) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ, Qre)
          else
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ)
          endif
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(4, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif
      enddo
      if (mod(myrank, 2) == 0) then
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, QJ, Q, ke0, entropy0)
      else
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
        if (myrank == rerank + 1 .and. kind(id_rescale) == 4) then
          call write_Qm(ny, t2, y, Qm_cpu)
        endif
      endif
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJs, Rs, E, F, G, dx, dy, dz, xix, etay, zetaz, Jacobian)
      if (kind(id_LL) == 4) then
        deallocate(seed)
      endif
      if (kind(id_igr) == 4) then
        deallocate(sigma)
      endif
    endif
    if (kind(id_rescale) == 4) then
      if (myrank == 0 .or. myrank == rerank) then
        deallocate(Qre, Qm)
      elseif (myrank == rerank+1) then
        call write_Qm(ny, np, y, Qm_cpu)
        deallocate(Qm_cpu)
      endif
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_4th
 

  subroutine Gauss_RungeKutta(id_RungeKutta, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    integer(kind=8), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(5,nx,ny,nz)
    integer i, j, k, itr, max_itr, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireqs(2)
    integer istat(MPI_STATUS_SIZE), istats(MPI_STATUS_SIZE,2)
    ! rescale !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: step, flag_re = 0, flag_req
    real(8) :: c1, c2, a11, a12, a21, a22, b1, b2, err, tol = 1.d-16
    real(8), allocatable, device :: Qre(:), Qm(:)
    real(8), allocatable, pinned :: Qm_cpu(:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)
    real(8), allocatable, device :: R1(:,:,:,:), R2(:,:,:,:), R1_new(:,:,:,:), R2_new(:,:,:,:)
    real(8), allocatable, device :: dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)
    ! Landau !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(8), allocatable, device :: seed(:,:,:)
    ! IGR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(4), allocatable, device :: sigma(:,:,:)
    ! for plot
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! count GPU
    stat = cudaGetDeviceCount(ndevices)
    if (myrank == 0) then
      print '(2x, i2, a)', ndevices, " GPU devices are found"
    endif
    if (mod(myrank,2) == 0) then
      call check_gpu(mygpu)
      call allocate_device_mem(myrank, nx, ny, nz, dx, dy, dz, xix, etay, zetaz, Jacobian, ruvwp, T, mu, mut, qc2, E, F, G)
      allocate(QJ(5,nx,ny,nz), QJs(5,nx,ny,nz), R1(5,nx-2,ny-2,nz-2), R2(5,nx-2,ny-2,nz-2))
      allocate(R1_new(5,nx-2,ny-2,nz-2), R2_new(5,nx-2,ny-2,nz-2))
      print *, "myrank is ", myrank, " memory allocation has completed"
      call pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, overlap, &
                    dx, dy, dz, xix, etay, zetaz, Jacobian, QJ, ke0, entropy0)
      if (kind(id_LL) == 4) then
        allocate(seed(nx,ny,nz))
        call init_seed(nx, ny, nz, seed)
      endif
      if (kind(id_igr) == 4) then
        allocate(sigma(nx,ny,nz))
        call set_init_sigma(nx, ny, nz, dx_cpu, dy_cpu, dz_cpu, Jacobian_cpu, Q, sigma)
      endif
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
    endif
    if (kind(id_rescale) == 4) then
      call pre_rescale(myrank, flag_re, flag_req, ny, nz, Qre, Qm, Qm_cpu)
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    print *, "myrank is ", myrank, " start Runge-Kutta"
    c1  = 0.5d0 - sqrt(3.d0) / 6.d0
    c2  = 0.5d0 - sqrt(3.d0) / 6.d0
    a11 = 0.25d0
    a12 = 0.25d0 - sqrt(3.d0) / 6.d0
    a21 = 0.25d0 + sqrt(3.d0) / 6.d0
    a22 = 0.25d0
    b1  = 0.5d0
    b2  = 0.5d0
    max_itr = 100
    do t2 = 1, np
      do t1 = 1, nt
        if (mod(myrank,2) == 0) then
          ! calc R1
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, sigma)
          call calc_step1(nx, ny, nz, c1, dx, dy, dz, E, F, G, QJ, QJs)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, sigma)
          call calc_R(nx, ny, nz, dx, dy, dz, E, F, G, R1)
          ! calc R2
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, sigma)
          call calc_step1(nx, ny, nz, c2, dx, dy, dz, E, F, G, QJ, QJs)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, sigma)
          call calc_R(nx, ny, nz, dx, dy, dz, E, F, G, R2)
          do itr = 1, max_itr
            ! calc R1
            call calc_Gauss_step(nx, ny, nz, a11, a12, dx, dy, dz, R1, R2, QJ, QJs)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs)
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, sigma)
            call calc_R(nx, ny, nz, dx, dy, dz, E, F, G, R1_new)
            ! calc R2
            call calc_Gauss_step(nx, ny, nz, a21, a22, dx, dy, dz, R1, R2, QJ, QJs)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs)
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, sigma)
            call calc_R(nx, ny, nz, dx, dy, dz, E, F, G, R2_new)
            call calc_error(nx, ny, nz, R1, R2, R1_new, R2_new, err)
            if (err < tol) exit
            R1 = R1_new
            R2 = R2_new
          enddo
          if (err < tol) then
            print *, "Converged at itr=", itr
          else
            print *, "Didn't Converged error=", err
          endif
          call calc_Gauss_step_Q(nx, ny, nz, b1, b2, xix, etay, zetaz, R1, R2, QJ)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ, Qre)
        endif
      enddo
      if (mod(myrank, 2) == 0) then
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, QJ, Q, ke0, entropy0)
      else
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
      endif
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJs, R1, R2, R1_new, R2_new, E, F, G, dx, dy, dz, xix, etay, zetaz, Jacobian)
      if (kind(id_LL) == 4) then
        deallocate(seed)
      endif
      if (kind(id_igr) == 4) then
        deallocate(sigma)
      endif
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine Gauss_RungeKutta
end module calc_time_dev

