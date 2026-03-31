module calc_time_dev
  use cudafor
  use mpi
  use nvtx
  use mod_globals, only : id_exchange, id_visc, nt, np, nre2, rerank, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use calc_flux_base
  use calc_steps
  use calc_rescale
  use calc_para
  use set
  use preprocess
  use print
  implicit none
  interface RungeKutta
    module procedure RungeKutta_3rd, RungeKutta_3rd_rescale, RungeKutta_4th, RungeKutta_4th_rescale
  end interface
contains 
  subroutine RungeKutta_3rd(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    integer(2), intent(in) :: id_RungeKutta
    integer(2), intent(in) :: id_rescale
    integer, intent(in)    :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)    :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)    :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(5,nx,ny,nz)
    integer i, j, k, l, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireq2(2)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJ2(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)
    real(8), allocatable, device :: dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)
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
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    print *, "myrank is ", myrank, " start Runge-Kutta"
    do t2 = 1, np
      if (mod(myrank,2) == 0) then
        do t1 = 1, nt
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step1<<<blocks,threads>>>(nx, ny, nz, 1.d0, dx, dy, dz, E, F, G, QJ, QJ2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ2)

          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step2_3<<<blocks,threads>>>(nx, ny, nz, 0.75d0, 0.25d0, 0.25d0, 1.d0, dx, dy, dz, E, F, G, QJ, QJ2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ2)

          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step2_3<<<blocks,threads>>>(nx, ny, nz, 2.d0, 1.d0, 2.d0, 3.d0, dx, dy, dz, E, F, G, QJ2, QJ)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ)
        enddo
      endif
      if (mod(myrank, 2) == 0) then
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, QJ, Q, ke0, entropy0)
      else
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
      endif
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJ2, E, F, G, dx, dy, dz, xix, etay, zetaz, Jacobian)
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_3rd


  subroutine RungeKutta_3rd_rescale(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    integer(2), intent(in) :: id_RungeKutta
    integer(4), intent(in) :: id_rescale
    integer, intent(in)    :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)    :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)    :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(5,nx,ny,nz)
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
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
    endif
    call pre_rescale(myrank, flag_re, flag_req, ny, nz, Qre, Qm, Qm_cpu)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    print *, "myrank is ", myrank, " start Runge-Kutta"
    do t2 = 1, np
      do t1 = 1, nt
        if (mod(myrank,2) == 0) then
          call step_rescale(1, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJ, Qm, Qre)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step1<<<blocks,threads>>>(nx, ny, nz, 1.d0, dx, dy, dz, E, F, G, QJ, QJ2)
          call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ2, Qre)
        elseif (myrank == rerank+1) then
          call rescale_recv_send(1, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          call step_rescale(2, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJ2, Qm, Qre)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step2_3<<<blocks,threads>>>(nx, ny, nz, 0.75d0, 0.25d0, 0.25d0, 1.d0, dx, dy, dz, E, F, G, QJ, QJ2)
          call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ2, Qre)
        elseif (myrank == rerank+1) then
          call rescale_recv_send(2, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          call step_rescale(3, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJ2, Qm, Qre)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step2_3<<<blocks,threads>>>(nx, ny, nz, 2.d0, 1.d0, 2.d0, 3.d0, dx, dy, dz, E, F, G, QJ2, QJ)
          call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ, Qre)
        elseif (myrank == rerank+1) then
          call rescale_recv_send(3, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif
      enddo
      if (mod(myrank, 2) == 0) then
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, QJ, Q, ke0, entropy0)
      else
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
      endif
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJ2, E, F, G, dx, dy, dz, xix, etay, zetaz, Jacobian)
    endif
    if (myrank == 0 .or. myrank == rerank) then
      deallocate(Qre, Qm)
    elseif (myrank == rerank+1) then
      deallocate(Qm_cpu)
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_3rd_rescale


  subroutine RungeKutta_4th(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    integer(4), intent(in) :: id_RungeKutta
    integer(2), intent(in) :: id_rescale
    integer, intent(in)    :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)    :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)    :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(5,nx,ny,nz)
    integer i, j, k, l, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireq2(2)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)
    real(8), allocatable, device :: dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)
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
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
    endif
 
    do t2 = 1, np
      if (mod(myrank,2) == 0) then
        do t1 = 1, nt
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step<<<blocks,threads>>>(nx, ny, nz, 0.5d0, 1.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q2
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)

          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step<<<blocks,threads>>>(nx, ny, nz, 0.5d0, 2.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q3
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)

          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step<<<blocks,threads>>>(nx, ny, nz, 1.0d0, 2.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q4
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)

          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step4<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, E, F, G, Rs, QJ)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ)
        enddo
      endif
      if (mod(myrank, 2) == 0) then
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, QJ, Q, ke0, entropy0)
      else
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
      endif
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJs, Rs, E, F, G, dx, dy, dz, xix, etay, zetaz, Jacobian)
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_4th


  subroutine RungeKutta_4th_rescale(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    integer(4), intent(in) :: id_RungeKutta
    integer(4), intent(in) :: id_rescale
    integer, intent(in)    :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)    :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)    :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(5,nx,ny,nz)
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
          call step_rescale(1, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJ, Qm, Qre)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step<<<blocks,threads>>>(nx, ny, nz, 0.5d0, 1.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q2
          call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs, Qre)
        elseif (myrank == rerank+1) then
          call rescale_recv_send(1, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          call step_rescale(2, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJs, Qm, Qre)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step<<<blocks,threads>>>(nx, ny, nz, 0.5d0, 2.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q3
          call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs, Qre)
        elseif (myrank == rerank+1) then
          call rescale_recv_send(2, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          call step_rescale(3, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJs, Qm, Qre)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step<<<blocks,threads>>>(nx, ny, nz, 1.0d0, 2.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q4
          call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs, Qre)
        elseif (myrank == rerank+1) then
          call rescale_recv_send(3, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          call step_rescale(4, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJs, Qm, Qre)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step4<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, E, F, G, Rs, QJ)
          call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ, Qre)
        elseif (myrank == rerank+1) then
          call rescale_recv_send(4, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif
      enddo
      if (mod(myrank, 2) == 0) then
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, QJ, Q, ke0, entropy0)
      else
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
      endif
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJs, Rs, E, F, G, dx, dy, dz, xix, etay, zetaz, Jacobian)
    endif
    if (myrank == 0 .or. myrank == rerank) then
      deallocate(Qre, Qm)
    elseif (myrank == rerank+1) then
      deallocate(Qm_cpu)
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_4th_rescale
end module calc_time_dev

