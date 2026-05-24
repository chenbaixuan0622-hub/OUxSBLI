!> Module for time-stepping orchestration and time integration schemes
!> Implements 3rd-order and 4th-order Runge-Kutta time stepping with MPI/GPU support
module calc_time_dev
  use cudafor
  use mpi
  use nvtx
  use mod_globals, only : id_visc, nt, np, nre2, rerank, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use mod_constant, only : one_third
  use calc_flux_base
  use calc_steps
  use calc_rescale
  use calc_para
  use set
  use preprocess
  use print
  implicit none
  interface RungeKutta
    module procedure RungeKutta_3rd, RungeKutta_3rd_rescale, RungeKutta_4th, RungeKutta_4th_rescale, &
                     RungeKutta_3rd_zdec, RungeKutta_4th_zdec
  end interface
contains 

  !> 3rd-order TVD Runge-Kutta time stepping without rescaling
  !> Solves dQ/dt = RHS(Q) using total variation diminishing (TVD) 3-stage scheme
  !> Stages: Q(1) = Q^n + (dt)*RHS(Q^n); Q(2) = (3/4)*Q^n + (1/4)*Q(1) + (1/4)*dt*RHS(Q(1))
  !>         Q^(n+1) = (1/3)*Q^n + (2/3)*Q(2) + (2/3)*dt*RHS(Q(2))
  !> TVD property: |Q^(n+1)|_TV <= |Q^n|_TV prevents spurious oscillations near shocks
  !> GPU computation: Each rank manages one GPU asynchronously; MPI sync only for I/O
  subroutine RungeKutta_3rd(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    integer(2), intent(in) :: id_RungeKutta                     !< time integration method ID
    integer(2), intent(in) :: id_rescale                        !< rescaling method ID
    integer, intent(in)    :: myrank                            !< MPI rank
    integer, intent(in)    :: mygpu                             !< GPU index for this rank
    integer, intent(in)    :: nx                                !< x grid dimension
    integer, intent(in)    :: ny                                !< y grid dimension
    integer, intent(in)    :: nz                                !< z grid dimension
    real(8), intent(in)    :: x(nx)                             !< x coordinates
    real(8), intent(in)    :: dx_cpu(nx-1)                      !< inverse x spacing (host)
    real(8), intent(in)    :: y(ny)                             !< y coordinates
    real(8), intent(in)    :: dy_cpu(ny-1)                      !< inverse y spacing (host)
    real(8), intent(in)    :: z(nz)                             !< z coordinates
    real(8), intent(in)    :: dz_cpu(nz-1)                      !< inverse z spacing (host)
    real(8), intent(in)    :: Jacobian_cpu(nx,ny)               !< Jacobian determinant (host)
    real(8), intent(inout) :: Q(nx,5,ny,nz)                     !< conservative variables on host
    integer i, j, k, l, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireq2(2)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJ2(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), zetaz(:), Jacobian(:,:), dtdxdy(:,:), dtdydz(:,:), dtdzdx(:,:)
    ! for plot
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! count GPU
    stat = cudaGetDeviceCount(ndevices)
    print *, "rank", myrank, " has found ", ndevices, " GPU devices"
    if (mod(myrank,2) == 0) then
      call check_gpu(mygpu)
      call allocate_device_mem(myrank, nx, ny, nz, dtdxdy, dtdydz, dtdzdx, xix, etay, zetaz, Jacobian, ruvwp, T, mu, mut, qc2, E, F, G)
      allocate(QJ(nx,5,ny,nz), QJ2(nx,5,ny,nz), stat=ierr)
      if (ierr /= 0) then
        print *, "myrank is ", myrank, " memory allocation failed", ierr
      else
        print *, "myrank is ", myrank, " memory allocation has completed"
      endif
      call pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, overlap, &
                    dtdxdy, dtdydz, dtdzdx, xix, etay, zetaz, Jacobian, QJ, ke0, entropy0, .false.)
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
    endif

    !> Time integration main loop with MPI domain decomposition
    !> Each rank computes assigned 3D subdomain on dedicated GPU independently
    !> Synchronization points: (1) end of each time integration step for boundary exchange
    !>                         (2) after every np iterations for I/O and statistics
    !> CFL constraint: dt = CFL * min_grid_spacing / max_wave_speed
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    print *, "myrank is ", myrank, " start Runge-Kutta"
    do t2 = 1, np
      if (mod(myrank,2) == 0) then
        ! GPU-accelerated ranks perform time integration
        do t1 = 1, nt
          ! Step 1: Compute fluxes E, F, G from current state QJ
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, 1, nz-1)
          ! Step 2a: TVD RK3 Stage 1 - compute Q(1), store in QJ2
          call calc_step1<<<blocks,threads>>>(nx, ny, nz, 1.d0, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ, QJ2)
          ! Enforce boundary conditions at cell interfaces (extrapolation or characteristic-based)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ2)

          ! Step 2b: TVD RK3 Stage 2 - blend Q(1) with Q^n, store in QJ2
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G, 1, nz-1)
          call calc_step2_3<<<blocks,threads>>>(nx, ny, nz, 0.75d0, 0.25d0, 0.25d0, 1.d0, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ, QJ2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ2)

          ! Step 2c: TVD RK3 Stage 3 - final solution Q^(n+1), store in QJ (swap arrays)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G, 1, nz-1)
          call calc_step2_3<<<blocks,threads>>>(nx, ny, nz, 2.d0, 1.d0, 2.d0, one_third, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ2, QJ)
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
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJ2, E, F, G, xix, etay, zetaz, Jacobian, dtdxdy, dtdydz, dtdzdx)
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_3rd


  !> 3rd-order TVD Runge-Kutta with re-scaling for density clipping/limiting
  !> Combines TVD RK3 time stepping with optional density-based re-scaling for stability
  !> When density becomes negative or too small, re-scale conserved variables at marker plane
  !> Similar TVD RK3 stages as above, plus calls to step_rescale() for non-conservative correction
  subroutine RungeKutta_3rd_rescale(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    integer(2), intent(in) :: id_RungeKutta
    integer(4), intent(in) :: id_rescale
    integer, intent(in)    :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)    :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)    :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(nx,5,ny,nz)
    integer i, j, k, l, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireq2(2)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    ! rescal_cpu!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: step = 1, flag_re = 0, flag_req
    real(8), allocatable, device :: Qre(:), Qm(:)
    real(8), allocatable, pinned :: Qm_cpu(:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJ2(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), zetaz(:), Jacobian(:,:), dtdxdy(:,:), dtdydz(:,:), dtdzdx(:,:)
    ! for plot
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! count GPU
    stat = cudaGetDeviceCount(ndevices)
    print *, "rank", myrank, " has found ", ndevices, " GPU devices"
    if (mod(myrank,2) == 0) then
      call check_gpu(mygpu)
      call allocate_device_mem(myrank, nx, ny, nz, dtdxdy, dtdydz, dtdzdx, xix, etay, zetaz, Jacobian, ruvwp, T, mu, mut, qc2, E, F, G)
      allocate(QJ(nx,5,ny,nz), QJ2(nx,5,ny,nz), stat=ierr)
      if (ierr /= 0) then
        print *, "myrank is ", myrank, " memory allocation failed", ierr
      else
        print *, "myrank is ", myrank, " memory allocation has completed"
      endif
      call pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, overlap, &
                    dtdxdy, dtdydz, dtdzdx, xix, etay, zetaz, Jacobian, QJ, ke0, entropy0, .false.)
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
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, 1, nz-1)
          call calc_step1<<<blocks,threads>>>(nx, ny, nz, 1.d0, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ, QJ2)
          call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ2, Qre)
        elseif (myrank == rerank+1) then
          call rescale_recv_send(1, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          call step_rescale(2, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJ2, Qm, Qre)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G, 1, nz-1)
          call calc_step2_3<<<blocks,threads>>>(nx, ny, nz, 0.75d0, 0.25d0, 0.25d0, 1.d0, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ, QJ2)
          call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ2, Qre)
        elseif (myrank == rerank+1) then
          call rescale_recv_send(2, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          call step_rescale(3, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJ2, Qm, Qre)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G, 1, nz-1)
          call calc_step2_3<<<blocks,threads>>>(nx, ny, nz, 2.d0, 1.d0, 2.d0, one_third, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ2, QJ)
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
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJ2, E, F, G, xix, etay, zetaz, Jacobian, dtdxdy, dtdydz, dtdzdx)
    endif
    if (myrank == 0 .or. myrank == rerank) then
      deallocate(Qre, Qm)
    elseif (myrank == rerank+1) then
      deallocate(Qm_cpu)
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_3rd_rescale


  !> 4th-order classical Runge-Kutta time stepping without rescaling
  !> Solves dQ/dt = RHS(Q) using standard 4-stage RK4 scheme
  !> Stages: k1 = RHS(Q^n); k2 = RHS(Q^n + 0.5*dt*k1); k3 = RHS(Q^n + 0.5*dt*k2)
  !>         k4 = RHS(Q^n + dt*k3); Q^(n+1) = Q^n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
  !> More accurate than RK3 but lacks TVD property; requires smaller CFL (~0.8 vs 1.0)
  subroutine RungeKutta_4th(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    integer(4), intent(in) :: id_RungeKutta
    integer(2), intent(in) :: id_rescale
    integer, intent(in)    :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)    :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)    :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(nx,5,ny,nz)
    integer i, j, k, l, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireq2(2)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), zetaz(:), Jacobian(:,:), dtdxdy(:,:), dtdydz(:,:), dtdzdx(:,:)
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
      call allocate_device_mem(myrank, nx, ny, nz, dtdxdy, dtdydz, dtdzdx, xix, etay, zetaz, Jacobian, ruvwp, T, mu, mut, qc2, E, F, G)
      allocate(QJ(nx,5,ny,nz), QJs(nx,5,ny,nz), Rs(nx-2,5,ny-2,nz-2))
      print *, "myrank is ", myrank, " memory allocation has completed"
      call pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, overlap, &
                    dtdxdy, dtdydz, dtdzdx, xix, etay, zetaz, Jacobian, QJ, ke0, entropy0, .false.)
      Rs = 0.d0
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
    endif
 
    !> 4-4 RK main integration loop with residual accumulation
    !> Stages 1-3: compute fluxes and update intermediate solutions, accumulate residuals in Rs
    !> Stage 4: final flux computation and assembly of weighted sum Q^(n+1) = Q^n - (1/6)*sum(R_i)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    print *, "myrank is ", myrank, " start Runge-Kutta"
    do t2 = 1, np
      if (mod(myrank,2) == 0) then
        do t1 = 1, nt
          ! Stage 1: k1 = RHS(Q^n), coefficients: 0.5*dt applied, weight 1.0 to Rs
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, 1, nz-1)
          call calc_step<<<blocks,threads>>>(nx, ny, nz, 0.5d0, 1.d0, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ, QJs, Rs) ! QJs = Q2 = Q^n + 0.5*dt*k1
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)

          ! Stage 2: k2 = RHS(Q^n + 0.5*dt*k1), weight 2.0 to Rs for (2*k2 term)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, 1, nz-1)
          call calc_step<<<blocks,threads>>>(nx, ny, nz, 0.5d0, 2.d0, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ, QJs, Rs) ! QJs = Q3 = Q^n + 0.5*dt*k2
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)

          ! Stage 3: k3 = RHS(Q^n + 0.5*dt*k2), weight 2.0 for (2*k3 term)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, 1, nz-1)
          call calc_step<<<blocks,threads>>>(nx, ny, nz, 1.0d0, 2.d0, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ, QJs, Rs) ! QJs = Q4 = Q^n + dt*k3
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)

          ! Stage 4: k4 = RHS(Q^n + dt*k3), weight 1.0, assemble final Q^(n+1)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, 1, nz-1)
          call calc_step4<<<blocks,threads>>>(nx, ny, nz, dtdxdy, dtdydz, dtdzdx, E, F, G, Rs, QJ)  ! QJ = Q^(n+1)
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
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJs, Rs, E, F, G, xix, etay, zetaz, Jacobian, dtdxdy, dtdydz, dtdzdx)
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_4th


  !> 4th-order Runge-Kutta with re-scaling for density clipping at monitoring plane
  !> Combines high-order 4-4 RK accuracy with non-conservative re-scaling correction
  !> Step-rescale calls handle inter-rank communication for marking critical planes
  subroutine RungeKutta_4th_rescale(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    integer(4), intent(in) :: id_RungeKutta
    integer(4), intent(in) :: id_rescale
    integer, intent(in)    :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)    :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)    :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(nx,5,ny,nz)
    integer i, j, k, l, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireq2(2)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    ! rescale !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: step = 1, flag_re = 0, flag_req
    real(8), allocatable, device :: Qre(:), Qm(:)
    real(8), allocatable, pinned :: Qm_cpu(:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), zetaz(:), Jacobian(:,:), dtdxdy(:,:), dtdydz(:,:), dtdzdx(:,:)
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
      call allocate_device_mem(myrank, nx, ny, nz, dtdxdy, dtdydz, dtdzdx, xix, etay, zetaz, Jacobian, ruvwp, T, mu, mut, qc2, E, F, G)
      allocate(QJ(nx,5,ny,nz), QJs(nx,5,ny,nz), Rs(nx-2,5,ny-2,nz-2))
      print *, "myrank is ", myrank, " memory allocation has completed"
      call pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, overlap, &
                    dtdxdy, dtdydz, dtdzdx, xix, etay, zetaz, Jacobian, QJ, ke0, entropy0, .false.)
      Rs = 0.d0
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
          call step_rescale(1, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJ, Qm, Qre)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, 1, nz-1)
          call calc_step<<<blocks,threads>>>(nx, ny, nz, 0.5d0, 1.d0, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ, QJs, Rs) ! QJs = Q2
          call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs, Qre)
        elseif (myrank == rerank+1) then
          call rescale_recv_send(1, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          call step_rescale(2, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJs, Qm, Qre)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, 1, nz-1)
          call calc_step<<<blocks,threads>>>(nx, ny, nz, 0.5d0, 2.d0, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ, QJs, Rs) ! QJs = Q3
          call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs, Qre)
        elseif (myrank == rerank+1) then
          call rescale_recv_send(2, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          call step_rescale(3, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJs, Qm, Qre)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, 1, nz-1)
          call calc_step<<<blocks,threads>>>(nx, ny, nz, 1.0d0, 2.d0, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ, QJs, Rs) ! QJs = Q4
          call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs, Qre)
        elseif (myrank == rerank+1) then
          call rescale_recv_send(3, flag_re, nx, ny, nz, np*(t2-1)+t1, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          call step_rescale(4, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, QJs, Qm, Qre)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, 1, nz-1)
          call calc_step4<<<blocks,threads>>>(nx, ny, nz, dtdxdy, dtdydz, dtdzdx, E, F, G, Rs, QJ)
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
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJs, Rs, E, F, G, xix, etay, zetaz, Jacobian, dtdxdy, dtdydz, dtdzdx)
    endif
    if (myrank == 0 .or. myrank == rerank) then
      deallocate(Qre, Qm)
    elseif (myrank == rerank+1) then
      deallocate(Qm_cpu)
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_4th_rescale

  !> 3rd-order TVD Runge-Kutta with z-direction domain decomposition and
  !> overlapped MPI communication.  All MPI ranks compute on GPU.  At each RK
  !> stage: post non-blocking z-halo exchange → compute interior fluxes (GPU
  !> runs while MPI is in flight) → complete exchange → compute halo fluxes.
  !> Activated by setting id_RungeKutta to integer(8) in mod_globals.f90.
  subroutine RungeKutta_3rd_zdec(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, &
                                   x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    integer(8), intent(in) :: id_RungeKutta
    integer(2), intent(in) :: id_rescale
    integer, intent(in)    :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1), y(ny), dy_cpu(ny-1)
    real(8), intent(in)    :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(nx,5,ny,nz)
    integer overlap, ierr, nranks, ndevices, stat, t1, t2
    integer req_z(4)
    integer istat(MPI_STATUS_SIZE)
    ! GPU arrays (compute ranks only)
    real(8), allocatable, device :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJ2(:,:,:,:)
    real(8), allocatable, device :: E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), zetaz(:), Jacobian(:,:)
    real(8), allocatable, device :: dtdxdy(:,:), dtdydz(:,:), dtdzdx(:,:)
    real(8), allocatable, device :: Qs1d_lo(:), Qs1d_hi(:), Qr1d_lo(:), Qr1d_hi(:)
    real(8), allocatable, pinned :: send_lo(:), send_hi(:), recv_lo(:), recv_hi(:)
    real(4) :: ke0 = 1.e0, entropy0 = 1.e0

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    if (mod(myrank,2) == 0) then
      stat = cudaGetDeviceCount(ndevices)
      if (myrank == 0) print '(2x, i2, a)', ndevices, " GPU devices are found"
      call check_gpu(mygpu)
      call allocate_device_mem(myrank, nx, ny, nz, dtdxdy, dtdydz, dtdzdx, xix, etay, zetaz, &
                               Jacobian, ruvwp, T, mu, mut, qc2, E, F, G)
      allocate(QJ(nx,5,ny,nz), QJ2(nx,5,ny,nz), stat=ierr)
      if (ierr /= 0) then
        print *, "myrank is ", myrank, " memory allocation failed", ierr
      else
        print *, "myrank is ", myrank, " memory allocation has completed"
      endif
      call pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, &
                    Jacobian_cpu, Q, overlap, dtdxdy, dtdydz, dtdzdx, xix, etay, zetaz, &
                    Jacobian, QJ, ke0, entropy0, .false.)
      allocate(Qs1d_lo(nx*ny*overlap*5), Qs1d_hi(nx*ny*overlap*5), &
               Qr1d_lo(nx*ny*overlap*5), Qr1d_hi(nx*ny*overlap*5), stat=ierr)
      allocate(send_lo(nx*ny*overlap*5), send_hi(nx*ny*overlap*5), &
               recv_lo(nx*ny*overlap*5), recv_hi(nx*ny*overlap*5), stat=ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      print *, "myrank is ", myrank, " start Runge-Kutta (z-decomp)"
      do t2 = 1, np
        do t1 = 1, nt
          call start_exchange_z(myrank, nranks, overlap, nx, ny, nz, QJ, &
                                     Qs1d_lo, Qs1d_hi, Qr1d_lo, Qr1d_hi, &
                                     send_lo, send_hi, recv_lo, recv_hi, req_z)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, &
                        ruvwp, T, mu, mut, qc2, E, F, G, 2*overlap, nz-2*overlap)
          call finish_exchange_z(overlap, nx, ny, nz, recv_lo, recv_hi, Qr1d_lo, Qr1d_hi, QJ, req_z)
          call calc_EFG_halo(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, &
                                ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step1<<<blocks,threads>>>(nx, ny, nz, 1.d0, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ, QJ2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ2)

          call start_exchange_z(myrank, nranks, overlap, nx, ny, nz, QJ2, &
                                     Qs1d_lo, Qs1d_hi, Qr1d_lo, Qr1d_hi, &
                                     send_lo, send_hi, recv_lo, recv_hi, req_z)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, &
                        ruvwp, T, mu, mut, qc2, E, F, G, 2*overlap, nz-2*overlap)
          call finish_exchange_z(overlap, nx, ny, nz, recv_lo, recv_hi, Qr1d_lo, Qr1d_hi, QJ2, req_z)
          call calc_EFG_halo(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, &
                                ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step2_3<<<blocks,threads>>>(nx, ny, nz, 0.75d0, 0.25d0, 0.25d0, 1.d0, &
                                                dtdxdy, dtdydz, dtdzdx, E, F, G, QJ, QJ2)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ2)

          call start_exchange_z(myrank, nranks, overlap, nx, ny, nz, QJ2, &
                                     Qs1d_lo, Qs1d_hi, Qr1d_lo, Qr1d_hi, &
                                     send_lo, send_hi, recv_lo, recv_hi, req_z)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, &
                        ruvwp, T, mu, mut, qc2, E, F, G, 2*overlap, nz-2*overlap)
          call finish_exchange_z(overlap, nx, ny, nz, recv_lo, recv_hi, Qr1d_lo, Qr1d_hi, QJ2, req_z)
          call calc_EFG_halo(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, &
                                ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step2_3<<<blocks,threads>>>(nx, ny, nz, 2.d0, 1.d0, 2.d0, one_third, &
                                                dtdxdy, dtdydz, dtdzdx, E, F, G, QJ2, QJ)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ)
        enddo
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, QJ, Q, ke0, entropy0)
      enddo
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJ2, E, F, G, xix, etay, zetaz, Jacobian, &
                 dtdxdy, dtdydz, dtdzdx)
      deallocate(Qs1d_lo, Qs1d_hi, Qr1d_lo, Qr1d_hi, send_lo, send_hi, recv_lo, recv_hi)
      print *, "myrank is ", myrank, " deallocate GPU memory"
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank, MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank, MPI_COMM_WORLD, istat, ierr)
      call send_recv_for_print_odd(myrank, nranks, 0, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      do t2 = 1, np
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
      enddo
    endif
  end subroutine RungeKutta_3rd_zdec

  !> 4th-order classical Runge-Kutta with z-direction domain decomposition and
  !> overlapped MPI communication.  All MPI ranks compute on GPU.  At each RK
  !> stage: post non-blocking z-halo exchange → compute interior fluxes (GPU
  !> runs while MPI is in flight) → complete exchange → compute halo fluxes.
  !> Activated by setting id_RungeKutta to integer(8) and id_rescale to integer(4)
  !> in mod_globals.f90.  The integer(4) rescale kind distinguishes this procedure
  !> from RungeKutta_3rd_zdec (which uses integer(2) for id_rescale) in the generic
  !> interface resolution.
  subroutine RungeKutta_4th_zdec(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, &
                                   x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    integer(8), intent(in) :: id_RungeKutta
    integer(4), intent(in) :: id_rescale
    integer, intent(in)    :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1), y(ny), dy_cpu(ny-1)
    real(8), intent(in)    :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(nx,5,ny,nz)
    integer overlap, ierr, nranks, ndevices, stat, t1, t2
    integer req_z(4)
    integer istat(MPI_STATUS_SIZE)
    ! GPU arrays (compute ranks only)
    real(8), allocatable, device :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:)
    real(8), allocatable, device :: E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), zetaz(:), Jacobian(:,:)
    real(8), allocatable, device :: dtdxdy(:,:), dtdydz(:,:), dtdzdx(:,:)
    real(8), allocatable, device :: Qs1d_lo(:), Qs1d_hi(:), Qr1d_lo(:), Qr1d_hi(:)
    real(8), allocatable, pinned :: send_lo(:), send_hi(:), recv_lo(:), recv_hi(:)
    real(4) :: ke0 = 1.e0, entropy0 = 1.e0

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    if (mod(myrank,2) == 0) then
      stat = cudaGetDeviceCount(ndevices)
      if (myrank == 0) print '(2x, i2, a)', ndevices, " GPU devices are found"
      call check_gpu(mygpu)
      call allocate_device_mem(myrank, nx, ny, nz, dtdxdy, dtdydz, dtdzdx, xix, etay, zetaz, &
                               Jacobian, ruvwp, T, mu, mut, qc2, E, F, G)
      allocate(QJ(nx,5,ny,nz), QJs(nx,5,ny,nz), Rs(nx-2,5,ny-2,nz-2), stat=ierr)
      if (ierr /= 0) then
        print *, "myrank is ", myrank, " memory allocation failed", ierr
      else
        print *, "myrank is ", myrank, " memory allocation has completed"
      endif
      call pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, &
                    Jacobian_cpu, Q, overlap, dtdxdy, dtdydz, dtdzdx, xix, etay, zetaz, &
                    Jacobian, QJ, ke0, entropy0, .true.)
      allocate(Qs1d_lo(nx*ny*overlap*5), Qs1d_hi(nx*ny*overlap*5), &
               Qr1d_lo(nx*ny*overlap*5), Qr1d_hi(nx*ny*overlap*5), stat=ierr)
      allocate(send_lo(nx*ny*overlap*5), send_hi(nx*ny*overlap*5), &
               recv_lo(nx*ny*overlap*5), recv_hi(nx*ny*overlap*5), stat=ierr)
      Rs = 0.d0
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      print *, "myrank is ", myrank, " start Runge-Kutta (z-decomp, RK4)"
      do t2 = 1, np
        do t1 = 1, nt
          call start_exchange_z(myrank, nranks, overlap, nx, ny, nz, QJ, &
                                     Qs1d_lo, Qs1d_hi, Qr1d_lo, Qr1d_hi, &
                                     send_lo, send_hi, recv_lo, recv_hi, req_z)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, &
                        ruvwp, T, mu, mut, qc2, E, F, G, 2*overlap, nz-2*overlap)
          call finish_exchange_z(overlap, nx, ny, nz, recv_lo, recv_hi, Qr1d_lo, Qr1d_hi, QJ, req_z)
          call calc_EFG_halo(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, &
                                ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step<<<blocks,threads>>>(nx, ny, nz, 0.5d0, 1.d0, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ, QJs, Rs)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)

          call start_exchange_z(myrank, nranks, overlap, nx, ny, nz, QJs, &
                                     Qs1d_lo, Qs1d_hi, Qr1d_lo, Qr1d_hi, &
                                     send_lo, send_hi, recv_lo, recv_hi, req_z)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, &
                        ruvwp, T, mu, mut, qc2, E, F, G, 2*overlap, nz-2*overlap)
          call finish_exchange_z(overlap, nx, ny, nz, recv_lo, recv_hi, Qr1d_lo, Qr1d_hi, QJs, req_z)
          call calc_EFG_halo(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, &
                                ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step<<<blocks,threads>>>(nx, ny, nz, 0.5d0, 2.d0, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ, QJs, Rs)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)

          call start_exchange_z(myrank, nranks, overlap, nx, ny, nz, QJs, &
                                     Qs1d_lo, Qs1d_hi, Qr1d_lo, Qr1d_hi, &
                                     send_lo, send_hi, recv_lo, recv_hi, req_z)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, &
                        ruvwp, T, mu, mut, qc2, E, F, G, 2*overlap, nz-2*overlap)
          call finish_exchange_z(overlap, nx, ny, nz, recv_lo, recv_hi, Qr1d_lo, Qr1d_hi, QJs, req_z)
          call calc_EFG_halo(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, &
                                ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step<<<blocks,threads>>>(nx, ny, nz, 1.0d0, 2.d0, dtdxdy, dtdydz, dtdzdx, E, F, G, QJ, QJs, Rs)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)

          call start_exchange_z(myrank, nranks, overlap, nx, ny, nz, QJs, &
                                     Qs1d_lo, Qs1d_hi, Qr1d_lo, Qr1d_hi, &
                                     send_lo, send_hi, recv_lo, recv_hi, req_z)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, &
                        ruvwp, T, mu, mut, qc2, E, F, G, 2*overlap, nz-2*overlap)
          call finish_exchange_z(overlap, nx, ny, nz, recv_lo, recv_hi, Qr1d_lo, Qr1d_hi, QJs, req_z)
          call calc_EFG_halo(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, &
                                ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step4<<<blocks,threads>>>(nx, ny, nz, dtdxdy, dtdydz, dtdzdx, E, F, G, Rs, QJ)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ)
          Rs = 0.d0
        enddo
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, QJ, Q, ke0, entropy0)
      enddo
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJs, Rs, E, F, G, xix, etay, zetaz, Jacobian, &
                 dtdxdy, dtdydz, dtdzdx)
      deallocate(Qs1d_lo, Qs1d_hi, Qr1d_lo, Qr1d_hi, send_lo, send_hi, recv_lo, recv_hi)
      print *, "myrank is ", myrank, " deallocate GPU memory"
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank, MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank, MPI_COMM_WORLD, istat, ierr)
      call send_recv_for_print_odd(myrank, nranks, 0, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      do t2 = 1, np
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
      enddo
    endif
  end subroutine RungeKutta_4th_zdec
end module calc_time_dev
