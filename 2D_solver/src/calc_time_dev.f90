module calc_time_dev
  use cudafor
  use curand
  use curand_device
  use mpi
  use nvtx
  use mod_globals, only : id_visc, id_LL, id_igr, id_force, nt, np, rerank, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, &
  & blocksEv, blocksFv, threadsEv, threadsFv
  use calc_flux_base
  use calc_steps
  use calc_rand
  use set
  use preprocess
  use print
  implicit none
  interface RungeKutta
    module procedure RungeKutta_3rd, RungeKutta_4th, Gauss_RungeKutta
  end interface
contains
  subroutine RungeKutta_3rd(id_RungeKutta, myrank, mygpu, nx, ny, nz, x_cpu, dx_cpu, y_cpu, dy_cpu, z_cpu, dz_cpu, Jacobian_cpu, Q)
    integer(kind=2), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)         :: x_cpu(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y_cpu(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z_cpu(nz), dz_cpu(1), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(4,nx,ny,nz)
    integer i, j, l, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireq2(2)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: ruvp(:,:,:), QJ(:,:,:), QJ2(:,:,:), E(:,:,:), F(:,:,:)
    real(8), allocatable, device :: T(:,:), mu(:,:)
    real(8), allocatable, device :: x(:), y(:), dx(:), dy(:), xix(:), etay(:), Jacobian(:,:)
    ! IGR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(4), allocatable, device :: sigma(:,:)
    real(4) ke0, entropy0
    ! LL & forcing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: Fout(:,:,:)
    type(curandGenerator) gen
    type(curandStateXORWOW), allocatable, device :: state(:)
    istat = curandCreateGenerator(gen, CURAND_RNG_PSEUDO_DEFAULT)
    istat = curandSetPseudoRandomGeneratorSeed(gen, 12345_8)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! count GPU
    stat = cudaGetDeviceCount(ndevices)
    print *, "rank", myrank, " has found ", ndevices, " GPU devices"
    if (mod(myrank,2) == 0) then
      call check_gpu(mygpu)
      call allocate_device_mem(myrank, nx, ny, x, y, dx, dy, xix, etay, Jacobian, ruvp, T, mu, E, F, Fout)
      allocate(QJ(4,nx,ny), QJ2(4,nx,ny), stat=ierr)
      if (ierr /= 0) then
        print *, "myrank is ", myrank, " memory allocation failed", ierr
      else
        print *, "myrank is ", myrank, " memory allocation has completed"
      endif
      call pre_calc(nx, ny, myrank, nranks, x_cpu, dx_cpu, y_cpu, dy_cpu, Jacobian_cpu, Q(:,:,:,1), overlap, &
                    x, y, dx, dy, xix, etay, Jacobian, QJ)
      if (kind(id_LL) == 4) then
        call init_state(nx, ny, state)
      endif
      if (kind(id_igr) == 4) then
        allocate(sigma(nx,ny))
        call set_init_sigma(nx, ny, dx_cpu, dy_cpu, Jacobian_cpu, Q(:,:,:,1), sigma)
      endif
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    print *, "myrank is ", myrank, " start Runge-Kutta"
    do t2 = 1, np
      do t1 = 1, nt
        if (mod(myrank,2) == 0) then
          call nvtxStartRange("calc flux", 1)
          if (kind(id_LL) == 4) then
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJ, ruvp, T, mu, E, F, Fout, gen, sigma, state)
          else
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJ, ruvp, T, mu, E, F, Fout, gen, sigma)
          endif
          !print *, "myrank is ", myrank, " calc EFG"
          call nvtxEndRange
          call nvtxStartRange("calc step", 2)
          call calc_step1(id_force, nx, ny, 1.d0, dx, dy, E, F, Fout, QJ, QJ2)
          !print *, "myrank is ", myrank, " calc step"
          call nvtxEndRange
          call set_bc(myrank, nx, ny, x, y, Jacobian, QJ2)
          !print *, "myrank is ", myrank, " set bc"
          call nvtxEndRange
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_LL) == 4) then
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJ2, ruvp, T, mu, E, F, Fout, gen, sigma, state)
          else
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJ2, ruvp, T, mu, E, F, Fout, gen, sigma)
          endif
          call calc_step2_3(id_force, nx, ny, 0.75d0, 0.25d0, 0.25d0, 1.d0, dx, dy, E, F, Fout, QJ, QJ2)
          call set_bc(myrank, nx, ny, x, y, Jacobian, QJ2)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_LL) == 4) then
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJ2, ruvp, T, mu, E, F, Fout, gen, sigma, state)
          else
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJ2, ruvp, T, mu, E, F, Fout, gen, sigma)
          endif
          call calc_step2_3(id_force, nx, ny, 2.d0, 1.d0, 2.d0, 3.d0, dx, dy, E, F, Fout, QJ2, QJ)
          call set_bc(myrank, nx, ny, x, y, Jacobian, QJ)
        endif
      enddo
      if (mod(myrank, 2) == 0) then
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x_cpu, y_cpu, z_cpu, Jacobian_cpu, QJ, Q(:,:,:,1), ke0, entropy0)
      else
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x_cpu, y_cpu, z_cpu, Jacobian_cpu, Q(:,:,:,1), ke0, entropy0)
      endif
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(ruvp, T, mu, QJ, QJ2, E, F, Fout, x, y, dx, dy, xix, etay, Jacobian)
      if (kind(id_LL) == 4) then
        deallocate(state)
      endif
      if (kind(id_igr) == 4) then
        deallocate(sigma)
      endif
    endif
    istat = curandDestroyGenerator(gen)
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_3rd


  subroutine RungeKutta_4th(id_RungeKutta, myrank, mygpu, nx, ny, nz, x_cpu, dx_cpu, y_cpu, dy_cpu, z_cpu, dz_cpu, Jacobian_cpu, Q)
    integer(kind=4), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)         :: x_cpu(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y_cpu(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z_cpu(nz), dz_cpu(1), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(4,nx,ny,nz)
    integer i, j, l, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireq2(2)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: ruvp(:,:,:), QJ(:,:,:), QJs(:,:,:), Rs(:,:,:), E(:,:,:), F(:,:,:)
    real(8), allocatable, device :: T(:,:), mu(:,:)
    real(8), allocatable, device :: x(:), y(:), dx(:), dy(:), xix(:), etay(:), Jacobian(:,:)
    ! IGR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(4), allocatable, device :: sigma(:,:)
    real(4) ke0, entropy0
    ! LL & forcing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: Fout(:,:,:)
    type(curandGenerator) gen
    type(curandStateXORWOW), allocatable, device :: state(:)
    istat = curandCreateGenerator(gen, CURAND_RNG_PSEUDO_DEFAULT)
    istat = curandSetPseudoRandomGeneratorSeed(gen, 12345_8)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! count GPU
    stat = cudaGetDeviceCount(ndevices)
    if (myrank == 0) then
      print '(2x, i2, a)', ndevices, " GPU devices are found"
    endif
    if (mod(myrank,2) == 0) then
      call check_gpu(mygpu)
      call allocate_device_mem(myrank, nx, ny, x, y, dx, dy, xix, etay, Jacobian, ruvp, T, mu, E, F, Fout)
      allocate(QJ(4,nx,ny), QJs(4,nx,ny), Rs(4,nx-2,ny-2))
      print *, "myrank is ", myrank, " memory allocation has completed"
      call pre_calc(nx, ny, myrank, nranks, x_cpu, dx_cpu, y_cpu, dy_cpu, Jacobian_cpu, Q(:,:,:,1), overlap, &
                    x, y, dx, dy, xix, etay, Jacobian, QJ)
      Rs = 0.d0
      if (kind(id_LL) == 4) then
        call init_state(nx, ny, state)
      endif
      if (kind(id_igr) == 4) then
        allocate(sigma(nx,ny))
        call set_init_sigma(nx, ny, dx_cpu, dy_cpu, Jacobian_cpu, Q(:,:,:,1), sigma)
      endif
    endif
    
    do t2 = 1, np
      do t1 = 1, nt
        if (mod(myrank,2) == 0) then
          if (kind(id_LL) == 4) then
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJ, ruvp, T, mu, E, F, Fout, gen, sigma, state)
          else
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJ, ruvp, T, mu, E, F, Fout, gen, sigma)
          endif
          call calc_step(id_force, nx, ny, 0.5d0, 1.d0, dx, dy, E, F, Fout, QJ, QJs, Rs) ! QJs = Q2
          call set_bc(myrank, nx, ny, x, y, Jacobian, QJs)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_LL) == 4) then
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJs, ruvp, T, mu, E, F, Fout, gen, sigma, state)
          else
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJs, ruvp, T, mu, E, F, Fout, gen, sigma)
          endif
          call calc_step(id_force, nx, ny, 0.5d0, 2.d0, dx, dy, E, F, Fout, QJ, QJs, Rs) ! QJs = Q3
          call set_bc(myrank, nx, ny, x, y, Jacobian, QJs)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_LL) == 4) then
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJs, ruvp, T, mu, E, F, Fout, gen, sigma, state)
          else
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJs, ruvp, T, mu, E, F, Fout, gen, sigma)
          endif
          call calc_step(id_force, nx, ny, 1.0d0, 2.d0, dx, dy, E, F, Fout, QJ, QJs, Rs) ! QJs = Q4
          call set_bc(myrank, nx, ny, x, y, Jacobian, QJs)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_LL) == 4) then
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJs, ruvp, T, mu, E, F, Fout, gen, sigma, state)
          else
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJs, ruvp, T, mu, E, F, Fout, gen, sigma)
          endif
          call calc_step4(id_force, nx, ny, dx, dy, E, F, Fout, Rs, QJ)
          call set_bc(myrank, nx, ny, x, y, Jacobian, QJ)
        endif
      enddo
      if (mod(myrank, 2) == 0) then
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x_cpu, y_cpu, z_cpu, Jacobian_cpu, QJ, Q(:,:,:,1), ke0, entropy0)
      else
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x_cpu, y_cpu, z_cpu, Jacobian_cpu, Q(:,:,:,1), ke0, entropy0)
      endif
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(ruvp, T, mu, QJ, QJs, Rs, E, F, Fout, x, y, dx, dy, xix, etay, Jacobian)
      if (kind(id_LL) == 4) then
        deallocate(state)
      endif
      if (kind(id_igr) == 4) then
        deallocate(sigma)
      endif
    endif
    istat = curandDestroyGenerator(gen)
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_4th
 

  subroutine Gauss_RungeKutta(id_RungeKutta, myrank, mygpu, nx, ny, nz, x_cpu, dx_cpu, y_cpu, dy_cpu, z_cpu, dz_cpu, Jacobian_cpu, Q)
    integer(kind=8), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)         :: x_cpu(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y_cpu(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z_cpu(nz), dz_cpu(1), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(4,nx,ny,nz)
    integer i, j, itr, max_itr, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireqs(2)
    integer istat(MPI_STATUS_SIZE), istats(MPI_STATUS_SIZE,2)
    integer :: step
    real(8) :: c1, c2, a11, a12, a21, a22, b1, b2, err, tol = 1.d-16
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: ruvp(:,:,:), QJ(:,:,:), QJs(:,:,:), E(:,:,:), F(:,:,:)
    real(8), allocatable, device :: T(:,:), mu(:,:) 
    real(8), allocatable, device :: R1(:,:,:), R2(:,:,:), R1_new(:,:,:), R2_new(:,:,:)
    real(8), allocatable, device :: x(:), y(:), dx(:), dy(:), xix(:), etay(:), Jacobian(:,:)
    ! IGR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(4), allocatable, device :: sigma(:,:)
    real(4) ke0, entropy0
    ! LL & forcing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: Fout(:,:,:)
    type(curandGenerator) gen
    type(curandStateXORWOW), allocatable, device :: state(:)
    istat = curandCreateGenerator(gen, CURAND_RNG_PSEUDO_DEFAULT)
    istat = curandSetPseudoRandomGeneratorSeed(gen, 12345_8)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! count GPU
    stat = cudaGetDeviceCount(ndevices)
    if (myrank == 0) then
      print '(2x, i2, a)', ndevices, " GPU devices are found"
    endif
    if (mod(myrank,2) == 0) then
      call check_gpu(mygpu)
      call allocate_device_mem(myrank, nx, ny, x, y, dx, dy, xix, etay, Jacobian, ruvp, T, mu, E, F, Fout)
      allocate(QJ(4,nx,ny), QJs(4,nx,ny), R1(4,nx-2,ny-2), R2(4,nx-2,ny-2))
      allocate(R1_new(4,nx-2,ny-2), R2_new(4,nx-2,ny-2))
      print *, "myrank is ", myrank, " memory allocation has completed"
      call pre_calc(nx, ny, myrank, nranks, x_cpu, dx_cpu, y_cpu, dy_cpu, Jacobian_cpu, Q(:,:,:,1), overlap, &
                    x, y, dx, dy, xix, etay, Jacobian, QJ)
      if (kind(id_LL) == 4) then
        call init_state(nx, ny, state)
      endif
      if (kind(id_igr) == 4) then
        allocate(sigma(nx,ny))
        call set_init_sigma(nx, ny, dx_cpu, dy_cpu, Jacobian_cpu, Q(:,:,:,1), sigma)
      endif
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
          call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJ, ruvp, T, mu, E, F, Fout, gen, sigma)
          call calc_step1(id_force, nx, ny, c1, dx, dy, E, F, Fout, QJ, QJs)
          call set_bc(myrank, nx, ny, x, y, Jacobian, QJs)
          call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJs, ruvp, T, mu, E, F, Fout, gen, sigma)
          call calc_R(id_force, nx, ny, dx, dy, E, F, Fout, R1)
          ! calc R2
          call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJ, ruvp, T, mu, E, F, Fout, gen, sigma)
          call calc_step1(id_force, nx, ny, c2, dx, dy, E, F, Fout, QJ, QJs)
          call set_bc(myrank, nx, ny, x, y, Jacobian, QJs)
          call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJs, ruvp, T, mu, E, F, Fout, gen, sigma)
          call calc_R(id_force, nx, ny, dx, dy, E, F, Fout, R2)
          do itr = 1, max_itr
            ! calc R1
            call calc_Gauss_step(nx, ny, a11, a12, dx, dy, R1, R2, QJ, QJs)
            call set_bc(myrank, nx, ny, x, y, Jacobian, QJs)
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJs, ruvp, T, mu, E, F, Fout, gen, sigma)
            call calc_R(id_force, nx, ny, dx, dy, E, F, Fout, R1_new)
            ! calc R2
            call calc_Gauss_step(nx, ny, a21, a22, dx, dy, R1, R2, QJ, QJs)
            call set_bc(myrank, nx, ny, x, y, Jacobian, QJs)
            call calc_EF(id_visc, id_force, nx, ny, x, y, xix, etay, Jacobian, QJs, ruvp, T, mu, E, F, Fout, gen, sigma)
            call calc_R(id_force, nx, ny, dx, dy, E, F, Fout, R2_new)
            call calc_error(nx, ny, R1, R2, R1_new, R2_new, err)
            if (err < tol) exit
            R1 = R1_new
            R2 = R2_new
          enddo
          if (err < tol) then
            print *, "Converged at itr=", itr
          else
            print *, "Didn't Converged error=", err
          endif
          call calc_Gauss_step_Q(nx, ny, b1, b2, xix, etay, R1, R2, QJ)
          call set_bc(myrank, nx, ny, x, y, Jacobian, QJ)
        endif
      enddo
      if (mod(myrank, 2) == 0) then
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x_cpu, y_cpu, z_cpu, Jacobian_cpu, QJ, Q(:,:,:,1), ke0, entropy0)
      else
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x_cpu, y_cpu, z_cpu, Jacobian_cpu, Q(:,:,:,1), ke0, entropy0)
      endif
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(ruvp, T, mu, QJ, QJs, R1, R2, R1_new, R2_new, E, F, Fout, x, y, dx, dy, xix, etay, Jacobian)
      if (kind(id_LL) == 4) then
        deallocate(state)
      endif
      if (kind(id_igr) == 4) then
        deallocate(sigma)
      endif
    endif
    istat = curandDestroyGenerator(gen)
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine Gauss_RungeKutta
end module calc_time_dev

