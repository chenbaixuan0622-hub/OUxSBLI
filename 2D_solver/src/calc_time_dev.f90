module calc_time_dev
  use cudafor
  use mpi
  use nvtx
  use mod_globals, only : accuracy, id_accuracy, id_forcing, id_exchange, nt, np
  use calc_flux_base
  use calc_steps
  use set
  use print
  implicit none
  interface RungeKutta
    module procedure RungeKutta_3rd, RungeKutta_4th, Gauss_RungeKutta
  end interface
contains
  subroutine check_gpu(mygpu)
    integer, intent(in) :: mygpu
    integer ilen, stat
    type(cudaDeviceProp) prop
    stat = cudaSetDevice(mygpu)
    stat = cudaGetDeviceProperties(prop, mygpu)
    ilen = verify(prop%name, ' ', .true.)
    print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", mygpu, ") is available"
  end subroutine check_gpu


  subroutine pre_calc(nx, ny, myrank, nranks, x, dx_cpu, y, dy_cpu, Jacobian_cpu, Q, overlap, xix, etay, Jacobian, QJ)
    integer, intent(in)    :: nx, ny, myrank, nranks
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1), y(ny), dy_cpu(ny-1), Jacobian_cpu(ny)
    real(8), intent(inout) :: Q(4,nx,ny,1)
    integer, intent(out)   :: overlap
    real(8), intent(out), device :: xix(nx-1), etay(ny-1), Jacobian(ny)
    real(8), intent(out), device :: QJ(4,nx,ny)
    real(8) xix_cpu(nx-1), etay_cpu(ny-1)
    real(4) rho1d(nx*ny), p1d(nx*ny), v1d(nx*ny*3)
    integer i, j, k, ierr
    ! set Q / Jacobian
    do j = 1, ny
      do i = 1, nx
        do k = 1, 4
          Q(k,i,j,1) = Q(k,i,j,1) / Jacobian_cpu(j)
    enddo;enddo;enddo
    ! copy on GPU
    xix_cpu  = 1.d0 / dx_cpu
    etay_cpu = 1.d0 / dy_cpu
    xix      = xix_cpu
    etay     = etay_cpu
    Jacobian = Jacobian_cpu
    QJ = Q(:,:,:,1)
    ! for multi GPU
    if (kind(id_accuracy) == 8) then
      overlap = 3
    elseif (kind(id_accuracy) == 4) then
      overlap = 2
    else
      overlap = 1
    endif
    call make_1d_for_print(nx, ny, Jacobian_cpu, Q(:,:,:,1), rho1d, p1d, v1d)
    call print_vtk(0, nx, ny, x, y, rho1d, p1d, v1d)
  end subroutine pre_calc


  subroutine RungeKutta_3rd(id_RungeKutta, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    use mod_globals, only : id_visc
    integer(kind=2), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z(nz), dz_cpu(nz-1), Jacobian_cpu(ny)
    real(8), intent(inout)      :: Q(4,nx,ny,nz)
    integer i, j, k, l, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireq2(2)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: QJ(:,:,:), QJ2(:,:,:), E(:,:,:), F(:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), Jacobian(:)
    ! forcing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: fx(:,:), fy(:,:)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! count GPU
    stat = cudaGetDeviceCount(ndevices)
    print *, "rank", myrank, " has found ", ndevices, " GPU devices"
    if (mod(myrank,2) == 0) then
      call check_gpu(mygpu)
      allocate(QJ(4,nx,ny), QJ2(4,nx,ny), E(4,nx-1,ny-2), F(4,nx-2,ny-1))
      allocate(xix(nx-1), etay(ny-1), Jacobian(ny))
      print *, "myrank is ", myrank, " memory allocation has completed"
      call pre_calc(nx, ny, myrank, nranks, x, dx_cpu, y, dy_cpu, Jacobian_cpu, Q, overlap, xix, etay, Jacobian, QJ)
      if (kind(id_forcing) == 4) then
        allocate(fx(nx-2,ny-2), fy(nx-2,ny-2))
        fx(:,:) = 0.d0
        fy(:,:) = 0.d0
      endif
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    print *, "myrank is ", myrank, " start Runge-Kutta"
    do t2 = 1, np
      do t1 = 1, nt
        if (mod(myrank,2) == 0) then
          if (kind(id_forcing) == 2) then
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJ, E, F)
            call calc_step(nx, ny, 1.d0, 0.d0, xix, etay, E, F, QJ, QJ2)
          else
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJ, E, F, fx, fy)
            call calc_step_forcing(nx, ny, 1.d0, 0.d0, xix, etay, E, F, fx, fy, QJ, QJ2)
          endif
          call set_bc(myrank, nx, ny, Jacobian, QJ2)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_forcing) == 2) then
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJ2, E, F)
            call calc_step2_3(nx, ny, 0.75d0, 0.25d0, 0.25d0, 1.d0, xix, etay, E, F, QJ, QJ2)
          else
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJ2, E, F, fx, fy)
            call calc_step2_3_forcing(nx, ny, 0.75d0, 0.25d0, 0.25d0, 1.d0, xix, etay, E, F, fx, fy, QJ, QJ2)
          endif
          call set_bc(myrank, nx, ny, Jacobian, QJ2)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_forcing) == 2) then
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJ2, E, F)
            call calc_step2_3(nx, ny, 2.d0, 1.d0, 2.d0, 3.d0, xix, etay, E, F, QJ2, QJ)
          else
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJ2, E, F, fx, fy)
            call calc_step2_3_forcing(nx, ny, 2.d0, 1.d0, 2.d0, 3.d0, xix, etay, E, F, fx, fy, QJ2, QJ)
          endif
          call set_bc(myrank, nx, ny, Jacobian, QJ)
        endif
      enddo
      call send_recv_for_print(myrank, nranks, t2, nx, ny, x, y, Jacobian_cpu, QJ, Q(:,:,:,1))
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(QJ, QJ2, E, F, xix, etay, Jacobian)
      if (kind(id_forcing) == 4) then
        deallocate(fx, fy)
      endif
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_3rd


  subroutine RungeKutta_4th(id_RungeKutta, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    use mod_globals, only : id_visc
    integer(kind=4), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z(nz), dz_cpu(nz-1), Jacobian_cpu(ny)
    real(8), intent(inout)      :: Q(4,nx,ny,nz)
    integer i, j, k, l, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireq2(2)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: QJ(:,:,:), QJs(:,:,:), Rs(:,:,:), E(:,:,:), F(:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), Jacobian(:)
    ! forcing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: fx(:,:), fy(:,:)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! count GPU
    stat = cudaGetDeviceCount(ndevices)
    if (myrank == 0) then
      print '(2x, i2, a)', ndevices, " GPU devices are found"
    endif
    if (mod(myrank,2) == 0) then
      call check_gpu(mygpu)
      allocate(QJ(4,nx,ny), QJs(4,nx,ny), Rs(4,nx-2,ny-2), E(4,nx-1,ny-2), F(4,nx-2,ny-1))
      allocate(xix(nx-1), etay(ny-1), Jacobian(ny))
      print *, "myrank is ", myrank, " memory allocation has completed"
      call pre_calc(nx, ny, myrank, nranks, x, dx_cpu, y, dy_cpu, Jacobian_cpu, Q, overlap, xix, etay, Jacobian, QJ)
      Rs = 0.d0
      if (kind(id_forcing) == 4) then
        allocate(fx(nx-2,ny-2), fy(nx-2,ny-2))
        fx(:,:) = 0.d0
        fy(:,:) = 0.d0
      endif
    endif
    
    do t2 = 1, np
      do t1 = 1, nt
        if (mod(myrank,2) == 0) then
          if (kind(id_forcing) == 2) then
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJ, E, F)
            call calc_step(nx, ny, 0.5d0, 1.d0, xix, etay, E, F, QJ, QJs, Rs) ! QJs = Q2
          else
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJ, E, F, fx, fy)
            call calc_step_forcing(nx, ny, 0.5d0, 1.d0, xix, etay, E, F, fx, fy, QJ, QJs, Rs) ! QJs = Q2
          endif
          call set_bc(myrank, nx, ny, Jacobian, QJs)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_forcing) == 2) then
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJs, E, F)
            call calc_step(nx, ny,  0.5d0, 2.d0, xix, etay, E, F, QJ, QJs, Rs) ! QJs = Q3
          else
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJs, E, F, fx, fy)
            call calc_step_forcing(nx, ny, 0.5d0, 2.d0, xix, etay, E, F, fx, fy, QJ, QJs, Rs) ! QJs = Q3
          endif
          call set_bc(myrank, nx, ny, Jacobian, QJs)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_forcing) == 2) then
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJs, E, F)
            call calc_step(nx, ny, 1.0d0, 2.d0, xix, etay, E, F, QJ, QJs, Rs) ! QJs = Q4
          else
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJs, E, F, fx, fy)
            call calc_step_forcing(nx, ny, 1.0d0, 2.d0, xix, etay, E, F, fx, fy, QJ, QJs, Rs) ! QJs = Q4
          endif
          call set_bc(myrank, nx, ny, Jacobian, QJs)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_forcing) == 2) then
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJs, E, F)
            call calc_step4(nx, ny, xix, etay, E, F, Rs, QJ)
          else
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJs, E, F, fx, fy)
            call calc_step4_forcing(nx, ny, xix, etay, E, F, fx, fy, Rs, QJ)
          endif
          call set_bc(myrank, nx, ny, Jacobian, QJ)
        endif
      enddo
      call send_recv_for_print(myrank, nranks, t2, nx, ny, x, y, Jacobian_cpu, QJ, Q(:,:,:,1))
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(QJ, QJs, Rs, E, F, xix, etay, Jacobian)
      if (kind(id_forcing) == 4) then
        deallocate(fx, fy)
      endif
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_4th
  
  subroutine Gauss_RungeKutta(id_RungeKutta, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    use mod_globals, only : id_visc
    integer(kind=8), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z(nz), dz_cpu(nz-1), Jacobian_cpu(ny)
    real(8), intent(inout)      :: Q(4,nx,ny,nz)
    integer i, j, k, itr, max_itr, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireqs(2)
    integer istat(MPI_STATUS_SIZE), istats(MPI_STATUS_SIZE,2)
    real(8) :: c1, c2, a11, a12, a21, a22, b1, b2, err, tol = 1.d-16
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: QJ(:,:,:), QJs(:,:,:), E(:,:,:), F(:,:,:)
    real(8), allocatable, device :: R1(:,:,:), R2(:,:,:), R1_new(:,:,:), R2_new(:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), Jacobian(:)
    ! forcing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: fx(:,:), fy(:,:)

    stat = cudaGetDeviceCount(ndevices)
    if (myrank == 0) then
      print '(2x, i2, a)', ndevices, " GPU devices are found"
    endif
    if (mod(myrank,2) == 0) then
      call check_gpu(mygpu)
      allocate(QJ(4,nx,ny), QJs(4,nx,ny), R1(4,nx-2,ny-2), R2(4,nx-2,ny-2))
      allocate(R1_new(4,nx-2,ny-2), R2_new(4,nx-2,ny-2))
      allocate(E(4,nx-1,ny-2), F(4,nx-2,ny-1))
      allocate(xix(nx-1), etay(ny-1), Jacobian(ny))
      print *, "myrank is ", myrank, " memory allocation has completed"
      call pre_calc(nx, ny, myrank, nranks, x, dx_cpu, y, dy_cpu, Jacobian_cpu, Q, overlap, xix, etay, Jacobian, QJ)
      if (kind(id_forcing) == 4) then
        allocate(fx(nx-2,ny-2), fy(nx-2,ny-2))
        fx(:,:) = 0.d0
        fy(:,:) = 0.d0
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
          call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJ, E, F)
          call calc_step(nx, ny, c1, 0.d0, xix, etay, E, F, QJ, QJs)
          call set_bc(myrank, nx, ny, Jacobian, QJs)
          call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJs, E, F)
          call calc_R(nx, ny, xix, etay, E, F, R1)
          ! calc R2
          call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJ, E, F)
          call calc_step(nx, ny, c2, 0.d0, xix, etay, E, F, QJ, QJs)
          call set_bc(myrank, nx, ny, Jacobian, QJs)
          call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJs, E, F)
          call calc_R(nx, ny, xix, etay, E, F, R2)
          do itr = 1, max_itr
            ! calc R1
            call calc_Gauss_step(nx, ny, a11, a12, xix, etay, R1, R2, QJ, QJs)
            call set_bc(myrank, nx, ny, Jacobian, QJs)
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJs, E, F)
            call calc_R(nx, ny, xix, etay, E, F, R1_new)
            ! calc R2
            call calc_Gauss_step(nx, ny, a21, a22, xix, etay, R1, R2, QJ, QJs)
            call set_bc(myrank, nx, ny, Jacobian, QJs)
            call calc_EFG(id_visc, nx, ny, xix, etay, Jacobian, QJs, E, F)
            call calc_R(nx, ny, xix, etay, E, F, R2_new)
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
          call set_bc(myrank, nx, ny, Jacobian, QJ)
        endif
      enddo
      call send_recv_for_print(myrank, nranks, t2, nx, ny, x, y, Jacobian_cpu, QJ, Q(:,:,:,1))
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(QJ, QJs, R1, R2, R1_new, R2_new, E, F, xix, etay, Jacobian)
      if (kind(id_forcing) == 4) then
        deallocate(fx, fy)
      endif
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine Gauss_RungeKutta
end module calc_time_dev

