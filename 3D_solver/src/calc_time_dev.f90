module calc_time_dev
  use cudafor
  use mpi
  use nvtx
  use mod_globals, only : accuracy, id_accuracy, id_scheme, id_turbulence, id_rescale, id_exchange, nt, np, nre2, rerank, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use calc_physical_quantities
  use calc_steps
  use calc_hybrid
  use calc_flux
  use calc_visc
  use calc_les
  use calc_rescale
  use calc_para
  use set
  use print
  implicit none
  interface calc_EFG
    module procedure calc_EFG_Euler, calc_EFG_visc, calc_EFG_LES
  end interface calc_EFG

  interface RungeKutta
    module procedure RungeKutta_3rd, RungeKutta_4th
  end interface
contains
  subroutine calc_EFG_Euler(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G)
    integer(kind=2), intent(in), value :: id_visc
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(in), device        :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device        :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device        :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device        :: Jacobian(ny)
    real(8), intent(in), device        :: QJ(nx,ny,nz,5) ! Q / Jacobian
    real(8), intent(out), device       :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device       :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device       :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p, sensor
    integer stat
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, rho, u, v, w, p)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, u, v, w, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(nx, ny, nz, rho, u, v, w, p, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(nx, ny, nz, rho, u, v, w, p, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(nx, ny, nz, rho, u, v, w, p, sensor, G)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_Euler


  subroutine calc_EFG_visc(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G)
    integer(kind=4), intent(in), value :: id_visc
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(in), device        :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device        :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device        :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device        :: Jacobian(ny)
    real(8), intent(in), device        :: QJ(nx,ny,nz,5) ! Q / Jacobian
    real(8), intent(out), device       :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device       :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device       :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p, sensor
    integer stat
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, rho, u, v, w, p)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, u, v, w, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(nx, ny, nz, rho, u, v, w, p, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(nx, ny, nz, rho, u, v, w, p, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(nx, ny, nz, rho, u, v, w, p, sensor, G)

    stat = cudaDeviceSynchronize()
    call calc_Ev<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, rho, u, v, w, p, E)
    call calc_Fv<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, rho, u, v, w, p, F)
    call calc_Gv<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, rho, u, v, w, p, G)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_visc
 

  subroutine calc_EFG_LES(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, E, F, G)
    integer(kind=8), intent(in), value :: id_visc
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(in), device        :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device        :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device        :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device        :: Jacobian(ny)
    real(8), intent(in), device        :: QJ(nx,ny,nz,5) ! Q / Jacobian
    real(8), intent(out), device       :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device       :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device       :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p, sensor, mut, qc2
    integer stat
    mut = 0.d0
    qc2 = 0.d0
    
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, rho, u, v, w, p)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, u, v, w, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(nx, ny, nz, rho, u, v, w, p, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(nx, ny, nz, rho, u, v, w, p, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(nx, ny, nz, rho, u, v, w, p, sensor, G)
    call calc_mut<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, rho, u, v, w, mut, qc2)
    stat = cudaDeviceSynchronize()
    call set_bc_mut(nx,ny,nz,mut,qc2)

    call calc_Ev_LES<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, rho, u, v, w, p, mut, qc2, E)
    call calc_Fv_LES<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, rho, u, v, w, p, mut, qc2, F)
    call calc_Gv_LES<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, rho, u, v, w, p, mut, qc2, G)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_LES


  subroutine RungeKutta_3rd(id_RungeKutta, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    use mod_globals, only : id_visc
    integer(kind=2), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z(nz), dz_cpu(nz-1), Jacobian_cpu(ny)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    integer i, j, k, l, t1, t2, overlap, itr, ierr, ilen, nranks, ndevices, stat, ireq, ireq2(2), ireq3(3)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2), istat3(MPI_STATUS_SIZE,3)
    real(8) xix_cpu(nx-1), etay_cpu(ny-1), zetaz_cpu(nz-1)
    ! rescale !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: step, flag_re = 0
    real(8), allocatable, device :: Qre(:), Qm(:)
    real(8), allocatable, pinned :: Qm_cpu(:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)         :: prop
    real(8), allocatable, device :: QJ(:,:,:,:), QJ2(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), zetaz(:), Jacobian(:)
    ! for plot
    real(4), allocatable :: rho1d(:), p1d(:), v1d(:)
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! check GPU
    stat = cudaGetDeviceCount(ndevices)
    print *, "rank", myrank, " has found ", ndevices, " GPU devices"
    ! plot
    allocate(rho1d(nx*ny*nz), p1d(nx*ny*nz), v1d(nx*ny*nz*3))
    if (mod(myrank,2) == 0) then
      stat = cudaSetDevice(mygpu)
      stat = cudaGetDeviceProperties(prop, mygpu)
      ilen = verify(prop%name, ' ', .true.)
      print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", mygpu, ") is available"
      allocate(QJ(nx,ny,nz,5), QJ2(nx,ny,nz,5), E(nx-1,ny-2,nz-2,5), F(nx-2,ny-1,nz-2,5), G(nx-2,ny-2,nz-1,5))
      allocate(xix(nx-1), etay(ny-1), zetaz(nz-1), Jacobian(ny))
      print *, "myrank is ", myrank, " memory allocation has completed"
      ! set Q / Jacobian
      do l = 1, 5
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              Q(i,j,k,l) = Q(i,j,k,l) / Jacobian_cpu(j)
      enddo;enddo;enddo;enddo
      ! copy on GPU
      xix_cpu   = 1.d0 / dx_cpu
      etay_cpu  = 1.d0 / dy_cpu
      zetaz_cpu = 1.d0 / dz_cpu
      xix      = xix_cpu
      etay     = etay_cpu
      zetaz    = zetaz_cpu
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
      call make_1d_for_print(nx, ny, nz, Jacobian_cpu, Q, rho1d, p1d, v1d)
      call print_vtk(0, nx, ny, nz, myrank+1, nranks, x, y, z, rho1d, p1d, v1d, ke0, entropy0)
      call MPI_SEND(ke0,      1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
      call MPI_SEND(entropy0, 1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
    endif
    ! rescale
    if (kind(id_rescale) == 4) then
      if (myrank == 0 .or. myrank == rerank) then
        allocate(Qre(ny*(nz-6)*5), Qm(ny*5))
      elseif (myrank == rerank+1) then
        stat = cudaSetDevice(0)
        stat = cudaGetDeviceProperties(prop, 0)
        ilen = verify(prop%name, ' ', .true.)
        print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", 0, ") calculates rescaling"
        allocate(Qm_cpu(ny*5))
      endif
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    print *, "myrank is ", myrank, " start Runge-Kutta"
    do t2 = 1, np
      do t1 = 1, nt
        step = np * (t2-1) + t1
        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJ, Qm, Qre)
          endif
          call nvtxStartRange("calc flux", 1)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, E, F, G)
          call nvtxEndRange
          print *, "myrank=", myrank, "calc flux"
          call nvtxStartRange("calc step", 2)
          call calc_step(nx, ny, nz, 1.d0, 0.d0, xix, etay, zetaz, E, F, G, QJ, QJ2)
          call nvtxEndRange
          print *, "myrank=", myrank, "calc step"
          if (ndevices >= 2 .and. kind(id_exchange) == 4) then
            call nvtxStartRange("exchange", 3)
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ2)
            call nvtxEndRange
            print *, "myrank=", myrank, "exchange"
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          endif
          call nvtxStartRange("set bc", 4)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ2, Qre)
          call nvtxEndRange
          print *, "myrank=", myrank, "set bc"
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call nvtxStartRange("calc rescale", 5)
          call rescale_recv_send(1, flag_re, nx, ny, nz, step, y, Jacobian_cpu, Qm_cpu)
          call nvtxEndRange
          print *, "myrank=", myrank, "calc rescale"
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJ2, Qm, Qre)
          endif
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, E, F, G)
          call calc_step2(nx, ny, nz, 0.75d0, 0.25d0, 0.25d0, 1.d0, xix, etay, zetaz, E, F, G, QJ, QJ2)
          if (ndevices >= 2 .and. kind(id_exchange) == 4) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ2)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          endif
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ2, Qre)
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(2, flag_re, nx, ny, nz, step, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJ2, Qm, Qre)
          endif
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, E, F, G)
          call calc_step3(nx, ny, nz, xix, etay, zetaz, E, F, G, QJ2, QJ)
          if (ndevices >= 2 .and. kind(id_exchange) == 4) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          endif
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ, Qre)
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(3, flag_re, nx, ny, nz, step, y, Jacobian_cpu, Qm_cpu)
        endif
      enddo

      ! send and recv device arrays
      if (mod(myrank,2) == 0) then
        Q = QJ
        call make_1d_for_print(nx, ny, nz, Jacobian_cpu, Q, rho1d, p1d, v1d)
        call MPI_ISEND(rho1d, nx*ny*nz,   MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ireq3(1), ierr) 
        call MPI_ISEND(p1d,   nx*ny*nz,   MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ireq3(2), ierr) 
        call MPI_ISEND(v1d,   nx*ny*nz*3, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ireq3(3), ierr) 
        call MPI_WAITALL(3, ireq3, istat3, ierr)
      else
        call MPI_IRECV(rho1d, nx*ny*nz,   MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, ireq3(1), ierr)
        call MPI_IRECV(p1d,   nx*ny*nz,   MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, ireq3(2), ierr)
        call MPI_IRECV(v1d,   nx*ny*nz*3, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, ireq3(3), ierr)
        call MPI_WAITALL(3, ireq3, istat3, ierr)
        call print_vtk(t2, nx, ny, nz, myrank, nranks, x, y, z, rho1d, p1d, v1d, ke0, entropy0)
      endif
    enddo

    deallocate(rho1d, p1d, v1d)
    if (mod(myrank,2) == 0) then
      deallocate(QJ, QJ2, E, F, G, xix, etay, zetaz, Jacobian)
    endif
    if (kind(id_rescale) == 4) then
      if (myrank == 0 .or. myrank == rerank) then
        deallocate(Qre, Qm)
      elseif (myrank == rerank+1) then
        deallocate(Qm_cpu)
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
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    integer i, j, k, l, t1, t2, itr, overlap, ierr, ilen, nranks, ndevices, stat, ireq, ireq2(2), ireq3(3)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2), istat3(MPI_STATUS_SIZE,3)
    real(8) xix_cpu(nx-1), etay_cpu(ny-1), zetaz_cpu(nz-1)
    ! rescale !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: step, flag_re = 0
    real(8), allocatable, device :: Qre(:), Qm(:)
    real(8), allocatable, pinned :: Qm_cpu(:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)         :: prop
    real(8), allocatable, device :: QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), zetaz(:), Jacobian(:)
    ! for plot
    real(4), allocatable :: rho1d(:), p1d(:), v1d(:)
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! check GPU
    stat = cudaGetDeviceCount(ndevices)
    if (myrank == 0) then
      print '(2x, i2, a)', ndevices, " GPU devices are found"
    endif
    ! plot
    allocate(rho1d(nx*ny*nz), p1d(nx*ny*nz), v1d(nx*ny*nz*3))
    if (mod(myrank,2) == 0) then
      stat = cudaSetDevice(mygpu)
      stat = cudaGetDeviceProperties(prop, mygpu)
      ilen = verify(prop%name, ' ', .true.)
      print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", mygpu, ") is available"
      allocate(QJ(nx,ny,nz,5), QJs(nx,ny,nz,5), Rs(nx-2,ny-2,nz-2,5), E(nx-1,ny-2,nz-2,5), F(nx-2,ny-1,nz-2,5), G(nx-2,ny-2,nz-1,5))
      allocate(xix(nx-1), etay(ny-1), zetaz(nz-1), Jacobian(ny))
      ! set Q / Jacobian
      do l = 1, 5
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              Q(i,j,k,l) = Q(i,j,k,l) / Jacobian_cpu(j)
      enddo;enddo;enddo;enddo
      ! copy on GPU
      xix_cpu   = 1.d0 / dx_cpu
      etay_cpu  = 1.d0 / dy_cpu
      zetaz_cpu = 1.d0 / dz_cpu
      xix       = xix_cpu
      etay      = etay_cpu
      zetaz     = zetaz_cpu
      Jacobian  = Jacobian_cpu
      QJ        = Q
      Rs        = 0.d0
      ! for multi GPU
      if (kind(id_accuracy) == 8) then
        overlap = 3
      elseif (kind(id_accuracy) == 4) then
        overlap = 2
      else
        overlap = 1
      endif
      call make_1d_for_print(nx, ny, nz, Jacobian_cpu, Q, rho1d, p1d, v1d)
      call print_vtk(0, nx, ny, nz, myrank+1, nranks, x, y, z, rho1d, p1d, v1d, ke0, entropy0)
      call MPI_SEND(ke0,      1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
      call MPI_SEND(entropy0, 1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
    endif
    ! rescale
    if (kind(id_rescale) == 4) then
      if (myrank == 0 .or. myrank == rerank) then
        allocate(Qre(ny*(nz-6)*5), Qm(ny*5))
      elseif (myrank == rerank+1) then
        stat = cudaSetDevice(0)
        stat = cudaGetDeviceProperties(prop, 0)
        ilen = verify(prop%name, ' ', .true.)
        print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", 0, ") calculates rescaling"
        allocate(Qm_cpu(ny*5))
      endif
    endif
    
    do t2 = 1, np
      do t1 = 1, nt
        step = np * (t2-1) + t1
        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJ, Qm, Qre)
          endif
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, E, F, G)
          call calc_step(nx, ny, nz, 0.5d0, 1.d0, xix, etay, zetaz, E, F, G, QJ, QJs, Rs) ! QJs = Q2
          if (ndevices >= 2) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJs)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          endif
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs, Qre)
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(1, flag_re, nx, ny, nz, step, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJs, Qm, Qre)
          endif
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, E, F, G)
          call calc_step(nx, ny, nz, 0.5d0, 2.d0, xix, etay, zetaz, E, F, G, QJ, QJs, Rs) ! QJs = Q3
          if (ndevices >= 2) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJs)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          endif
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs, Qre)
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(2, flag_re, nx, ny, nz, step, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJs, Qm, Qre)
          endif
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, E, F, G)
          call calc_step(nx, ny, nz, 1.0d0, 2.d0, xix, etay, zetaz, E, F, G, QJ, QJs, Rs) ! QJs = Q4
          if (ndevices >= 2) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJs)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          endif
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs, Qre)
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(3, flag_re, nx, ny, nz, step, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJs, Qm, Qre)
          endif
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, E, F, G)
          call calc_step4(nx, ny, nz, xix, etay, zetaz, E, F, G, Rs, QJ)
          if (ndevices >= 2) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
          endif
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ, Qre)
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(4, flag_re, nx, ny, nz, step, y, Jacobian_cpu, Qm_cpu)
        endif
      enddo

      ! send and recv device arrays
      if (mod(myrank,2) == 0) then
        Q = QJ
        call make_1d_for_print(nx, ny, nz, Jacobian_cpu, Q, rho1d, p1d, v1d)
        call MPI_ISEND(rho1d, nx*ny*nz,   MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ireq3(1), ierr) 
        call MPI_ISEND(p1d,   nx*ny*nz,   MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ireq3(2), ierr) 
        call MPI_ISEND(v1d,   nx*ny*nz*3, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ireq3(3), ierr) 
        call MPI_WAITALL(3, ireq3, istat3, ierr)
      else
        call MPI_IRECV(rho1d, nx*ny*nz,   MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, ireq3(1), ierr)
        call MPI_IRECV(p1d,   nx*ny*nz,   MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, ireq3(2), ierr)
        call MPI_IRECV(v1d,   nx*ny*nz*3, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, ireq3(3), ierr)
        call MPI_WAITALL(3, ireq3, istat3, ierr)
        call print_vtk(t2, nx, ny, nz, myrank, nranks, x, y, z, rho1d, p1d, v1d, ke0, entropy0)
      endif
    enddo

    deallocate(rho1d, p1d, v1d)
    if (mod(myrank,2) == 0) then
      deallocate(QJ, QJs, Rs, E, F, G, xix, etay, zetaz, Jacobian)
    endif
    if (kind(id_rescale) == 4) then
      if (myrank == 0 .or. myrank == rerank) then
        deallocate(Qre, Qm)
      elseif (myrank == rerank+1) then
        deallocate(Qm_cpu)
      endif
    endif
  end subroutine RungeKutta_4th
end module calc_time_dev

