module calc_time_dev
  use cudafor
  use mpi
  use mod_globals, only : accuracy, id_visc, nt, np, &
  & blocks, threads, blocksE, threadsE
  use calc_physical_quantities
  use calc_steps
  use calc_flux
  use calc_visc
  use set
  use print1d
  implicit none
  interface RungeKutta
    module procedure RungeKutta_3rd, RungeKutta_4th
  end interface
contains
  subroutine calc_EFG(nx,dx,QJ,E,sensor)
    integer, intent(in), value                   :: nx
    real(8), intent(in), dimension(nx-1), device :: dx ! 1 / dx
    real(8), intent(in), dimension(nx,3), device :: QJ ! Q / Jacobian
    real(8), intent(out), device                 :: E(nx-accuracy+1,3), sensor(nx-1)
    real(8), dimension(nx), device :: rho, u, p, T
    integer stat
    call calc_quantities(nx,QJ,rho,u,p,T)

    call calc_E<<<blocksE,threadsE>>>(nx,rho,u,p,E,sensor)

    if (id_visc == 1) then
      call calc_Ev<<<blocksE,threadsE>>>(nx,dx,rho,u,T,p,E)
    endif
  end subroutine calc_EFG

  subroutine RungeKutta_3rd(id_RungeKutta,myrank,nx,x,dx_cpu,xix_cpu,Q)
    integer(kind=2), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, nx
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1), xix_cpu(nx-1)
    real(8), intent(inout)      :: Q(nx,3)
    integer i, t1, t2, itr, ilen, ierr, stat, request, status(MPI_STATUS_SIZE)
    real(8), dimension(nx-1)    :: sensor_cpu
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)          :: prop
    real(8), allocatable, device  :: QJ(:,:), QJ2(:,:), QJ3(:,:), E(:,:)
    real(8), allocatable, device  :: dx(:), xix(:), sensor(:) 

    ! check GPU
    stat = cudaSetDevice(0)
    stat = cudaGetDeviceProperties(prop,0)
    ilen = verify(prop%name, ' ', .true.)
    print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", 0, ") is available"

    if (myrank == 0) then
      allocate(QJ(nx,3),QJ2(nx,3),QJ3(nx,3),E(nx-1,3),dx(nx-1),xix(nx-1),sensor(nx-1))
      sensor_cpu(:) = 0.d0

      ! print initial condition
      call print1(0,nx,real(x),real(Q),real(sensor_cpu))

      ! copy on GPU
      QJ     = Q
      xix    = xix_cpu
      dx     = dx_cpu
      sensor = sensor_cpu
    endif

    do t2 = 1, np
      if (myrank == 0) then
        do t1 = 1, nt
          call calc_EFG(nx,xix,QJ,E,sensor)
          call calc_step(nx,1.d0,0.d0,xix,E,QJ,QJ2)
          call set_bc(nx,QJ2)

          call calc_EFG(nx,xix,QJ2,E,sensor)
          call calc_step2(nx,0.75d0,0.25d0,0.25d0,1.d0,xix,E,QJ,QJ2,QJ3)
          call set_bc(nx,QJ3)

          call calc_EFG(nx,xix,QJ3,E,sensor)
          call calc_step3(nx,xix,E,QJ3,QJ)
          call set_bc(nx,QJ)
        enddo
      endif

      ! send and recv device arrays
      if (myrank == 0) then
        Q = QJ
        sensor_cpu = sensor
        call MPI_SEND(Q,          nx*3, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
        call MPI_SEND(sensor_cpu, nx-1, MPI_REAL8, 1, 1, MPI_COMM_WORLD, ierr) 
      elseif (myrank == 1) then
        call MPI_RECV(Q,          nx*3, MPI_REAL8, 0, 0, MPI_COMM_WORLD, status, ierr)
        call MPI_RECV(sensor_cpu, nx*3, MPI_REAL8, 0, 1, MPI_COMM_WORLD, status, ierr)
        call print1(t2,nx,real(x),real(Q),real(sensor_cpu))
      endif
    enddo
    
    if (myrank == 0) then
      deallocate(QJ,QJ2,QJ3,E,dx,xix,sensor)
    endif
  end subroutine RungeKutta_3rd

  subroutine RungeKutta_4th(id_RungeKutta,myrank,nx,x,dx_cpu,xix_cpu,Q)
    integer(kind=4), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, nx
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1), xix_cpu(nx-1)
    real(8), intent(inout)      :: Q(nx,3)
    integer i, j, k, t1, t2, itr, ierr, ilen, stat, status(MPI_STATUS_SIZE)
    real(8), dimension(nx-1)    :: sensor_cpu
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)          :: prop
    real(8), allocatable, device  :: QJ(:,:), QJs(:,:), Rs(:,:), E(:,:)
    real(8), allocatable, device  :: dx(:), xix(:), sensor(:)

    ! check GPU
    stat = cudaSetDevice(0)
    stat = cudaGetDeviceProperties(prop,0)
    ilen = verify(prop%name, ' ', .true.)

    if (myrank == 0) then
      allocate(QJ(nx,3),QJs(nx,3),Rs(nx,3),E(nx-1,3),dx(nx-1),xix(nx-1),sensor(nx-1))
      sensor_cpu(:) = 0.d0

      ! print initial condition
      call print1(0,nx,real(x),real(Q),real(sensor_cpu))

      ! copy on GPU
      QJ      = Q
      Rs(:,:) = 0.d0
      xix     = xix_cpu
      dx      = dx_cpu
      sensor  = sensor_cpu
    endif

    do t2 = 1, np
      if (myrank == 0) then
        do t1 = 1, nt
          call calc_EFG(nx,xix,QJ,E,sensor)
          call calc_step(nx,0.5d0,1.d0,xix,E,QJ,QJs,Rs) ! QJs = Q2
          call set_bc(nx,QJs)

          call calc_EFG(nx,xix,QJs,E,sensor)
          call calc_step(nx,0.5d0,2.d0,xix,E,QJ,QJs,Rs) ! QJs = Q3
          call set_bc(nx,QJs)

          call calc_EFG(nx,xix,QJs,E,sensor)
          call calc_step(nx,1.0d0,2.d0,xix,E,QJ,QJs,Rs) ! QJs = Q4
          call set_bc(nx,QJs)

          call calc_EFG(nx,xix,QJs,E,sensor)
          call calc_step4(nx,xix,E,Rs,QJ)
          call set_bc(nx,QJ)
        enddo
      endif

      ! send and recv device arrays
      if (myrank == 0) then
        Q = QJ
        sensor_cpu = sensor
        call MPI_SEND(Q,          nx*3, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
        call MPI_SEND(sensor_cpu, nx-1, MPI_REAL8, 1, 1, MPI_COMM_WORLD, ierr) 
      elseif (myrank == 1) then
        call MPI_RECV(Q,          nx*3, MPI_REAL8, 0, 0, MPI_COMM_WORLD, status, ierr)
        call MPI_RECV(sensor_cpu, nx*3, MPI_REAL8, 0, 1, MPI_COMM_WORLD, status, ierr)
        call print1(t2,nx,real(x),real(Q),real(sensor_cpu))
      endif
    enddo

    if (myrank == 0) then
      deallocate(QJ,QJs,Rs,E,dx,xix,sensor)
    endif
  end subroutine RungeKutta_4th
end module calc_time_dev

