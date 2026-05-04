!> 3rd-order TVD Runge-Kutta time integration for curvilinear solver
module calc_time_dev_curv
  use cudafor
  use mpi
  use mod_globals, only : id_visc, nt, np, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, dt
  use mod_constant, only : one_third
  use calc_flux_base_curv
  use calc_steps_curv
  use set
  use preprocess_curv
  use print_curv
  implicit none
contains
  !> 3rd-order TVD Runge-Kutta time stepping for curvilinear grids
  !> Solves dQ/dt = RHS(Q) using curvilinear flux computation
  subroutine RungeKutta_curv(id_RungeKutta, id_rescale, myrank, mygpu, &
                             nx, ny, nz, x_phys, y_phys, z, dz_val, &
                             Jac_cpu, n_xi_x_cpu, n_xi_y_cpu, n_eta_x_cpu, n_eta_y_cpu, &
                             xi_x_cpu, xi_y_cpu, eta_x_cpu, eta_y_cpu, Q)
    integer(2), intent(in) :: id_RungeKutta
    integer(2), intent(in) :: id_rescale
    integer, intent(in)    :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)    :: x_phys(nx,ny), y_phys(nx,ny), z(nz), dz_val
    real(8), intent(in)    :: Jac_cpu(nx,ny)
    real(8), intent(in)    :: n_xi_x_cpu(nx-1,ny-2), n_xi_y_cpu(nx-1,ny-2)
    real(8), intent(in)    :: n_eta_x_cpu(nx-2,ny-1), n_eta_y_cpu(nx-2,ny-1)
    real(8), intent(in)    :: xi_x_cpu(nx,ny), xi_y_cpu(nx,ny)
    real(8), intent(in)    :: eta_x_cpu(nx,ny), eta_y_cpu(nx,ny)
    real(8), intent(inout) :: Q(nx,5,ny,nz)
    integer :: i, j, k, l, t1, t2, overlap, ierr, nranks, ndevices, stat
    integer :: istat(MPI_STATUS_SIZE)
    real(8) :: dt_xi, dt_eta
    ! GPU device arrays
    real(8), allocatable, device :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJ2(:,:,:,:)
    real(8), allocatable, device :: E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: T(:,:,:), mu(:,:,:)
    real(8), allocatable, device :: dt_Szeta(:,:)
    real(8), allocatable, device :: n_xi_x(:,:), n_xi_y(:,:)
    real(8), allocatable, device :: n_eta_x(:,:), n_eta_y(:,:)
    real(8), allocatable, device :: xi_x(:,:), xi_y(:,:)
    real(8), allocatable, device :: eta_x(:,:), eta_y(:,:)
    real(8), allocatable, device :: Jacobian(:,:)
    ! For I/O
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0
    
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    stat = cudaGetDeviceCount(ndevices)
    print *, "rank", myrank, " has found ", ndevices, " GPU devices"
    
    if (mod(myrank,2) == 0) then
      call check_gpu(mygpu)
      call allocate_device_mem_curv(myrank, nx, ny, nz, &
          dt_Szeta, n_xi_x, n_xi_y, n_eta_x, n_eta_y, &
          xi_x, xi_y, eta_x, eta_y, Jacobian, &
          ruvwp, T, mu, E, F, G)
      allocate(QJ(nx,5,ny,nz), QJ2(nx,5,ny,nz), stat=ierr)
      if (ierr /= 0) then
        print *, "myrank is ", myrank, " memory allocation failed", ierr
      else
        print *, "myrank is ", myrank, " memory allocation has completed"
      endif
      
      call pre_calc_curv(nx, ny, nz, myrank, nranks, dz_val, &
          Jac_cpu, n_xi_x_cpu, n_xi_y_cpu, n_eta_x_cpu, n_eta_y_cpu, &
          xi_x_cpu, xi_y_cpu, eta_x_cpu, eta_y_cpu, &
          x_phys, y_phys, z, Q, overlap, &
          dt_Szeta, n_xi_x, n_xi_y, n_eta_x, n_eta_y, &
          xi_x, xi_y, eta_x, eta_y, Jacobian, QJ, ke0, entropy0)
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank, MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank, MPI_COMM_WORLD, istat, ierr)
    endif
    
    ! Compute dt_xi, dt_eta scalars: dt * dz (uniform z spacing)
    dt_xi  = dt * dz_val
    dt_eta = dt * dz_val
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    print *, "myrank is ", myrank, " start Runge-Kutta (curvilinear)"
    
    do t2 = 1, np
      if (mod(myrank,2) == 0) then
        do t1 = 1, nt
          ! Stage 1: Compute E, F, G from current state QJ
          call calc_EFG_curv(id_visc, nx, ny, nz, dz_val, &
              n_xi_x, n_xi_y, n_eta_x, n_eta_y, xi_x, xi_y, eta_x, eta_y, Jacobian, &
              QJ, ruvwp, T, mu, E, F, G)

          ! TVD RK3 Stage 1: Q(1) = Q^n - (dt/J) * (E_flux_div + F_flux_div + G_flux_div)
          call calc_step1_curv<<<blocks,threads>>>(nx, ny, nz, 1.d0, dt_xi, dt_eta, dt_Szeta, E, F, G, QJ, QJ2)
          call set_bc(myrank, nx, ny, nz, Jacobian, eta_x, eta_y, QJ2)

          ! Stage 2: Compute E, F, G from Q(1)
          call calc_EFG_curv(id_visc, nx, ny, nz, dz_val, &
              n_xi_x, n_xi_y, n_eta_x, n_eta_y, xi_x, xi_y, eta_x, eta_y, Jacobian, &
              QJ2, ruvwp, T, mu, E, F, G)

          ! TVD RK3 Stage 2: Q(2) = (3/4)*Q^n + (1/4)*Q(1) - (1/4)*(dt/J)*flux_div
          call calc_step2_3_curv<<<blocks,threads>>>(nx, ny, nz, 0.75d0, 0.25d0, 0.25d0, 1.d0, &
              dt_xi, dt_eta, dt_Szeta, E, F, G, QJ, QJ2)
          call set_bc(myrank, nx, ny, nz, Jacobian, eta_x, eta_y, QJ2)

          ! Stage 3: Compute E, F, G from Q(2)
          call calc_EFG_curv(id_visc, nx, ny, nz, dz_val, &
              n_xi_x, n_xi_y, n_eta_x, n_eta_y, xi_x, xi_y, eta_x, eta_y, Jacobian, &
              QJ2, ruvwp, T, mu, E, F, G)
          
          ! TVD RK3 Stage 3: Q^(n+1) = (1/3)*Q^n + (2/3)*Q(2) - (2/3)*(dt/J)*flux_div
          call calc_step2_3_curv<<<blocks,threads>>>(nx, ny, nz, 2.d0, 1.d0, 2.d0, one_third, &
              dt_xi, dt_eta, dt_Szeta, E, F, G, QJ2, QJ)
          call set_bc(myrank, nx, ny, nz, Jacobian, eta_x, eta_y, QJ)
        enddo
      endif
      
      ! I/O synchronization
      if (mod(myrank, 2) == 0) then
        call send_recv_for_print_even_curv(myrank, nranks, t2, nx, ny, nz, x_phys, y_phys, z, Jac_cpu, QJ, Q, ke0, entropy0)
      else
        call send_recv_for_print_odd_curv(myrank, nranks, t2, nx, ny, nz, x_phys, y_phys, z, Jac_cpu, Q, ke0, entropy0)
      endif
    enddo
    
    if (mod(myrank,2) == 0) then
      deallocate(ruvwp, T, mu, QJ, QJ2, E, F, G)
      deallocate(dt_Szeta, n_xi_x, n_xi_y, n_eta_x, n_eta_y, xi_x, xi_y, eta_x, eta_y, Jacobian)
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_curv
end module calc_time_dev_curv
