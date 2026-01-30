module calc_time_dev
  use cudafor
  use mpi
  use mod_globals, only : id_visc, id_LL, id_igr, nt, np, dt, dx, over_dx, dtdx, blocks, threads
  use mod_constant, only : gamma_1, mu0_T0_S_over_T0_2_3
  use calc_flux
  use calc_visc
  use calc_igr
  use calc_steps
  use set
  use print
  implicit none
contains
  subroutine calc_R(nx, Q, Z, sigma, R)
    integer, intent(in), value     :: nx
    real(8), intent(in), device    :: Q(3,nx)
    real(8), intent(inout), device :: Z(2,nx)
    real(4), intent(inout), device :: sigma(nx)
    real(8), intent(out), device   :: R(3,nx-2)
    real(8), dimension(nx), device :: rho, u, p, T, mu
    integer stat, i, j
    real(8), dimension(3,nx-1), device :: E
    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      rho(i) = Q(1,i)
      u(i)   = Q(2,i) / rho(i)
      p(i)   = gamma_1 * (Q(3,i) - 0.5d0 * rho(i) * u(i)**2)
      T(i)   = p(i) / (rho(i) * 287.03d0)
      mu(i)  = mu0_T0_S_over_T0_2_3 / (T(i) + 111.d0) * T(i)**1.5d0
    enddo
    if (kind(id_igr) == 4) then
      call calc_sigma(nx, Q, sigma)
      !$cuf kernel do(1)<<<*,*>>>
      do i = 2, nx-1
        if (-u(i-1) + u(i+1) < 0.d0) then
          p(i) = p(i) + dble(sigma(i))
        endif
      enddo
    endif
    call calc_KEEP_G<<<blocks,threads,0>>>(nx, rho, u, p, T, sigma, E)
    if (kind(id_visc) == 4 .and. kind(id_LL) == 2) then
      call calc_Ev<<<blocks,threads,1>>>(nx, rho, u, p, T, mu, E)
    elseif (kind(id_visc) == 4 .and. kind(id_LL) == 4) then
      call calc_Ev_LL<<<blocks,threads,1>>>(nx, rho, u, p, T, mu, Z, E)
    endif
    stat = cudaDeviceSynchronize() 
    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx-2
      do j = 1, 3
        R(j,i) = dtdx * (-E(j,i) + E(j,i+1))
    enddo;enddo
  end subroutine calc_R


  subroutine RungeKutta(myrank, nx, x, Q_cpu)
    integer, intent(in)    :: myrank, nx
    real(8), intent(in)    :: x(nx)
    real(8), intent(inout) :: Q_cpu(3,nx)
    integer t1, t2, ndevices, ilen, ierr, stat, ireq, istat(MPI_STATUS_SIZE)
    type(cudaDeviceProp)         :: prop
    real(8), allocatable, device :: Q(:,:), Q2(:,:), Z(:,:), R(:,:)
    real(4), allocatable, device :: sigma(:)
    real(4), allocatable         :: sigma_cpu(:)

    allocate(sigma_cpu(nx))
    if (myrank == 0) then
      stat = cudaGetDeviceCount(ndevices)
      print '(2x, i2, a)', ndevices, " GPU devices are found"
      stat = cudaSetDevice(myrank)
      stat = cudaGetDeviceProperties(prop,myrank)
      ilen = verify(prop%name, ' ', .true.)
      print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", myrank, ") is available"

      allocate(Q(3,nx), Q2(3,nx), Z(2,nx), sigma(nx), R(3,nx-2))

      sigma     = 0.e0
      sigma_cpu = 0.e0
      call print_1d(0, nx, real(x), real(Q_cpu), sigma_cpu)
      Q = Q_cpu
    endif
    do t2 = 1, np
      if (myrank == 0) then
        do t1 = 1, nt
          call calc_R(nx, Q, Z, sigma, R)
          call calc_step1(nx, Q, R, Q2)
          call set_bc(nx, Q2)

          call calc_R(nx, Q2, Z, sigma, R)
          call calc_step2(nx, Q, R, Q2)
          call set_bc(nx, Q2)

          call calc_R(nx, Q2, Z, sigma, R)
          call calc_step3(nx, Q2, R, Q)
          call set_bc(nx, Q)
        enddo
      endif

      ! send and recv device arrays
      if (myrank == 0) then
        Q_cpu     = Q
        sigma_cpu = sigma
        call MPI_SEND(Q_cpu,   nx*3, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
        call MPI_SEND(sigma_cpu, nx, MPI_REAL4, 1, 1, MPI_COMM_WORLD, ierr) 
      elseif (myrank == 1) then
        call MPI_RECV(Q_cpu,   nx*3, MPI_REAL8, 0, 0, MPI_COMM_WORLD, istat, ierr)
        call MPI_RECV(sigma_cpu, nx, MPI_REAL4, 0, 1, MPI_COMM_WORLD, istat, ierr)
        call print_1D(t2, nx, real(x), real(Q_cpu), sigma_cpu)
      endif
    enddo
    
    if (myrank == 0) then
      deallocate(Q, Q2, Z, R, sigma)
    endif
    deallocate(sigma_cpu)
  end subroutine RungeKutta
end module calc_time_dev

