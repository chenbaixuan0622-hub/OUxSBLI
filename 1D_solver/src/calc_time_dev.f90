module calc_time_dev
  use cudafor
  use mpi
  use mod_globals, only : id_visc, id_LL, gamma, Pr, R, nt, np, dt, dx, over_dx, dtdx, blocks, threads
  use mod_constant, only : one_third, gamma_1, over_gamma_1, R_over_gamma_1, Cp, Cp_over_Pr, mu0_T0_S_over_T0_2_3
  use set
  use print
  implicit none
contains
  attributes(global) subroutine calc_E(nx, rho, u, p, T, E)
    integer, intent(in), value   :: nx
    real(8), intent(in), device  :: rho(nx), u(nx), p(nx), T(nx)
    real(8), intent(out), device :: E(3,nx-1)
    integer i
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (i > nx-1) return
    E(1,i) = 0.25d0 * (rho(i) + rho(i+1)) * (u(i) + u(i+1))
    E(2,i) = E(1,i) * 0.5d0 * (u(i) + u(i+1)) + 0.5d0 * (p(i) + p(i+1))
    E(3,i) = 0.5d0 * E(1,i) * u(i) * u(i+1) &
           + E(1,i) * 0.5d0 * (T(i) + T(i+1)) * R_over_gamma_1 &
           + 0.5d0 * (u(i) * p(i+1) + u(i+1) * p(i))
  end subroutine calc_E

  
  attributes(global) subroutine calc_Ev(nx, rho, u, p, T, mu, E)
    integer, intent(in), value   :: nx
    real(8), intent(in), device  :: rho(nx), u(nx), p(nx), T(nx), mu(nx)
    real(8), intent(out), device :: E(3,nx-1)
    real(8) mua, txx, utxx, kTx
    integer i
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (i > nx-1) return
    mua  = 0.5d0 * (mu(i) + mu(i+1))
    txx  = 4.d0 * mua * (-u(i) + u(i+1)) * one_third * over_dx
    utxx = 0.5d0 * (u(i) + u(i+1)) * txx
    kTx  = Cp_over_Pr * mua * (-T(i) + T(i+1)) * over_dx
    E(2,i) = E(2,i) - txx
    E(3,i) = E(3,i) - (utxx + kTx)
  end subroutine calc_Ev


  attributes(global) subroutine calc_Ev_LL(nx, rho, u, p, T, mu, Z, E)
    integer, intent(in), value   :: nx
    real(8), intent(in), device  :: rho(nx), u(nx), p(nx), T(nx), mu(nx)
    real(8), intent(in), device  :: Z(2,nx)  
    real(8), intent(out), device :: E(3,nx-1)
    real(8) :: kb = 1.380650d-23
    real(8) mua, txx, utxx, kappa, kTx, s, q
    integer i
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (i > nx-1) return
    mua   = 0.5d0 * (mu(i) + mu(i+1))
    txx   = 4.d0 * mua * (-u(i) + u(i+1)) * one_third * over_dx
    utxx  = 0.5d0 * (u(i) + u(i+1)) * txx
    kappa = Cp_over_Pr * mua
    kTx   = kappa * (-T(i) + T(i+1)) * over_dx
    s     = sqrt(4.d0 * kb * mua * (T(i) + T(i+1)) / (3.d0 * dt * dx)) * 0.5d0 * (Z(1,i) + Z(1,i+1))
    q     = sqrt(kb * kappa * (T(i)**2 + T(i+1)**2) / (dt * dx)) * 0.5d0 * (Z(2,i) + Z(2,i+1))
    E(2,i) = E(2,i) - (txx + sqrt(2.d0) * s)
    E(3,i) = E(3,i) - (utxx + kTx + sqrt(2.d0) * (q + 0.5d0 * (u(i) + u(i+1)) * s))
  end subroutine calc_Ev_LL


  subroutine calc_R(nx, Q, Z, R)
    integer, intent(in), value     :: nx
    real(8), intent(in), device    :: Q(3,nx)
    real(8), intent(inout), device :: Z(2,nx)
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
    call calc_E<<<blocks,threads,0>>>(nx, rho, u, p, T, E)
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


  subroutine calc_step1(nx, Q, R, Q2)
    integer, intent(in), value   :: nx
    real(8), intent(in), device  :: Q(3,nx)
    real(8), intent(in), device  :: R(3,nx-2)
    real(8), intent(out), device :: Q2(3,nx)
    integer i, j
    !$cuf kernel do(1)<<<*,*>>>
    do i = 2, nx-1
      do j = 1, 3
        Q2(j,i) = Q(j,i) - R(j,i-1)
    enddo;enddo
  end subroutine calc_step1


  subroutine calc_step2(nx, Q, R, Q2)
    integer, intent(in), value     :: nx
    real(8), intent(in), device    :: Q(3,nx)
    real(8), intent(in), device    :: R(3,nx-2)
    real(8), intent(inout), device :: Q2(3,nx)
    integer i, j
    !$cuf kernel do(1)<<<*,*>>>
    do i = 2, nx-1
      do j = 1, 3
        Q2(j,i) = 0.25d0 * (3.d0 * Q(j,i) + Q2(j,i) - R(j,i-1))
    enddo;enddo
  end subroutine calc_step2
  

  subroutine calc_step3(nx, Q2, R, Q)
    integer, intent(in), value     :: nx
    real(8), intent(in), device    :: Q2(3,nx)
    real(8), intent(in), device    :: R(3,nx-2)
    real(8), intent(inout), device :: Q(3,nx)
    integer i, j
    !$cuf kernel do(1)<<<*,*>>>
    do i = 2, nx-1
      do j = 1, 3
        Q(j,i) = (Q(j,i) + 2.d0 * Q2(j,i) - 2.d0 * R(j,i-1)) * one_third
    enddo;enddo
  end subroutine calc_step3
  

  subroutine RungeKutta(myrank, nx, x, Q_cpu)
    integer, intent(in)    :: myrank, nx
    real(8), intent(in)    :: x(nx)
    real(8), intent(inout) :: Q_cpu(3,nx)
    integer t1, t2, ndevices, ilen, ierr, stat, ireq, istat(MPI_STATUS_SIZE)
    type(cudaDeviceProp)         :: prop
    real(8), allocatable, device :: Q(:,:), Q2(:,:), Z(:,:), R(:,:)

    if (myrank == 0) then
      stat = cudaGetDeviceCount(ndevices)
      print '(2x, i2, a)', ndevices, " GPU devices are found"
      stat = cudaSetDevice(myrank)
      stat = cudaGetDeviceProperties(prop,myrank)
      ilen = verify(prop%name, ' ', .true.)
      print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", myrank, ") is available"

      allocate(Q(3,nx), Q2(3,nx), Z(2,nx), R(3,nx-2))

      call print_1d(0, nx, real(x), real(Q_cpu))
      Q = Q_cpu
    endif

    do t2 = 1, np
      if (myrank == 0) then
        do t1 = 1, nt
          call calc_R(nx, Q, Z, R)
          call calc_step1(nx, Q, R, Q2)
          call set_bc(nx, Q2)

          call calc_R(nx, Q2, Z, R)
          call calc_step2(nx, Q, R, Q2)
          call set_bc(nx, Q2)

          call calc_R(nx, Q2, Z, R)
          call calc_step3(nx, Q2, R, Q)
          call set_bc(nx, Q)
        enddo
      endif

      ! send and recv device arrays
      if (myrank == 0) then
        Q_cpu = Q
        call MPI_SEND(Q_cpu, nx*3, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
      elseif (myrank == 1) then
        call MPI_RECV(Q_cpu, nx*3, MPI_REAL8, 0, 0, MPI_COMM_WORLD, istat, ierr)
        call print_1D(t2, nx, real(x), real(Q_cpu))
      endif
    enddo
    
    if (myrank == 0) then
      deallocate(Q, Q2, Z, R)
    endif
  end subroutine RungeKutta
end module calc_time_dev

