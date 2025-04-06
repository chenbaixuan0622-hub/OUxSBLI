module calc_time_dev
  use mpi
  use mod_globals, only : gamma, Pr, R_gas, nt, np, dx, dt
  use calc_muscl
  use calc_step
  use set
  use print
  implicit none
contains
  subroutine calc_RHS_rho(nx, rho_m, u, RHS_rho)
    integer, intent(in)  :: nx
    real(8), intent(in)  :: rho_m(nx-1), u(nx+1)
    real(8), intent(out) :: RHS_rho(nx-2)
    integer i
    do i = 1, nx-2
      RHS_rho(i) = (- rho_m(i)   * u(i+1) &
                    + rho_m(i+1) * u(i+2)) / dx
    enddo
  end subroutine calc_RHS_rho


  subroutine calc_RHS_u(nx, rho_m, u, p, RHS_u)
    integer, intent(in)  :: nx
    real(8), intent(in)  :: rho_m(nx-1), u(nx+1), p(nx)
    real(8), intent(out) :: RHS_u(nx-1)
    real(8) rho2(2)
    integer i
    do i = 1, nx-1
      RHS_u(i) = 0.5d0 * (0.5d0 * (u(i)   + u(i+1)) * (-u(i)   + u(i+1)) / dx &
                        + 0.5d0 * (u(i+1) + u(i+2)) * (-u(i+1) + u(i+2)) / dx) &
                 + (-p(i) + p(i+1)) / (dx * rho_m(i))
    enddo
  end subroutine calc_RHS_u


  subroutine calc_RHS_p(nx, u, p, RHS_p)
    integer, intent(in)  :: nx
    real(8), intent(in)  :: u(nx+1), p(nx)
    real(8), intent(out) :: RHS_p(nx-2)
    integer i
    do i = 1, nx-2
      RHS_p(i) = 0.5d0 * (u(i+1) * (-p(i)   + p(i+1)) / dx &
                        + u(i+2) * (-p(i+1) + p(i+2)) / dx) &
                + gamma * p(i+1) * (-u(i+1) + u(i+2)) / dx
    enddo
  end subroutine calc_RHS_p


  function Sutherland(T) result(mu)
    real(8), intent(in) :: T
    real(8) :: mu, mu0 = 1.716d-5, T0 = 273.2d0, S = 111.d0
    mu = mu0 * ((T0 + S) / (T + S)) * (T / T0) ** 1.5d0
  end function Sutherland

  
  subroutine calc_tau(nx, rho, u, p, tau, kTx)
    integer, intent(in)  :: nx
    real(8), intent(in)  :: rho(nx), u(nx+1), p(nx)
    real(8), intent(out) :: tau(nx), kTx(nx-1)
    integer i
    real(8) :: Cp = gamma * R_gas / (gamma - 1.d0)
    real(8) mu, T(nx), kappa(nx)
    do i = 1, nx
      T(i) = p(i) / (R_gas * rho(i))
      mu = Sutherland(T(i))
      tau(i) = 4.d0 * mu * (-u(i) + u(i+1)) / (3.d0 * dx)
      kappa(i) = mu * Cp / Pr
    enddo
    do i = 1, nx-1
      kTx(i) = 0.5d0 * (kappa(i) + kappa(i+1)) * (-T(i) + T(i+1)) / dx
    enddo
  end subroutine calc_tau


  subroutine calc_RHS_u_visc(nx, rho_m, tau, RHS_u)
    integer, intent(in)    :: nx
    real(8), intent(in)    :: rho_m(nx-1), tau(nx)
    real(8), intent(inout) :: RHS_u(nx-1)
    integer i
    do i = 1, nx-1
      RHS_u(i) = RHS_u(i) - (-tau(i) + tau(i+1)) / (dx * rho_m(i))
    enddo
  end subroutine calc_RHS_u_visc


  subroutine calc_RHS_p_visc(nx, u, tau, kTx, RHS_p)
    integer, intent(in)    :: nx
    real(8), intent(in)    :: u(nx+1), tau(nx), kTx(nx-1)
    real(8), intent(inout) :: RHS_p(nx-2)
    integer i
    do i = 1, nx-2
      RHS_p(i) = RHS_p(i) - (gamma - 1.d0) * (tau(i+1) * (-u(i+1) + u(i+2)) / dx &
                                              + 0.5d0 * (kTx(i) + kTx(i+1)))
    enddo
  end subroutine calc_RHS_p_visc


  subroutine calc_R(nx, rho, u, p, RHS_rho, RHS_u, RHS_p)
    integer, intent(in)  :: nx
    real(8), intent(in)  :: rho(nx), u(nx+1), p(nx)
    real(8), intent(out) :: RHS_rho(nx-2), RHS_u(nx-1), RHS_p(nx-2)
    real(8) rho_m(nx-1)
    integer stat, i, j
    real(8) tau(nx), kTx(nx-1)
    
    call calc_rho_muscl(nx, rho, rho_m)
    call calc_RHS_rho(nx, rho_m, u, RHS_rho)
    call calc_RHS_u(nx, rho_m, u, p, RHS_u)
    call calc_RHS_p(nx, u, p, RHS_p)
    ! for viscous flow
    call calc_tau(nx, rho, u, p, tau, kTx)
    call calc_RHS_u_visc(nx, rho_m, tau, RHS_u)
    call calc_RHS_p_visc(nx, u, tau, kTx, RHS_p)
  end subroutine calc_R


  subroutine RungeKutta(myrank, nx, x, rho, u, p)
    integer, intent(in)    :: myrank, nx
    real(8), intent(in)    :: x(nx)
    real(8), intent(inout) :: rho(nx), u(nx+1), p(nx)
    integer t1, t2, ndevices, ilen, ierr, stat, ireq, istat(MPI_STATUS_SIZE)
    real(8), allocatable :: rho2(:), u2(:), p2(:) 
    real(8), allocatable :: RHS_rho(:), RHS_u(:), RHS_p(:) 

    if (myrank == 0) then
      allocate(rho2(nx), u2(nx+1), p2(nx))
      allocate(RHS_rho(nx-2), RHS_u(nx-1), RHS_p(nx-2))
      RHS_rho(:) = 0.d0
      RHS_u(:)   = 0.d0
      RHS_p(:)   = 0.d0

      call print_1d(0, nx, real(x), real(rho), real(u), real(p))
    endif

    do t2 = 1, np
      if (myrank == 0) then
        do t1 = 1, nt
          call calc_R(nx, rho, u, p, RHS_rho, RHS_u, RHS_p)
          call calc_step1(nx, dt, rho, RHS_rho, rho2)
          call calc_step1(nx+1, dt, u,   RHS_u,   u2)
          call calc_step1(nx,   dt, p,   RHS_p,   p2)
          call set_bc(nx, rho2, u2, p2)

          call calc_R(nx, rho2, u2, p2, RHS_rho, RHS_u, RHS_p)
          call calc_step2(nx, dt, rho, RHS_rho, rho2)
          call calc_step2(nx+1, dt, u,   RHS_u,   u2)
          call calc_step2(nx,   dt, p,   RHS_p,   p2)
          call set_bc(nx, rho2, u2, p2)

          call calc_R(nx, rho2, u2, p2, RHS_rho, RHS_u, RHS_p)
          call calc_step3(nx, dt, rho2, RHS_rho, rho)
          call calc_step3(nx+1, dt, u2,   RHS_u,   u)
          call calc_step3(nx,   dt, p2,   RHS_p,   p)
          call set_bc(nx, rho, u, p)
        enddo
      endif

      ! send and recv device arrays
      if (myrank == 0) then
        call MPI_SEND(rho, nx, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
        call MPI_SEND(u, nx+1, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
        call MPI_SEND(p,   nx, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
      elseif (myrank == 1) then
        call MPI_RECV(rho, nx, MPI_REAL8, 0, 0, MPI_COMM_WORLD, istat, ierr)
        call MPI_RECV(u, nx+1, MPI_REAL8, 0, 0, MPI_COMM_WORLD, istat, ierr)
        call MPI_RECV(p,   nx, MPI_REAL8, 0, 0, MPI_COMM_WORLD, istat, ierr)
        call print_1d(t2, nx, real(x), real(rho), real(u), real(p))
      endif
    enddo
    
    if (myrank == 0) then
      deallocate(rho2, u2, p2, RHS_rho, RHS_u, RHS_p)
    endif
  end subroutine RungeKutta
end module calc_time_dev

