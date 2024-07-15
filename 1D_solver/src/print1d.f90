module print1d
  use mod_globals, only : id_scheme, nt, dt, gamma, R, Lx, rho0, p0, T0, mu0
  implicit none
contains
  subroutine print1(step,nx,x,Q,sensor)
    integer, intent(in) :: step, nx
    real(4), intent(in) :: x(nx)
    real(4), intent(in) :: Q(nx,3)
    real(4), intent(in) :: sensor(nx-1)
    integer i
    real(4) time, T, M, sensorm
    real(4), dimension(nx) :: rho, u, p
    character(len=40) filename, scheme
    time = real(nt * step) * real(dt)
    do i = 1, nx
      rho(i) = Q(i,1)
      u(i)   = Q(i,2) / rho(i)
      p(i)   = (real(gamma) - 1.e0) * (Q(i,3) - 0.5e0 * rho(i) * u(i)**2)
    enddo

    if (id_scheme == 0) then
      scheme = "/KEEP/"
    elseif (id_scheme == 1) then
      scheme = "/KEP/"
    elseif (id_scheme == 2) then
      scheme = "/KEEPKEP/"
    endif

    ! rho, u, p, T, Mach, mu
    write(filename, "(a, i5.5, a)") "data/Q", int(step), ".d"
    open(10,file=filename)
    do i = 1, nx
      T = p(i) / (real(R) * rho(i))
      M = u(i) / sqrt(real(gamma * R) * T)
      if (2 <= i .and. i <= nx-1) then
        sensorm = 0.5e0 * (sensor(i-1) + sensor(i))
      else
        sensorm = 0.e0
      endif
      write(10,"(7e12.4)") x(i) / Lx, rho(i) / rho0, u(i) / sqrt(p0 / rho0), p(i) / p0, T / T0, M, sensorm
    enddo
    close(10)
    ! kinetic energy
    write(filename, "(a)") "data/kinetic_energy.d"
    open(10,file=filename, position="append")
    write(10,"(2e12.4)") time, sum(rho(:) * u(:)**2)
    close(10)
    ! entropy
    write(filename, "(a)") "data/entropy.d"
    open(10,file=filename, position="append")
    write(10,"(2e12.4)") time, sum(rho(:) * log(p(:) * (rho(:)**real(-gamma))))
    close(10)
  end subroutine print1
end module print1d

