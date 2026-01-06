module set
  use cufft
  use mod_globals, only : id_accuracy, nx, ny, gamma, R
  use mod_constant, only : Cp
  use set_bc_common
  use set_init_common
  use set_coordinate
  use calc_forcing
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(1)
    call set_grid_cyclic(id_accuracy, nx, ny, Lx, Ly, xc, yc, dx, dy)
  end subroutine set_grid

  
  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    use mod_globals, only : rho0, p0, urms
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(4,nx,ny)
    real(8), allocatable :: u(:,:), v(:,:), p(:,:)
    real(8) rms, fact
    integer i, j, offset, accuracy, k0
    if (kind(id_accuracy) == 8) then
      accuracy = 6
      offset   = 3
    elseif (kind(id_accuracy) == 4) then
      accuracy = 4
      offset   = 2
    elseif (kind(id_accuracy) == 2) then
      accuracy = 2
      offset   = 1
    endif
    allocate(u(nx-accuracy,ny-accuracy), v(nx-accuracy,ny-accuracy), p(nx-accuracy,ny-accuracy))
    k0 = 5
    call calc_spectra(accuracy, nx, ny, k0, 1.d0, rho0, p0, u, v, p)
    ! calc RMS
    rms = 0.d0
    do j = 1, ny-accuracy
      do i = 1, nx-accuracy
        rms = rms + u(i,j)*u(i,j) + v(i,j)*v(i,j)
    enddo;enddo
    rms = sqrt(rms / dble((nx-accuracy)*(ny-accuracy)))
    fact = urms / rms
    do j = 1+offset, ny-offset
      do i = 1+offset, nx-offset
        Q(1,i,j) = rho0
        Q(2,i,j) = rho0 * fact * u(i-offset,j-offset)
        Q(3,i,j) = rho0 * fact * v(i-offset,j-offset)
        Q(4,i,j) = (fact * p(i-offset,j-offset) + p0) / (gamma - 1.d0) &
                  + 0.5d0 * (Q(2,i,j)**2 + Q(3,i,j)**2) / Q(1,i,j)
    enddo;enddo
    call set_bc_cyclic_init(id_accuracy, nx, ny, Q)
    deallocate(u, v, p)
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, x, y, Jacobian, Q)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: x(nx), y(ny)
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: Q(4,nx,ny)
    call set_bc_cyclic(id_accuracy, nx, ny, Q)
  end subroutine set_bc


  subroutine calc_force_stream(nx, ny, x, y, dx, dy, Q, Fout, rand, phase)
    use mod_globals, only : rho0, urms, kf, kmin, kmax, Lx, Ly, dt, pi
    integer, intent(in), value     :: nx, ny
    real(8), intent(inout), device :: x(nx), y(ny)
    real(8), intent(in), device    :: dx(nx-1), dy(ny-1)
    real(8), intent(in), device    :: Q(4,nx,ny)
    real(8), intent(out), device   :: Fout(3,nx-2,ny-2)
    real(4), intent(in), device    :: rand(kmax,kmax)
    real(8), intent(inout), device :: phase(kmax,kmax)
    real(8), parameter :: tau   = 0.001d0 * Lx / (dble(kf) * urms)
    real(8), parameter :: expT  = exp(-dt / tau)
    real(8), parameter :: sigma = 2.d0 * pi
    real(8), device :: amp(kmax,kmax), u(nx,ny), v(nx,ny)
    real(8) kabs, phix, phiy, sum, coef
    integer i, j, kx, ky
    amp = 0.d0; sum = 0.d0
    !$cuf kernel do(2)<<<*,*>>>
    do ky = 1, kmax
      do kx = 1, kmax
        ! Ornstein-Uhlenbeck
        phase(kx,ky) = expT * phase(kx,ky) + sigma * sqrt(1.d0 - expT*expT) * dble(2.e0 * rand(kx,ky) - 1.e0)
        kabs = sqrt(dble(kx*kx + ky*ky))
        if (kabs < dble(kmin) .or. kabs > dble(kmax)) cycle
        amp(kx,ky) = exp(-(kabs/dble(kf))**2) / (kabs*kabs)
        sum = sum + amp(kx,ky)
    enddo;enddo
    u = 0.d0; v = 0.d0
    coef = urms / (sum * dble(kf))
    !$cuf kernel do(2)<<<*,*>>>
    do j = 2, ny-1
      do i = 2, nx-1
        phix = dble(i-2) * (2.d0 * pi / dble(nx-2))
        phiy = dble(j-2) * (2.d0 * pi / dble(ny-2))
        do ky = 1, kmax
          do kx = 1, kmax
            u(i,j) = u(i,j) + coef * amp(kx,ky) * dble(ky) * sin(dble(kx)*phix + dble(ky)*phiy + phase(kx,ky))
            v(i,j) = v(i,j) - coef * amp(kx,ky) * dble(kx) * sin(dble(kx)*phix + dble(ky)*phiy + phase(kx,ky))
        enddo;enddo
        Fout(1,i-1,j-1) = Q(1,i,j) * (u(i,j) - urms * 0.05d0 * Q(2,i,j))
        Fout(2,i-1,j-1) = Q(1,i,j) * (v(i,j) - urms * 0.05d0 * Q(3,i,j))
        Fout(3,i-1,j-1) = Q(2,i,j) * Fout(1,i-1,j-1) + Q(3,i,j) * Fout(2,i-1,j-1)
    enddo;enddo
  end subroutine calc_force_stream
  

  subroutine calc_force(nx, ny, x, y, dx, dy, Q, Fout)
    use mod_globals, only : rho0, urms, kf, kmin, kmax, Lx, Ly
    integer, intent(in), value     :: nx, ny
    real(8), intent(inout), device :: x(nx), y(ny)
    real(8), intent(in), device    :: dx(nx-1), dy(ny-1)
    real(8), intent(in), device    :: Q(4,nx,ny)
    real(8), intent(out), device   :: Fout(3,nx-2,ny-2)
    real(8), device :: F_x(nx-2,ny-2)
    real(8), device :: F_y(nx-2,ny-2)
    integer i, j
    call calc_forcing_main(2, plan_fwd, plan_inv, nx, ny, kmin, kmax, 0.1d0, Q, F_x, F_y)
    !$cuf kernel do(2)<<<*,*>>>
    do j = 2, ny-1
      do i = 2, nx-1
        Fout(1,i-1,j-1) = Q(1,i,j) * (F_x(i-1,j-1) - urms * 0.05d0 * Q(2,i,j))
        Fout(2,i-1,j-1) = Q(1,i,j) * (F_y(i-1,j-1) - urms * 0.05d0 * Q(3,i,j))
        Fout(3,i-1,j-1) = Q(2,i,j) * Fout(1,i-1,j-1) + Q(3,i,j) * Fout(2,i-1,j-1)
    enddo;enddo
  end subroutine calc_force
end module set

