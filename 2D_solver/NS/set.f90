module set
  use mod_globals, only : id_accuracy, nx, ny, gamma, R
  use mod_constant, only : Cp, over_gamma_1
  use set_bc_common
  use set_coordinate
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(1)
    integer i, j
    dx(:) = Lx / dble(nx-1)
    dy(:) = Ly / dble(ny-1)
    x(1) = 0.d0
    do i = 1, nx-1
      x(i+1) = x(i) + dx(i)
    enddo
    y(1) = -0.5d0 * Ly
    do j = 1, ny-1
      y(j+1) = y(j) + dy(j)
    enddo
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    use mod_globals, only : rhod, ud, pd, rhou, uu, pu
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(4,nx,ny)
    integer i, j
    do j = 1, ny
      do i = 1, int(nx * 0.25)
        Q(1,i,j) = rhod
        Q(2,i,j) = rhod * ud
        Q(3,i,j) = 0.d0
        Q(4,i,j) = pd * over_gamma_1 + 0.5d0 * rhod * ud**2
      enddo
      do i = int(nx * 0.25), nx
        Q(1,i,j) = rhou
        Q(2,i,j) = rhou * uu
        Q(3,i,j) = 0.d0
        Q(4,i,j) = pu * over_gamma_1 + 0.5d0 * rhou * uu**2
      enddo
    enddo
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, x, y, Jacobian, QJ)
    use mod_globals, only : id_accuracy, rhod, ud, pd, rhou, uu, pu
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: x(nx), y(ny)
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ(4,nx,ny)
    real(8) over_Jacobian
    integer i, j
    !$cuf kernel do(1)<<<*,*>>>
    do j = 1, ny
      ! shock-downstream
      over_Jacobian = 1.d0 / Jacobian(1,j)
      QJ(1,1,j) = rhod * over_Jacobian
      QJ(2,1,j) = rhod * ud * over_Jacobian
      QJ(3,1,j) = 0.d0
      QJ(4,1,j) = (pd * over_gamma_1 + 0.5d0 * rhod * ud**2) * over_Jacobian
      ! shock-upstream
      over_Jacobian = 1.d0 / Jacobian(nx,j)
      QJ(1,nx,j) = rhou * over_Jacobian
      QJ(2,nx,j) = rhou * uu * over_Jacobian
      QJ(3,nx,j) = 0.d0
      QJ(4,nx,j) = (pu * over_gamma_1 + 0.5d0 * rhou * uu**2) * over_Jacobian
    enddo
    if (kind(id_accuracy) == 2) then
      !$cuf kernel do(1)<<<*,*>>>
      do i = 2, nx-1
        do j = 1, 4
          QJ(j,i,1)  = QJ(j,i,ny-1)
          QJ(j,i,ny) = QJ(j,i,2)
      enddo;enddo
    elseif (kind(id_accuracy) == 4) then
      !$cuf kernel do(1)<<<*,*>>>
      do i = 2, nx-1
        do j = 1, 4
          QJ(j,i,1)  = QJ(j,i,ny-2)
          QJ(j,i,2)  = QJ(j,i,ny-1)
          QJ(j,i,ny-1) = QJ(j,i,3)
          QJ(j,i,ny)   = QJ(j,i,4)
      enddo;enddo
    else
      !$cuf kernel do(1)<<<*,*>>>
      do i = 2, nx-1
        do j = 1, 4
          QJ(j,i,1)  = QJ(j,i,ny-3)
          QJ(j,i,2)  = QJ(j,i,ny-2)
          QJ(j,i,3)  = QJ(j,i,ny-1)
          QJ(j,i,ny-2) = QJ(j,i,4)
          QJ(j,i,ny-1) = QJ(j,i,5)
          QJ(j,i,ny)   = QJ(j,i,6)
      enddo;enddo
    endif
  end subroutine set_bc


  attributes(global) subroutine calc_force(nx, ny, x, y, dx, dy, Q, Fout)
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: x(nx), y(ny), dx(nx-1), dy(ny-1)
    real(8), intent(in), device  :: Q(4,nx,ny)
    real(8), intent(out), device :: Fout(3,nx-2,ny-2)
  end subroutine calc_force
end module set

