module set
  use mod_globals, only : nx, ny, gamma, R, rho0, p0, rho1, p1
  use mod_constant, only : id_accuracy, gamma_1, over_gamma_1
  use set_bc_common
  use set_coordinate
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, Lx, Ly, x, y, dx, dy)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: Lx, Ly
    real(8), intent(out) :: x(nx), y(ny), dx(nx-1), dy(ny-1)
    real(8) dx0, dy0
    integer i, j
    dx0 = Lx / dble(nx-1)
    dy0 = Ly / dble(ny-1)
    dx  = dx0
    dy  = dy0
    x(1) = 0.d0
    do i = 1, nx-1
      x(i+1) = x(i) + dx0
    enddo
    y(1) = 0.d0
    do j = 1, ny-1
      y(j+1) = y(j) + dy0
    enddo
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, x, y, Q)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: x(nx), y(ny)
    real(8), intent(out) :: Q(nx,4,ny)
    integer i, j
    do j = 1, ny
      do i = 1, nx / 2
        Q(i,1,j) = rho0
        Q(i,2,j) = 0.d0
        Q(i,3,j) = 0.d0
        Q(i,4,j) = p0 / (gamma - 1.d0)
      enddo
      do i = nx / 2 + 1, nx
        Q(i,1,j) = rho1
        Q(i,2,j) = 0.d0
        Q(i,3,j) = 0.d0
        Q(i,4,j) = p1 / (gamma - 1.d0)
      enddo
    enddo
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, Jacobian, QJ)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ(nx,4,ny)
    real(8) p_wall
    integer i, j, l
    ! inlet outlet
    !$cuf kernel do<<<*,*>>>
    do j = 2, ny-1
      do l = 1, 4
        QJ(1,l,j)  = QJ(2,l,j)
        QJ(nx,l,j) = QJ(nx-1,l,j)
    enddo;enddo

    ! no wall
    !!$cuf kernel do<<<*,*>>>
    !do i = 1, nx
    !  do l = 1, 4
    !    QJ(i,l,1)  = QJ(i,l,2)
    !    QJ(i,l,ny) = QJ(i,l,ny-1)
    !enddo;enddo

    ! no-slip wall & symetric boundary condition
    !$cuf kernel do<<<*,*>>>
    do i = 1, nx
      ! no-slip wall
      p_wall = gamma_1 * (QJ(i,4,2) - 0.5d0 * (QJ(i,2,2)**2 + QJ(i,3,2)**2) / QJ(i,1,2))
      QJ(i,1,1) = QJ(i,1,2)
      QJ(i,2,1) = 0.d0
      QJ(i,3,1) = 0.d0
      QJ(i,4,1) = p_wall * over_gamma_1
      ! symetric boundary condition
      QJ(i,1,ny) =  QJ(i,1,ny-1)
      QJ(i,2,ny) =  QJ(i,2,ny-1)
      QJ(i,3,ny) = -QJ(i,3,ny-1)
      QJ(i,4,ny) =  QJ(i,4,ny-1)
    enddo
  end subroutine set_bc
end module set

