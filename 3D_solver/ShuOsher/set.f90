module set
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, gamma, R, rhol, ul, pl
  implicit none
contains
  subroutine set_grid(nx,ny,nz,x,y,z,dx,dy,dz)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx), dy(ny), dz(nz)
    integer i, j, k
    real(8) dx1, dy1, dz1
    dx1 = Lx / dble(nx-1)
    dy1 = Ly / dble(ny-1)
    dz1 = Lz / dble(nz-1)
    x(1) = 0.d0
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo
    y(1) = 0.d0
    do j = 1, ny-1
      dy(j) = dy1
      y(j+1) = y(j) + dy(j)
    enddo
    z(1) = 0.d0
    do k = 1, nz-1
      dz(k) = dz1
      z(k+1) = z(k) + dz(k)
    enddo
  end subroutine set_grid
  
  subroutine set_init(nx,ny,nz,x,y,z,Q)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer i, j, k
    real(8) :: pi = acos(-1.d0)
    do i = 1, nx
      if (x(i) < 0.125d0) then
        Q(i,:,:,1) = rhol
        Q(i,:,:,2) = rhol * ul
        Q(i,:,:,3) = 0.d0
        Q(i,:,:,4) = 0.d0
        Q(i,:,:,5) = pl / (gamma - 1.d0) + 0.5d0 * rhol * ul**2
      else
        Q(i,:,:,1) = 1.d0 + 0.2d0 * sin(8.d0 * 2.d0 * pi * x(i))
        Q(i,:,:,2) = 0.d0
        Q(i,:,:,3) = 0.d0
        Q(i,:,:,4) = 0.d0
        Q(i,:,:,5) = 1.d0 / (gamma - 1.d0)
      endif
    enddo
  end subroutine set_init
  
  subroutine set_bc(nx,ny,nz,Jacobian,QJ,Qre)
    use mod_globals, only : id_accuracy
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: Jacobian(ny)
    real(8), intent(inout), device :: QJ(nx,ny,nz,5) ! Q / J
    real(8), intent(in), device    :: Qre(2,ny,nz,5)
    integer i, j, k, l, jc, kc, offset
    real(8), device :: Qc(nx,5)
    real(8) pw
    jc     = 4
    kc     = 4
    offset = 3

    ! inlet and outlet
    !$cuf kernel do(2) <<<*,*>>>
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        QJ(1,j,k,1)    = rhol / Jacobian(j)
        QJ(1,j,k,2)    = QJ(1,j,k,1) * ul
        QJ(1,j,k,3)    = 0.d0
        QJ(1,j,k,4)    = 0.d0
        QJ(1,j,k,5)    = (pl / (gamma - 1.d0) + 0.5d0 * rhol * ul**2) / Jacobian(j)

        QJ(nx,j,k,1)   = QJ(nx-1,j,k,1)
        QJ(nx,j,k,2)   = 0.d0
        QJ(nx,j,k,3)   = 0.d0
        QJ(nx,j,k,4)   = 0.d0
        pw             = (gamma - 1.d0) * (QJ(nx-1,j,k,5) - 0.5d0 * (QJ(nx-1,j,k,2)) / QJ(nx-1,j,k,1))
        QJ(nx,j,k,5)   = pw / (gamma - 1.d0)
    enddo;enddo

    !$cuf kernel do(2)<<<*,*>>>
    do l = 1, 5
      do i = 1, nx
        Qc(i,l) = QJ(i,jc,kc,l)
    enddo;enddo

    !$cuf kernel do(4)<<<*,*>>>
    do l = 1, 5
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            QJ(i,j,k,l) = Qc(i,l)
    enddo;enddo;enddo;enddo
  end subroutine set_bc

  subroutine set_bc_mut(nx,ny,nz,mut,qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2) <<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        mut(1,j,k) = mut(nx-1,j,k)
        mut(nx,j,k) = mut(2,j,k)
        qc2(1,j,k) = qc2(nx-1,j,k)
        qc2(nx,j,k) = qc2(2,j,k)
    enddo;enddo

    !$cuf kernel do(2) <<<*,*>>>
    do k = 2, nz-1
      do i = 2, nx-1
        mut(i,1,k) = mut(i,ny-1,k)
        mut(i,ny,k) = mut(i,2,k)
        qc2(i,1,k) = qc2(i,ny-1,k)
        qc2(i,ny,k) = qc2(i,2,k)
    enddo;enddo

    !$cuf kernel do(1) <<<*,*>>>
    do k = 2, nz-1
      mut(1,1,k) = mut(nx-1,ny-1,k)
      mut(nx,1,k) = mut(2,ny-1,k)
      mut(1,ny,k) = mut(nx-1,2,k)
      mut(nx,ny,k) = mut(2,2,k)
      qc2(1,1,k) = qc2(nx-1,ny-1,k)
      qc2(nx,1,k) = qc2(2,ny-1,k)
      qc2(1,ny,k) = qc2(nx-1,2,k)
      qc2(nx,ny,k) = qc2(2,2,k)
    enddo

    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        mut(i,j,1) = mut(i,j,nz-1)
        mut(i,j,nz) = mut(i,j,2)
        qc2(i,j,1) = qc2(i,j,nz-1)
        qc2(i,j,nz) = qc2(i,j,2)
    enddo;enddo
  end subroutine set_bc_mut
end module set

