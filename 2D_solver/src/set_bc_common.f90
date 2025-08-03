module set_bc_common
  use mod_globals, only : nx, ny, nz
  implicit none
  interface set_bc_cyclic
    module procedure set_bc_cyclic2, set_bc_cyclic4, set_bc_cyclic6
  end interface set_bc_cyclic
contains
  subroutine set_bc_cyclic2_init(id_accuracy, nx, ny, Q)
    integer(kind=2), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny
    real(8), intent(inout)             :: Q(4,nx,ny)
    integer i, j
    do j = 2, ny-1
      Q(:,1,j)  = Q(:,nx-1,j)
      Q(:,nx,j) = Q(:,2,j)
    enddo
    do i = 2, nx-1
      Q(:,i,1)  = Q(:,i,ny-1)
      Q(:,i,ny) = Q(:,i,2)
    enddo
    Q(:,1,1)   = Q(:,nx-1,ny-1)
    Q(:,nx,1)  = Q(:,2,ny-1)
    Q(:,1,ny)  = Q(:,nx-1,2)
    Q(:,nx,ny) = Q(:,2,2)
  end subroutine set_bc_cyclic2_init 

  subroutine set_bc_cyclic4_init(id_accuracy, nx, ny, Q)
    integer(kind=4), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny
    real(8), intent(inout)             :: Q(4,nx,ny)
    integer i, j
    do j = 3, ny-2
      Q(:,1:2,j) = Q(:,nx-3:nx-2,j)
      Q(:,nx-1:nx,j) = Q(:,3:4,j)
    enddo
    do i = 3, nx-2
      Q(:,i,1:2) = Q(:,i,ny-3:ny-2)
      Q(:,i,ny-1:ny) = Q(:,i,3:4)
    enddo
    Q(:,1:2,1:2) = Q(:,nx-3:nx-2,ny-3:ny-2)
    Q(:,nx-1:nx,1:2) = Q(:,3:4,ny-3:ny-2)
    Q(:,1:2,ny-1:ny) = Q(:,nx-3:nx-2,3:4)
    Q(:,nx-1:nx,ny-1:ny) = Q(:,3:4,3:4)
  end subroutine set_bc_cyclic4_init
  
  subroutine set_bc_cyclic6_init(id_accuracy, nx, ny, Q)
    integer(kind=8), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny
    real(8), intent(inout)             :: Q(4,nx,ny)
    integer i, j
    do j = 4, ny-3
      Q(:,1:3,j) = Q(:,nx-5:nx-3,j)
      Q(:,nx-2:nx,j) = Q(:,4:6,j)
    enddo
    do i = 4, nx-3
      Q(:,i,1:3) = Q(:,i,ny-5:ny-3)
      Q(:,i,ny-2:ny) = Q(:,i,4:6)
    enddo
    Q(:,1:3,1:3) = Q(:,nx-5:nx-3,ny-5:ny-3)
    Q(:,nx-2:nx,1:3) = Q(:,4:6,ny-5:ny-3)
    Q(:,1:3,ny-2:ny) = Q(:,nx-5:nx-3,4:6)
    Q(:,nx-2:nx,ny-2:ny) = Q(:,4:6,4:6)
  end subroutine set_bc_cyclic6_init
  
  subroutine set_bc_cyclic2(id_accuracy, nx, ny, Q)
    integer(kind=2), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny
    real(8), intent(inout), device     :: Q(4,nx,ny)
    integer i, j, l
    !$cuf kernel do(2)<<<*,*>>>
    do j = 2, ny-1
      do l = 1, 4
        Q(l,1,j)  = Q(l,nx-1,j)
        Q(l,nx,j) = Q(l,2,j)
    enddo;enddo
    !$cuf kernel do(2)<<<*,*>>>
    do i = 2, nx-1
      do l = 1, 4
        Q(l,i,1)  = Q(l,i,ny-1)
        Q(l,i,ny) = Q(l,i,2)
    enddo;enddo
    !$cuf kernel do(1)<<<*,*>>>
    do l = 1, 4
      Q(l,1,1)   = Q(l,nx-1,ny-1)
      Q(l,nx,1)  = Q(l,2,ny-1)
      Q(l,1,ny)  = Q(l,nx-1,2)
      Q(l,nx,ny) = Q(l,2,2)
    enddo
  end subroutine set_bc_cyclic2

  subroutine set_bc_cyclic4(id_accuracy, nx, ny, Q)
    integer(kind=4), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny
    real(8), intent(inout), device     :: Q(4,nx,ny)
    integer i, j, l
    !$cuf kernel do(2) <<<*,*>>>
    do j = 3, ny-2
      do l = 1, 4
        Q(l,1,j) = Q(l,nx-3,j)
        Q(l,2,j) = Q(l,nx-2,j)
        Q(l,nx-1,j) = Q(l,3,j)
        Q(l,nx,j) = Q(l,4,j)
    enddo;enddo
    !$cuf kernel do(2) <<<*,*>>>
    do i = 3, nx-2
      do l = 1, 4
        Q(l,i,1) = Q(l,i,ny-3)
        Q(l,i,2) = Q(l,i,ny-2)
        Q(l,i,ny-1) = Q(l,i,3)
        Q(l,i,ny) = Q(l,i,4)
    enddo;enddo
    !$cuf kernel do(1) <<<*,*>>>
    do l = 1, 4
        Q(l,1,1) = Q(l,nx-3,ny-3)
        Q(l,1,2) = Q(l,nx-3,ny-2)
        Q(l,2,1) = Q(l,nx-2,ny-3)
        Q(l,2,2) = Q(l,nx-2,ny-2)
        Q(l,nx-1,1) = Q(l,3,ny-3)
        Q(l,nx-1,2) = Q(l,3,ny-2)
        Q(l,nx,1)   = Q(l,4,ny-3)
        Q(l,nx,2)   = Q(l,4,ny-2)
        Q(l,1,ny-1) = Q(l,nx-3,3)
        Q(l,1,ny)   = Q(l,nx-3,4)
        Q(l,2,ny-1) = Q(l,nx-2,3)
        Q(l,2,ny)   = Q(l,nx-2,4)
        Q(l,nx-1,ny-1) = Q(l,3,3)
        Q(l,nx-1,ny)   = Q(l,3,4)
        Q(l,nx,ny-1)   = Q(l,4,3)
        Q(l,nx,ny)     = Q(l,4,4)
    enddo
  end subroutine set_bc_cyclic4

  subroutine set_bc_cyclic6(id_accuracy, nx, ny, Q)
    integer(kind=8), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny
    real(8), intent(inout), device     :: Q(4,nx,ny)
    integer i, j, l
    !$cuf kernel do(2)<<<*,*>>>
    do j = 4, ny-3
      do l = 1, 4
        Q(l,1,j) = Q(l,nx-5,j)
        Q(l,2,j) = Q(l,nx-4,j)
        Q(l,3,j) = Q(l,nx-3,j)
        Q(l,nx-2,j) = Q(l,4,j)
        Q(l,nx-1,j) = Q(l,5,j)
        Q(l,nx,j)   = Q(l,6,j)
    enddo;enddo
    !$cuf kernel do(2)<<<*,*>>>
    do i = 4, nx-3
      do l = 1, 4
        Q(l,i,1) = Q(l,i,ny-5)
        Q(l,i,2) = Q(l,i,ny-4)
        Q(l,i,3) = Q(l,i,ny-3)
        Q(l,i,ny-2) = Q(l,i,4)
        Q(l,i,ny-1) = Q(l,i,5)
        Q(l,i,ny)   = Q(l,i,6)
    enddo;enddo
    !$cuf kernel do(1)<<<*,*>>>
    do l = 1, 4
      Q(l,1,1) = Q(l,nx-5,ny-5)
      Q(l,1,2) = Q(l,nx-5,ny-4)
      Q(l,1,3) = Q(l,nx-5,ny-3)
      Q(l,2,1) = Q(l,nx-4,ny-5)
      Q(l,2,2) = Q(l,nx-4,ny-4)
      Q(l,2,3) = Q(l,nx-4,ny-3)
      Q(l,3,1) = Q(l,nx-3,ny-5)
      Q(l,3,2) = Q(l,nx-3,ny-4)
      Q(l,3,3) = Q(l,nx-3,ny-3)
      Q(l,nx-2,1) = Q(l,4,ny-5)
      Q(l,nx-2,2) = Q(l,4,ny-4)
      Q(l,nx-2,3) = Q(l,4,ny-3)
      Q(l,nx-1,1) = Q(l,5,ny-5)
      Q(l,nx-1,2) = Q(l,5,ny-4)
      Q(l,nx-1,3) = Q(l,5,ny-3)
      Q(l,nx,1)   = Q(l,6,ny-5)
      Q(l,nx,2)   = Q(l,6,ny-4)
      Q(l,nx,3)   = Q(l,6,ny-3)
      Q(l,1,ny-2) = Q(l,nx-5,4)
      Q(l,1,ny-1) = Q(l,nx-5,5)
      Q(l,1,ny)   = Q(l,nx-5,6)
      Q(l,2,ny-2) = Q(l,nx-4,4)
      Q(l,2,ny-1) = Q(l,nx-4,5)
      Q(l,2,ny)   = Q(l,nx-4,6)
      Q(l,3,ny-2) = Q(l,nx-3,4)
      Q(l,3,ny-1) = Q(l,nx-3,5)
      Q(l,3,ny)   = Q(l,nx-3,6)
      Q(l,nx-2,ny-2) = Q(l,4,4)
      Q(l,nx-2,ny-1) = Q(l,4,5)
      Q(l,nx-2,ny)   = Q(l,4,6)
      Q(l,nx-1,ny-2) = Q(l,5,4)
      Q(l,nx-1,ny-1) = Q(l,5,5)
      Q(l,nx-1,ny)   = Q(l,5,6)
      Q(l,nx,ny-2)   = Q(l,6,4)
      Q(l,nx,ny-1)   = Q(l,6,5)
      Q(l,nx,ny)     = Q(l,6,6)
    enddo
  end subroutine set_bc_cyclic6
end module set_bc_common

