module set_bc_common
  use mod_globals, only : nx, ny, nz
  implicit none
  interface set_bc_cyclic
    module procedure set_bc_cyclic2_init, set_bc_cyclic2, set_bc_cyclic4, &
                     set_bc_cyclic4_init, set_bc_cyclic6, set_bc_cyclic6_init
  end interface set_bc_cyclic
contains
  !> Cyclic boundary condition initialization for second-order accuracy
  !> Executed on CPU before main time-stepping loop
  subroutine set_bc_cyclic2_init(id_accuracy, nx, ny, nz, Q)
    integer(kind=2), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout)             :: Q(nx,5,ny,nz)
    integer i, j, k
    do k = 2, nz-1
      do j = 2, ny-1
        Q(1,:,j,k)  = Q(nx-1,:,j,k)
        Q(nx,:,j,k) = Q(2,:,j,k)
    enddo;enddo
    do k = 2, nz-1
      do i = 2, nx-1
        Q(i,:,1,k)  = Q(i,:,ny-1,k)
        Q(i,:,ny,k) = Q(i,:,2,k)
    enddo;enddo
    do k = 2, nz-1
      Q(1,:,1,k)   = Q(nx-1,:,ny-1,k)
      Q(nx,:,1,k)  = Q(2,:,ny-1,k)
      Q(1,:,ny,k)  = Q(nx-1,:,2,k)
      Q(nx,:,ny,k) = Q(2,:,2,k)
    enddo
    do j = 1, ny
      do i = 1, nx
        Q(i,:,j,1)  = Q(i,:,j,nz-1)
        Q(i,:,j,nz) = Q(i,:,j,2)
    enddo;enddo
  end subroutine set_bc_cyclic2_init 

  !> Cyclic boundary condition initialization for 4th-order accuracy
  !> Executed on CPU before main time-stepping loop
  subroutine set_bc_cyclic4_init(id_accuracy, nx, ny, nz, Q)
    integer(kind=4), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout)             :: Q(nx,5,ny,nz)
    integer i, j, k
    do k = 3, nz-2
      do j = 3, ny-2
        Q(1:2,:,j,k) = Q(nx-3:nx-2,:,j,k)
        Q(nx-1:nx,:,j,k) = Q(3:4,:,j,k)
    enddo;enddo
    do k = 3, nz-2
      do i = 3, nx-2
        Q(i,:,1:2,k) = Q(i,:,ny-3:ny-2,k)
        Q(i,:,ny-1:ny,k) = Q(i,:,3:4,k)
    enddo;enddo
    do k = 3, nz-2
      Q(1:2,:,1:2,k) = Q(nx-3:nx-2,:,ny-3:ny-2,k)
      Q(nx-1:nx,:,1:2,k) = Q(3:4,:,ny-3:ny-2,k)
      Q(1:2,:,ny-1:ny,k) = Q(nx-3:nx-2,:,3:4,k)
      Q(nx-1:nx,:,ny-1:ny,k) = Q(3:4,:,3:4,k)
    enddo
    do j = 1, ny
      do i = 1, nx
        Q(i,:,j,1:2) = Q(i,:,j,nz-3:nz-2)
        Q(i,:,j,nz-1:nz) = Q(i,:,j,3:4)
    enddo;enddo
  end subroutine set_bc_cyclic4_init

  !> Cyclic boundary condition initialization for 6th-order accuracy
  !> Executed on CPU before main time-stepping loop
  subroutine set_bc_cyclic6_init(id_accuracy, nx, ny, nz, Q)
    integer(kind=8), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout)             :: Q(nx,5,ny,nz)
    integer i, j, k
    do k = 4, nz-3
      do j = 4, ny-3
        Q(1:3,:,j,k) = Q(nx-5:nx-3,:,j,k)
        Q(nx-2:nx,:,j,k) = Q(4:6,:,j,k)
    enddo;enddo
    do k = 4, nz-3
      do i = 4, nx-3
        Q(i,:,1:3,k) = Q(i,:,ny-5:ny-3,k)
        Q(i,:,ny-2:ny,k) = Q(i,:,4:6,k)
    enddo;enddo
    do k = 4, nz-3
      Q(1:3,:,1:3,k) = Q(nx-5:nx-3,:,ny-5:ny-3,k)
      Q(nx-2:nx,:,1:3,k) = Q(4:6,:,ny-5:ny-3,k)
      Q(1:3,:,ny-2:ny,k) = Q(nx-5:nx-3,:,4:6,k)
      Q(nx-2:nx,:,ny-2:ny,k) = Q(4:6,:,4:6,k)
    enddo
    do j = 1, ny
      do i = 1, nx
        Q(i,:,j,1:3) = Q(i,:,j,nz-5:nz-3)
        Q(i,:,j,nz-2:nz) = Q(i,:,j,4:6)
    enddo;enddo
  end subroutine set_bc_cyclic6_init

  !> Cyclic boundary condition for second-order accuracy
  !> Executed on GPU during main time-stepping loop
  subroutine set_bc_cyclic2(id_accuracy, nx, ny, nz, Q)
    integer(kind=2), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout), device     :: Q(nx,5,ny,nz)
    integer i, j, k, l
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        do l = 1, 5
          Q(1,l,j,k)  = Q(nx-1,l,j,k)
          Q(nx,l,j,k) = Q(2,l,j,k)
    enddo;enddo;enddo
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do i = 2, nx-1
        do l = 1, 5
          Q(i,l,1,k)  = Q(i,l,ny-1,k)
          Q(i,l,ny,k) = Q(i,l,2,k)
    enddo;enddo;enddo
    !$cuf kernel do(1)<<<*,*>>>
    do k = 2, nz-1
      do l = 1, 5
        Q(1,l,1,k)   = Q(nx-1,l,ny-1,k)
        Q(nx,l,1,k)  = Q(2,l,ny-1,k)
        Q(1,l,ny,k)  = Q(nx-1,l,2,k)
        Q(nx,l,ny,k) = Q(2,l,2,k)
    enddo;enddo
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        do l = 1, 5
          Q(i,l,j,1)  = Q(i,l,j,nz-1)
          Q(i,l,j,nz) = Q(i,l,j,2)
    enddo;enddo;enddo
  end subroutine set_bc_cyclic2

  !> Cyclic boundary condition for 4th-order accuracy
  !> Executed on GPU during main time-stepping loop
  subroutine set_bc_cyclic4(id_accuracy, nx, ny, nz, Q)
    integer(kind=4), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout), device     :: Q(nx,5,ny,nz)
    integer i, j, k, l
    !$cuf kernel do(2) <<<*,*>>>
    do k = 3, nz-2
      do j = 3, ny-2
        do l = 1, 5
          Q(1,l,j,k) = Q(nx-3,l,j,k)
          Q(2,l,j,k) = Q(nx-2,l,j,k)
          Q(nx-1,l,j,k) = Q(3,l,j,k)
          Q(nx,l,j,k) = Q(4,l,j,k)
    enddo;enddo;enddo
    !$cuf kernel do(2) <<<*,*>>>
    do k = 3, nz-2
      do i = 3, nx-2
        do l = 1, 5
          Q(i,l,1,k) = Q(i,l,ny-3,k)
          Q(i,l,2,k) = Q(i,l,ny-2,k)
          Q(i,l,ny-1,k) = Q(i,l,3,k)
          Q(i,l,ny,k) = Q(i,l,4,k)
    enddo;enddo;enddo
    !$cuf kernel do(1) <<<*,*>>>
    do k = 3, nz-2
      do l = 1, 5
        Q(1,l,1,k) = Q(nx-3,l,ny-3,k)
        Q(1,l,2,k) = Q(nx-3,l,ny-2,k)
        Q(2,l,1,k) = Q(nx-2,l,ny-3,k)
        Q(2,l,2,k) = Q(nx-2,l,ny-2,k)
        Q(nx-1,l,1,k) = Q(3,l,ny-3,k)
        Q(nx-1,l,2,k) = Q(3,l,ny-2,k)
        Q(nx,l,1,k)   = Q(4,l,ny-3,k)
        Q(nx,l,2,k)   = Q(4,l,ny-2,k)
        Q(1,l,ny-1,k) = Q(nx-3,l,3,k)
        Q(1,l,ny,k)   = Q(nx-3,l,4,k)
        Q(2,l,ny-1,k) = Q(nx-2,l,3,k)
        Q(2,l,ny,k)   = Q(nx-2,l,4,k)
        Q(nx-1,l,ny-1,k) = Q(3,l,3,k)
        Q(nx-1,l,ny,k)   = Q(3,l,4,k)
        Q(nx,l,ny-1,k)   = Q(4,l,3,k)
        Q(nx,l,ny,k)     = Q(4,l,4,k)
    enddo;enddo
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        do l = 1, 5
          Q(i,l,j,1) = Q(i,l,j,nz-3)
          Q(i,l,j,2) = Q(i,l,j,nz-2)
          Q(i,l,j,nz-1) = Q(i,l,j,3)
          Q(i,l,j,nz) = Q(i,l,j,4)
    enddo;enddo;enddo
  end subroutine set_bc_cyclic4

  !> Cyclic boundary condition for 6th-order accuracy
  !> Executed on GPU during main time-stepping loop
  subroutine set_bc_cyclic6(id_accuracy, nx, ny, nz, Q)
    integer(kind=8), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout), device     :: Q(nx,5,ny,nz)
    integer i, j, k, l
    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz-3
      do j = 4, ny-3
        do l = 1, 5
          Q(1,l,j,k) = Q(nx-5,l,j,k)
          Q(2,l,j,k) = Q(nx-4,l,j,k)
          Q(3,l,j,k) = Q(nx-3,l,j,k)
          Q(nx-2,l,j,k) = Q(4,l,j,k)
          Q(nx-1,l,j,k) = Q(5,l,j,k)
          Q(nx,l,j,k)   = Q(6,l,j,k)
    enddo;enddo;enddo
    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz-3
      do i = 4, nx-3
        do l = 1, 5
          Q(i,l,1,k) = Q(i,l,ny-5,k)
          Q(i,l,2,k) = Q(i,l,ny-4,k)
          Q(i,l,3,k) = Q(i,l,ny-3,k)
          Q(i,l,ny-2,k) = Q(i,l,4,k)
          Q(i,l,ny-1,k) = Q(i,l,5,k)
          Q(i,l,ny,k)   = Q(i,l,6,k)
    enddo;enddo;enddo
    !$cuf kernel do(1)<<<*,*>>>
    do k = 4, nz-3
      do l = 1, 5
        Q(1,l,1,k) = Q(nx-5,l,ny-5,k)
        Q(1,l,2,k) = Q(nx-5,l,ny-4,k)
        Q(1,l,3,k) = Q(nx-5,l,ny-3,k)
        Q(2,l,1,k) = Q(nx-4,l,ny-5,k)
        Q(2,l,2,k) = Q(nx-4,l,ny-4,k)
        Q(2,l,3,k) = Q(nx-4,l,ny-3,k)
        Q(3,l,1,k) = Q(nx-3,l,ny-5,k)
        Q(3,l,2,k) = Q(nx-3,l,ny-4,k)
        Q(3,l,3,k) = Q(nx-3,l,ny-3,k)
        Q(nx-2,l,1,k) = Q(4,l,ny-5,k)
        Q(nx-2,l,2,k) = Q(4,l,ny-4,k)
        Q(nx-2,l,3,k) = Q(4,l,ny-3,k)
        Q(nx-1,l,1,k) = Q(5,l,ny-5,k)
        Q(nx-1,l,2,k) = Q(5,l,ny-4,k)
        Q(nx-1,l,3,k) = Q(5,l,ny-3,k)
        Q(nx,l,1,k)   = Q(6,l,ny-5,k)
        Q(nx,l,2,k)   = Q(6,l,ny-4,k)
        Q(nx,l,3,k)   = Q(6,l,ny-3,k)
        Q(1,l,ny-2,k) = Q(nx-5,l,4,k)
        Q(1,l,ny-1,k) = Q(nx-5,l,5,k)
        Q(1,l,ny,k)   = Q(nx-5,l,6,k)
        Q(2,l,ny-2,k) = Q(nx-4,l,4,k)
        Q(2,l,ny-1,k) = Q(nx-4,l,5,k)
        Q(2,l,ny,k)   = Q(nx-4,l,6,k)
        Q(3,l,ny-2,k) = Q(nx-3,l,4,k)
        Q(3,l,ny-1,k) = Q(nx-3,l,5,k)
        Q(3,l,ny,k)   = Q(nx-3,l,6,k)
        Q(nx-2,l,ny-2,k) = Q(4,l,4,k)
        Q(nx-2,l,ny-1,k) = Q(4,l,5,k)
        Q(nx-2,l,ny,k)   = Q(4,l,6,k)
        Q(nx-1,l,ny-2,k) = Q(5,l,4,k)
        Q(nx-1,l,ny-1,k) = Q(5,l,5,k)
        Q(nx-1,l,ny,k)   = Q(5,l,6,k)
        Q(nx,l,ny-2,k)   = Q(6,l,4,k)
        Q(nx,l,ny-1,k)   = Q(6,l,5,k)
        Q(nx,l,ny,k)     = Q(6,l,6,k)
    enddo;enddo
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        do l = 1, 5
          Q(i,l,j,1) = Q(i,l,j,nz-5)
          Q(i,l,j,2) = Q(i,l,j,nz-4)
          Q(i,l,j,3) = Q(i,l,j,nz-3)
          Q(i,l,j,nz-2) = Q(i,l,j,4)
          Q(i,l,j,nz-1) = Q(i,l,j,5)
          Q(i,l,j,nz)   = Q(i,l,j,6)
    enddo;enddo;enddo
  end subroutine set_bc_cyclic6

  !> Cyclic boundary condition in z-direction for 6th-order accuracy
  !> Executed on GPU during main time-stepping loop
  subroutine set_bc_cyclic_z(nx, ny, nz, QJ)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: QJ(nx,5,ny,nz)
    integer i, j, l
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        do l = 1, 5
          QJ(i,l,j,1) = QJ(i,l,j,nz-5)
          QJ(i,l,j,2) = QJ(i,l,j,nz-4)
          QJ(i,l,j,3) = QJ(i,l,j,nz-3)
          QJ(i,l,j,nz-2) = QJ(i,l,j,4)
          QJ(i,l,j,nz-1) = QJ(i,l,j,5)
          QJ(i,l,j,nz)   = QJ(i,l,j,6)
    enddo;enddo;enddo
  end subroutine set_bc_cyclic_z

  !> Boundary condition for SGS viscosity and turbulent kinetic energy
  !> Executed on GPU during main time-stepping loop
  subroutine set_bc_mut_common(nx, ny, nz, mut, qc2)
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
  end subroutine set_bc_mut_common
end module set_bc_common

