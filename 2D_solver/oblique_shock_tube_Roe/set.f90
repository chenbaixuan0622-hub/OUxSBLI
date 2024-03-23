module set
  use cudafor
  use mod_globals, only : nx, ny, gamma
  implicit none
contains
  subroutine set_init(Q)
    use mod_globals, only : rhol, rhor, pl, pr
    real(8), intent(out), dimension(nx,ny,4) :: Q
    integer i, j
    do j = 1, ny
      do i = 1, nx
        if (i > j) then
          Q(i,j,1) = rhol
          Q(i,j,4) = pl / (gamma - 1.d0)
        else
          Q(i,j,1) = rhor
          Q(i,j,4) = pr / (gamma - 1.d0)
        endif
      enddo
    enddo
    ! both (velocity = 0)
    Q(:,:,2) = 0.d0
    Q(:,:,3) = 0.d0
  end subroutine set_init

  subroutine set_bc(Q)
    real(8), intent(inout), device :: Q(nx,ny,4)
    integer i, j, k
    !$cuf kernel do <<<*,*>>>
    do k = 1, 4
      do j = 2, ny-1
        ! inlet
        Q(1,j,k) = Q(2,j,k)
        ! outlet
        Q(nx,j,k) = Q(nx-1,j,k)
      enddo
    enddo

    !$cuf kernel do <<<*,*>>>
    do k = 1, 4
      do i = 1, nx
        ! bottom-wall
        Q(i,1,k) = Q(i,2,k)
        ! top-wall
        Q(i,ny,k) = Q(i,ny-1,k)
      enddo
    enddo
  end subroutine set_bc
end module set

