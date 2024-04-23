module set
  use cudafor
  use mod_globals, only : nx, ny, dx, dy, gamma, u0, R, p0, rho0
  implicit none
contains
  subroutine calc_Blasius(eta,d,u,v)
    real(8), intent(in), value :: eta, d
    real(8), intent(out) :: u, v
    real(8) f, df
    real(8) fs(9), dfs(9)
    real(8) :: nu0 = 3.8206d-5
    integer i
    fs(:) = (/0.d0, 0.165d0, 0.65d0, 1.397d0, 2.306d0, 3.283d0, 4.28d0, 5.279d0, 6.279d0/)
    dfs(:) = (/0.d0, 0.3298d0, 0.6298d0, 0.8461d0, 0.9555d0, 0.9915d0, 0.999d0, 0.9999d0, 1.d0/)
    do i = 1, 8
      if (dble(i-1) <= eta .and. eta <= dble(i)) then
        f = fs(i) + (fs(i+1) - fs(i)) * (eta - dble(i-1))
        df = dfs(i) + (dfs(i+1) - dfs(i)) * (eta - dble(i-1))
      elseif (8.d0 < eta .and. eta <= 8.8d0) then
        f = fs(9) + 1.25d0 * (7.07923d0 - fs(9)) * (eta - 8.d0)
        df = 1.d0
      elseif (8.8d0 < eta) then
        f = 7.07923d0
        df = 1.d0
      endif
    enddo
    u = u0 * df
    v = 0.5d0 * (nu0 / d) * (eta * df - f)
    write(*,*) eta, f, df
  end subroutine calc_Blasius

  subroutine set_init(Q)
    real(8), intent(out), dimension(nx,ny,4) :: Q
    real(8) y, eta, u, v, p_wall
    real(8) :: d = 1.d-3
    real(8) :: nu0 = 3.8206d-5
    integer i, j
    Q(:,:,1) = rho0
    ! boundary layer
    do j = 3, ny
      y = dble(j-2) * dy
      eta = y / d
      !u = min(u0, u0 * (0.0015d0 * eta**4 - 0.0181d0 etayd**3 + 0.029d0 etayd**2 + 0.3192d0 etayd + 0.0003d0)) 
      call calc_Blasius(eta,d,u,v)
      Q(:,j,2) = rho0 * u
      Q(:,j,3) = rho0 * v
      Q(:,j,4) = p0 / (gamma - 1.d0) + 0.5d0 * (Q(:,j,2)**2 + Q(:,j,3)**2) / Q(:,j,1) 
    enddo

    ! noSlip
    do i = 1, nx
      ! wall
      Q(i,2,1) = Q(i,3,1)
      Q(i,2,2) = Q(i,3,2)
      Q(i,2,3) = 0.d0
      p_wall = (gamma - 1.d0) * (Q(i,3,4) - 0.5d0 * (Q(i,3,2)**2 + Q(i,3,3)**2) / Q(i,3,1))
      Q(i,2,4) = p_wall / (gamma - 1.d0)
      ! imaginray
      Q(i,1,1) = Q(i,3,1)
      Q(i,1,2) = Q(i,3,2) 
      Q(i,1,3) = -Q(i,3,3) 
      Q(i,1,4) = Q(i,3,4) 
    enddo
  end subroutine set_init

  subroutine set_bc(Q)
    real(8), intent(inout), device :: Q(nx,ny,4)
    integer i, j, k
    real(8) p_wall
    !$cuf kernel do <<<*,*>>>
    do k = 1, 4
      do j = 3, ny-1
        ! inlet
        Q(1,j,k) = Q(nx-3,j,k)
        Q(2,j,k) = Q(nx-2,j,k)
        ! outlet
        Q(nx-1,j,k) = Q(3,j,k)
        Q(nx,j,k) = Q(4,j,k)
      enddo
    enddo

    !$cuf kernel do <<<*,*>>>
    do i = 1, nx
      ! top
      Q(i,ny,1) = Q(i,ny-1,1)!rho0
      Q(i,ny,2) = Q(i,ny-1,1) * u0
      Q(i,ny,3) = Q(i,ny-1,3)!0.d0
      Q(i,ny,4) = Q(i,ny-1,4)!p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
    enddo
    
    ! noSlip
    !$cuf kernel do <<<*,*>>>
    do i = 1, nx
      ! wall
      Q(i,2,1) = Q(i,3,1)
      Q(i,2,2) = Q(i,3,2)
      Q(i,2,3) = 0.d0
      p_wall = (gamma - 1.d0) * (Q(i,3,4) - 0.5d0 * (Q(i,3,2)**2 + Q(i,3,3)**2) / Q(i,3,1))
      Q(i,2,4) = p_wall / (gamma - 1.d0)
      ! imaginray
      Q(i,1,1) = Q(i,3,1)
      Q(i,1,2) = Q(i,3,2) 
      Q(i,1,3) = -Q(i,3,3) 
      Q(i,1,4) = Q(i,3,4) 
    enddo
  end subroutine set_bc
end module set

