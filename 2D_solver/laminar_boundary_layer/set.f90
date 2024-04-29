module set
  use cudafor
  use mod_globals, only : nx, ny, dx, dy, gamma, u0, R, p0, rho0
  implicit none
contains
  subroutine calc_Blasius(eta,d,u,v)
    real(8), intent(in), value :: eta, d
    real(8), intent(out) :: u, v
    real(8) f, df, x
    real(8) fs(45), dfs(45)
    real(8) :: nu0 = 3.8206d-5
    integer i
    fs(:) = (/0.d0, 0.00664d0, 0.02656d0, 0.05974d0, 0.10611d0, 0.16557d0, 0.23795d0, &
    & 0.32298d0, 0.42032d0, 0.52952d0, 0.65003d0, 0.78120d0, 0.92230d0, 1.07252d0, &
    & 1.23099d0, 1.39682d0, 1.56911d0, 1.74696d0, 1.92954d0, 2.11605d0, 2.30576d0, &
    & 2.49806d0, 2.69238d0, 2.88826d0, 3.08534d0, 3.28329d0, 3.48189d0, 3.68094d0, &
    & 3.88031d0, 4.07990d0, 4.27964d0, 4.47948d0, 4.67938d0, 4.87931d0, 5.07928d0, &
    & 5.27926d0, 5.47925d0, 5.67924d0, 5.87924d0, 6.07923d0, 6.27923d0, 6.47923d0, &
    & 6.67923d0, 6.87923d0, 7.07923d0/)
    dfs(:) = (/0.d0, 0.06641d0, 0.13277d0, 0.19894d0, 0.26471d0, 0.32979d0, 0.39378d0, &
    & 0.45627d0, 0.51676d0, 0.57477d0, 0.62977d0, 0.68132d0, 0.72899d0, 0.77246d0, &
    & 0.81152d0, 0.84605d0, 0.87609d0, 0.90177d0, 0.92333d0, 0.94112d0, 0.95552d0, &
    & 0.96696d0, 0.97587d0, 0.98269d0, 0.98779d0, 0.99155d0, 0.99425d0, 0.99616d0, &
    & 0.99748d0, 0.99838d0, 0.99898d0, 0.99937d0, 0.99961d0, 0.99977d0, 0.99987d0, &
    & 0.99992d0, 0.99996d0, 0.99998d0, 0.99999d0, 1.00000d0, 1.00000d0, 1.00000d0, &
    & 1.00000d0, 1.00000d0, 1.00000d0/)
    do i = 1, 44
      x = 0.2d0 * dble(i-1) 
      if (x <= eta .and. eta <= x + 0.2d0) then
        f = fs(i) + 5.d0 * (fs(i+1) - fs(i)) * (eta - x)
        df = dfs(i) + 5.d0 * (dfs(i+1) - dfs(i)) * (eta - x)
      elseif (8.8d0 < eta) then
        f = 7.07923d0
        df = 1.d0
      endif
    enddo
    u = u0 * df
    v = 0.5d0 * (nu0 / d) * (eta * df - f)
    write(*,*) eta, f, df
  end subroutine calc_Blasius

  subroutine set_grid(x,y,z)
    real(8), intent(out), dimension(nx,ny) :: x, y, z
    integer i, j
    do j = 1, ny
      do i = 1, nx
        x(i,j) = dble(i-1) * dx
        y(i,j) = dble(j-1) * dy
        z(i,j) = 0.d0
      enddo
    enddo
  end subroutine set_grid

  subroutine set_init(Q,Vin)
    real(8), intent(out), dimension(nx,ny,4) :: Q
    real(8), intent(in), dimension(ny,2) :: Vin
    real(8) y, eta, u, v, p_wall
    real(8) :: d = 1.d-3
    real(8) :: nu0 = 3.8206d-5
    integer i, j
    Q(:,:,1) = rho0
    ! boundary layer
    do j = 2, ny
      y = dble(j-2) * dy
      eta = y / d
      call calc_Blasius(eta,d,u,v)
      Q(:,j,2) = rho0 * u
      Q(:,j,3) = rho0 * v
      Q(:,j,4) = p0 / (gamma - 1.d0) + 0.5d0 * (Q(:,j,2)**2 + Q(:,j,3)**2) / Q(:,j,1) 
    enddo

    ! noSlip
    do i = 1, nx
      ! wall
      Q(i,1,1) = Q(i,2,1)
      Q(i,1,2) = 0.d0
      Q(i,1,3) = 0.d0
      p_wall = (gamma - 1.d0) * (Q(i,2,4) - 0.5d0 * (Q(i,2,2)**2 + Q(i,2,3)**2) / Q(i,2,1))
      Q(i,1,4) = p_wall / (gamma - 1.d0)
    enddo
  end subroutine set_init

  subroutine set_bc(Q,Vin)
    real(8), intent(inout), device :: Q(nx,ny,4)
    real(8), intent(in), device :: Vin(ny,2)
    integer i, j, k
    real(8) p_wall
    !$cuf kernel do <<<*,*>>>
    do k = 1, 4
      do j = 2, ny-1
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
      Q(i,1,1) = Q(i,2,1)
      Q(i,1,2) = 0.d0
      Q(i,1,3) = 0.d0
      p_wall = (gamma - 1.d0) * (Q(i,2,4) - 0.5d0 * (Q(i,2,2)**2 + Q(i,2,3)**2) / Q(i,2,1))
      Q(i,1,4) = p_wall / (gamma - 1.d0)
    enddo
  end subroutine set_bc
end module set

