module set
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, gamma, u0, T0, &
  & beta, theta, Ms, Ms2, a1, rho0, rho2, p0, p2, u1, u2, v1, v2, u_magnitude, ux, uy
  implicit none
contains
  subroutine calc_Blasius(eta,d,u,v)
    real(8), intent(in), value  :: eta, d
    real(8), intent(out)        :: u, v
    real(8) f, df, x
    real(8) fs(45), dfs(45)
    real(8) :: nu0
    integer i
    nu0 = (1.716d-5 * ((273.2d0 + 111.d0) / (T0 + 111.d0)) * (T0 / 273.2d0)**1.5d0) / rho0 
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
      elseif (8.8d0 <= eta) then
        f = 7.07923d0
        df = 1.d0
      endif
    enddo
    u = u0 * df
    v = 0.d0
  end subroutine calc_Blasius

  subroutine set_grid(x,y,z,dx,dy)
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1)
    integer i, j, k
    real(8) :: dx1 = Lx / dble(nx-1)
    real(8) :: dy1 = Ly / dble(ny-1)
    real(8) :: dz1 = Lz / dble(nz-1)
    x(1) = 0.d0
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo

    y(1) = 0.d0
    do j = 1, ny-1
      dy(j) = max(0.25d0, 2.d0 * dble(j)/dble(ny)) * dy1
      y(j+1) = y(j) + dy(j)
    enddo

    do k = 1, nz
      z(k) = dble(k-1) * dz1
    enddo
  end subroutine set_grid

  subroutine set_init(xs,ys,zs,Q,Vin)
    real(8), intent(in)                         :: xs(nx), ys(ny), zs(nz)
    real(8), intent(out), dimension(nx,ny,nz,5) :: Q
    real(8), intent(in), dimension(ny,2)        :: Vin
    integer i, j, k
    integer :: No = int(0.3 * nx)
    real(8) :: d = 0.2d0 * 1.d-3
    real(8) :: eta, u, v, w, p_wall
    ! random
    real(8) :: ustd
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          eta = ys(j) / d
          call calc_Blasius(eta,d,u,v)
          if (ys(j) <= d) then
            call random_number(ustd)
            ustd = 2.d0 * ustd - 1.d0
          else
            ustd = 0.d0
          endif
          ustd = u * ustd
          u = u + ustd
          v = v + 0.5d0 * ustd
          w = 0.5d0 * ustd
          Q(i,j,k,1) = rho0
          Q(i,j,k,2) = Q(i,j,k,1) * u
          Q(i,j,k,3) = Q(i,j,k,1) * v
          Q(i,j,k,4) = Q(i,j,k,1) * w
          Q(i,j,k,5) = p0 / (gamma - 1.d0) + 0.5d0 * (Q(i,j,k,2)**2 + Q(i,j,k,3)**2 + Q(i,j,k,4)**2) / Q(i,j,k,1)
    enddo;enddo;enddo

    ! oblique shock
    do k = 2, nz-1
      do i = No+1, nx
        Q(i,ny,k,1) = rho2
        Q(i,ny,k,2) = rho2 * ux
        Q(i,ny,k,3) = rho2 * uy
        Q(i,ny,k,4) = 0.d0
        Q(i,ny,k,5) = p2 / (gamma - 1.d0) + 0.5d0 * rho2 * (ux**2 + uy**2)
    enddo;enddo

    ! bottom
    Q(:,1,:,1) = Q(:,2,:,1)
    Q(:,1,:,2) = 0.d0
    Q(:,1,:,3) = 0.d0
    Q(:,1,:,4) = 0.d0
    p_wall = (gamma - 1.d0) * (Q(2,2,2,5) - 0.5d0 * (Q(2,2,2,2)**2 + Q(2,2,2,3)**2 + Q(2,2,2,4)**2) / Q(2,2,2,1))
    Q(:,1,:,5) = p_wall / (gamma - 1.d0)
  end subroutine set_init
  
  subroutine set_bc(id_accuracy,Jacobian,QJ)
    integer(kind=2), intent(in), value  :: id_accuracy
    real(8), intent(in), device         :: Jacobian(nx,ny)
    real(8), intent(inout), device      :: QJ(nx,ny,nz,5)
    integer i, j, k, l
    integer :: No = int(0.3 * nx)
    integer :: Nre = int(0.25 * nx)
    real(8) :: p_wall
    ! Riemann invariants
    real(8) :: pin, cin, vin, Rp, Rm, rhob, vb, cb, pb
    real(8) :: v0 = 0.d0
    real(8) :: c0 = sqrt(gamma * p0 / rho0)
    !$cuf kernel do<<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        ! inlet
        QJ(1,j,k,1) = QJ(Nre,j,k,1)
        QJ(1,j,k,2) = QJ(Nre,j,k,2)
        QJ(1,j,k,3) = QJ(Nre,j,k,3)
        QJ(1,j,k,4) = QJ(Nre,j,k,4)
        QJ(1,j,k,5) = QJ(Nre,j,k,5)
        ! outlet
        QJ(nx,j,k,1) = QJ(nx-1,j,k,1)
        QJ(nx,j,k,2) = QJ(nx-1,j,k,2)
        QJ(nx,j,k,3) = QJ(nx-1,j,k,3)
        QJ(nx,j,k,4) = QJ(nx-1,j,k,4)
        QJ(nx,j,k,5) = QJ(nx-1,j,k,5)
    enddo;enddo

    ! Riemann boundary condition
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do i = 1, No
        ! Riemann invariants
        pin = (gamma - 1.d0) * (QJ(i,ny-1,k,5) - 0.5d0 * (QJ(i,ny-1,k,2)**2 + QJ(i,ny-1,k,3)**2 + QJ(i,ny-1,k,4)**2) / QJ(i,ny-1,k,1)) &
        & * Jacobian(i,ny-1)
        cin = sqrt(gamma * pin / (QJ(i,ny-1,k,1) * Jacobian(i,ny-1)))
        vin = QJ(i,ny-1,k,3) / QJ(i,ny-1,k,1)
        Rp = vin + 2.d0 * cin / (gamma - 1.d0)
        Rm = v0  - 2.d0 * c0  / (gamma - 1.d0)
        vb = v0 + (0.5d0 * (Rp + Rm) - v0)

        rhob = QJ(i,ny-1,k,1) * Jacobian(i,ny-1)
        QJ(i,ny,k,1) = rhob / Jacobian(i,ny)
        QJ(i,ny,k,2) = QJ(i,ny-1,k,1) * u0 
        QJ(i,ny,k,3) = QJ(i,ny-1,k,1) * vb
        QJ(i,ny,k,4) = 0.d0
        cb = 0.25d0 * (gamma - 1.d0) * (Rp - Rm)
        pb = (rhob * cb**2) / gamma
        QJ(i,ny,k,5) = (pb / (gamma - 1.d0)) / Jacobian(i,ny)  + 0.5d0 * (QJ(i,ny,k,2)**2 + QJ(i,ny,k,3)**2 + QJ(i,ny,k,4)**2) / QJ(i,ny,k,1)
    enddo;enddo

    ! oblique shock
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do i = No+1, nx
        QJ(i,ny,k,1) = rho2 / Jacobian(i,ny)
        QJ(i,ny,k,2) = rho2 * ux / Jacobian(i,ny)
        QJ(i,ny,k,3) = rho2 * uy / Jacobian(i,ny)
        QJ(i,ny,k,4) = 0.d0
        QJ(i,ny,k,5) = (p2 / (gamma - 1.d0) + 0.5d0 * rho2 * (ux**2 + uy**2)) / Jacobian(i,ny)
    enddo;enddo

    ! bottom
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do i = 1, nx
        ! NoSlip
        QJ(i,1,k,1) = QJ(i,2,k,1)
        QJ(i,1,k,2) = 0.d0
        QJ(i,1,k,3) = 0.d0
        QJ(i,1,k,4) = 0.d0
        p_wall = (gamma - 1.d0) * (QJ(i,2,k,5) - 0.5d0 * (QJ(i,2,k,2)**2 + QJ(i,2,k,3)**2 + QJ(i,2,k,4)**2) / QJ(i,2,k,1))
        QJ(i,1,k,5) = p_wall / (gamma - 1.d0)
    enddo;enddo

    ! cyclic
    !$cuf kernel do(3)<<<*,*>>>
    do l = 1, 5
      do j = 1, ny
        do i = 1, nx
          QJ(i,j,1,l) = QJ(i,j,nz-3,l)
          QJ(i,j,2,l) = QJ(i,j,nz-2,l)
          QJ(i,j,nz-1,l) = QJ(i,j,3,l)
          QJ(i,j,nz,l) = QJ(i,j,4,l)
    enddo;enddo;enddo
  end subroutine set_bc

  subroutine set_bc_mut(mut)
    real(8), intent(inout), device :: mut(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2) <<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        ! inlet
        mut(1,j,k) = mut(2,j,k)
        ! outlet
        mut(nx,j,k) = mut(nx-1,j,k)
    enddo;enddo

    !$cuf kernel do(2) <<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        ! wall
        mut(i,1,k) = 0.d0
        ! top
        mut(i,ny,k) = mut(i,ny-1,k)
    enddo;enddo

    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        ! span
        mut(i,j,1) = mut(i,j,2)
        mut(i,j,nz) = mut(i,j,nz-1)
    enddo;enddo
  end subroutine set_bc_mut
end module set

