module set
  use cudafor
  use mpi
  use mod_globals, only : gamma, R, rho0, u0, p0, T0, M0
  implicit none
contains
  subroutine set_blocks_threads(myrank,nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads)
    integer, intent(in), value :: myrank, nx, ny, nz
    type(dim3), intent(out)    :: blocksE, blocksF, blocksG, blocks, threadsE, threadsF, threadsG, threads
    ! 769, 257, 257
    blocksE = dim3((nx-1)/32,(ny-2)/5,(nz-2)/1)
    blocksF = dim3((nx-2)/1,(ny-1)/128,(nz-2)/1)
    blocksG = dim3((nx-2)/1,(ny-2)/5,(nz-1)/32)
    blocks  = dim3((nx-2)/13,(ny-2)/5,(nz-2)/1)
    threadsE = dim3(32,5,1)
    threadsF = dim3(1,128,1)
    threadsG = dim3(1,5,32)
    threads  = dim3(13,5,1)
  end subroutine set_blocks_threads

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

  subroutine set_grid(myrank,nx,ny,nz,Lx,Ly,Lz,x,y,z,dx,dy,dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    integer i, j, k
    real(8) dx1, dy1, dz1
    dx1 = 0.5d0 * Lx / dble(nx-1)
    dy1 = 10.d-3 / dble(256)
    dz1 = Lz / dble(nz-1)
    x(1) = dble(myrank/2) * 0.125d0 * Lx
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo

    y(1) = 0.d0
    do j = 1, ny-1
      dy(j) = min(1.d0, max(0.05d0, dble(j)/dble(128))) * dy1
      y(j+1) = y(j) + dy(j)
    enddo

    z(1) = 0.d0
    do k = 1, nz-1
      dz(k) = dz1
      z(k+1) = z(k) + dz(k)
    enddo
  end subroutine set_grid

  subroutine set_init(myrank,nx,ny,nz,xs,ys,zs,Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: xs(nx), ys(ny), zs(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer i, j, k
    real(8) :: d = 0.2d0 * 1.d-3
    real(8) :: d1= 2.d-3
    real(8) :: eta, rho, u, v, w, T, p_wall
    ! random
    real(8) :: std, ustd, Tstd
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          eta = ys(j) / d
          call calc_Blasius(eta,d,u,v)
          if (ys(j) <= d1) then
            call random_number(std)   ! 0 <= std <= 1
            std = 2.d0 * std - 1.d0   !-1 <= std <= 1
          else
            std = 0.d0
          endif
          ustd = 0.2d0 * u * std
          u = u + ustd
          v = v + 0.5d0 * ustd
          w = 0.5d0 * ustd
          Tstd = T0 * (gamma - 1.d0) * M0**2 / u0
          T = T0 + Tstd * std
          rho = p0 / (R * T)
          Q(i,j,k,1) = rho
          Q(i,j,k,2) = Q(i,j,k,1) * u
          Q(i,j,k,3) = Q(i,j,k,1) * v
          Q(i,j,k,4) = Q(i,j,k,1) * w
          Q(i,j,k,5) = p0 / (gamma - 1.d0) + 0.5d0 * (Q(i,j,k,2)**2 + Q(i,j,k,3)**2 + Q(i,j,k,4)**2) / Q(i,j,k,1)
    enddo;enddo;enddo

    ! bottom
    Q(:,1,:,1) = Q(:,2,:,1)
    Q(:,1,:,2) = 0.d0
    Q(:,1,:,3) = 0.d0
    Q(:,1,:,4) = 0.d0
    p_wall = (gamma - 1.d0) * (Q(2,2,2,5) - 0.5d0 * (Q(2,2,2,2)**2 + Q(2,2,2,3)**2 + Q(2,2,2,4)**2) / Q(2,2,2,1))
    Q(:,1,:,5) = p_wall / (gamma - 1.d0)
  end subroutine set_init
 
  subroutine set_bc(myrank,nx,ny,nz,Jacobian,QJ)
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ(nx,ny,nz,5) ! Q / Jacobian
    integer i, j, k, l, ierr, istat(MPI_STATUS_SIZE)
    real(8) :: p_wall
    ! Riemann invariants
    real(8) :: pin, cin, vin, Rp, Rm, rhob, vb, cb, pb
    real(8) :: v0 = 0.d0
    real(8) :: c0 = sqrt(gamma * p0 / rho0)
    real(8) :: Qsend(3,ny,nz,5), Qrecv(3,ny,nz,5)
    !$cuf kernel do(2)<<<*,*>>>
    do k = 3, nz-2
      do i = 2, nx-1
        ! top
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

        ! NoSlip
        QJ(i,1,k,1) = QJ(i,2,k,1)
        QJ(i,1,k,2) = 0.d0
        QJ(i,1,k,3) = 0.d0
        QJ(i,1,k,4) = 0.d0
        p_wall = (gamma - 1.d0) * (QJ(i,2,k,5) - 0.5d0 * (QJ(i,2,k,2)**2 + QJ(i,2,k,3)**2 + QJ(i,2,k,4)**2) / QJ(i,2,k,1))
        QJ(i,1,k,5) = p_wall / (gamma - 1.d0)
    enddo;enddo

    if (myrank == 0) then
      Qsend = QJ(nx-5:nx-3,:,4:nz-3,:)
      call MPI_SENDRECV(Qsend, 15*ny*(nz-6), MPI_REAL8, 2, 0, &
                        Qrecv, 15*ny*(nz-6), MPI_REAL8, 2, 1, MPI_COMM_WORLD, istat, ierr)
      QJ(nx-2:nx,:,4:nz-3,:) = Qrecv
    elseif (myrank == 2) then
      Qsend = QJ(4:6,:,4:nz-3,:)
      call MPI_SENDRECV(Qsend, 15*ny*(nz-6), MPI_REAL8, 0, 1,&
                        Qrecv, 15*ny*(nz-6), MPI_REAL8, 0, 0, MPI_COMM_WORLD, istat, ierr)
      QJ(1:3,:,4:nz-3,:) = Qrecv
      !$cuf kernel do(3)<<<*,*>>>
      do l = 1, 5
        do k = 4, nz-3
          do j = 1, ny
            QJ(nx,j,k,l) = QJ(nx-1,j,k,l)
      enddo;enddo;enddo
    endif

    ! cyclic
    !$cuf kernel do(3)<<<*,*>>>
    do l = 1, 5
      do j = 1, ny
        do i = 1, nx
          QJ(i,j,1,l) = QJ(i,j,nz-5,l)
          QJ(i,j,2,l) = QJ(i,j,nz-4,l)
          QJ(i,j,3,l) = QJ(i,j,nz-3,l)
          QJ(i,j,nz-2,l) = QJ(i,j,4,l)
          QJ(i,j,nz-1,l) = QJ(i,j,5,l)
          QJ(i,j,nz,l)   = QJ(i,j,6,l)
    enddo;enddo;enddo
  end subroutine set_bc
end module set

