module set
  use cudafor
  use mpi
  use mod_globals, only : Lx1, Ly1, Lz1, Lx2, Ly2, Lz2, gamma, R, u0, T0, &
  & beta, theta, M0, Ms, Ms2, a1, rho0, rho2, p0, p2, u1, u2, v1, v2, u_magnitude, ux, uy
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

  subroutine set_grid(myrank,nx,ny,nz,x,y,z,dx,dy,dz)
    integer, intent(in)   :: myrank, nx, ny, nz
    real(8), intent(out)  :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    integer i, j, k, ierr, status(MPI_STATUS_SIZE)
    real(8) dx1, dy1, dz1, dx2, dy2, dz2
    dy1 = 12d-3 / dble(256)!Ly1 / dble(256)
    if (myrank == 0) then
      dx1 = Lx1 / dble(nx-1)
      x(1) = 0.d0
      do i = 1, nx-1
        dx(i) = dx1
        x(i+1) = x(i) + dx(i)
      enddo

      call MPI_RECV(y(1), 1, MPI_REAL8, 2, 0, MPI_COMM_WORLD, status, ierr)
      do j = 1, ny-1
        dy(j) = dy1
        y(j+1) = y(j) + dy(j)
      enddo

      dz1 = Lz1 / dble(nz-1)
      do k = 1, nz
        dz(k) = dz1
        z(k) = dble(k-1) * dz1
      enddo
    elseif (myrank ==2) then
      dx2 = Lx2 / dble(nx-1)
      x(1) = 0.d0
      do i = 1, nx-1
        dx(i) = dx2
        x(i+1) = x(i) + dx(i)
      enddo

      y(1) = 0.d0
      do j = 1, ny-3
        dy(j) = min(1.d0, max(0.5d0, dble(j) / dble(ny-3))) * dy1
        y(j+1) = y(j) + dy(j)
        !print *, dy1, dy(j)
      enddo
      dy(ny-2) = dy1
      y(ny-1) = y(ny-2) + dy1
      dy(ny-1) = dy1
      y(ny) = y(ny-1) + dy1
      call MPI_SEND(y(ny-1), 1, MPI_REAL8, 0, 0, MPI_COMM_WORLD, ierr)
      
      dz2 = Lz2 / dble(nz-1)
      do k = 1, nz
        dz(k) = dz2
        z(k) = dble(k-1) * dz2
      enddo
    endif
  end subroutine set_grid

  subroutine set_init(myrank,nx,ny,nz,xs,ys,zs,Q)
    integer, intent(in)                         :: myrank, nx, ny, nz
    real(8), intent(in)                         :: xs(nx), ys(ny), zs(nz)
    real(8), intent(out), dimension(nx,ny,nz,5) :: Q
    integer i, j, k, No
    real(8) :: d = 0.2d0 * 1.d-3
    real(8) :: d1= 2.d-3
    real(8) :: eta, rho, u, v, w, T, p_wall
    ! random
    real(8) :: std, ustd, Tstd
    if (myrank == 0) then
      rho = p0 / (R * T0)
      Q(:,:,:,1) = rho
      Q(:,:,:,2) = rho * u0
      Q(:,:,:,3) = 0.d0
      Q(:,:,:,4) = 0.d0
      Q(:,:,:,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho * u0**2
    elseif (myrank == 2) then
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            eta = ys(j) / d
            call calc_Blasius(eta,d,u,v)
            if (ys(j) <= d1) then
              call random_number(std)
              std = 2.d0 * std - 1.d0
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
    endif
  end subroutine set_init
  
  subroutine set_blocks_threads(myrank,nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads)
    integer, intent(in)       :: myrank, nx, ny, nz
    type(dim3), intent(inout) :: blocksE, blocksF, blocksG, blocks, threadsE, threadsF, threadsG, threads
    integer :: accuracy = 2
    if (myrank == 0) then
      blocksE  = dim3((nx-accuracy+1)/32,(ny-accuracy)/1,(nz-accuracy)/3)
      blocksF  = dim3((nx-accuracy)/5,(ny-accuracy+1)/32,(nz-accuracy)/1)
      blocksG  = dim3((nx-accuracy)/5,(ny-accuracy)/1,(nz-accuracy+1)/32)
      blocks   = dim3((nx-accuracy)/5,(ny-accuracy)/1,(nz-accuracy)/9)
      threadsE = dim3(32,1,3)
      threadsF = dim3(5,32,1)
      threadsG = dim3(5,1,32)
      threads  = dim3(5,1,9)
    elseif (myrank == 2) then
      blocksE  = dim3((nx-accuracy+1)/32,(ny-accuracy)/3,(nz-accuracy)/1)
      blocksF  = dim3((nx-accuracy)/7,(ny-accuracy+1)/32,(nz-accuracy)/1)
      blocksG  = dim3((nx-accuracy)/7,(ny-accuracy)/1,(nz-accuracy+1)/32)
      blocks   = dim3((nx-accuracy)/7,(ny-accuracy)/3,(nz-accuracy)/1)
      threadsE = dim3(32,3,1)
      threadsF = dim3(7,32,1)
      threadsG = dim3(7,1,32)
      threads  = dim3(7,3,1)
    endif
  end subroutine set_blocks_threads

  subroutine GPU_SENDRECV(myrank,nx,ny,nz,Q,Qd)
    use mod_globals, only : nx1, ny1, nz1, nx2, ny2, nz2
    integer, intent(in)          :: myrank, nx, ny, nz
    real(8), intent(in), device  :: Q(nx,ny,nz,5) ! without Jacobian
    real(8), intent(out), device :: Qd(nx,2,nz,5)
    real(8), allocatable, pinned :: Qp1(:,:,:,:), Qp2(:,:,:,:)
    real(8), device :: Qi1(nx2,2,nz2,5), Qi2(nx1,2,nz1,5)
    integer i, k, ierr, status(MPI_STATUS_SIZE)
    if (myrank == 0) then
      allocate(Qp1(nx2,2,nz2,5),Qp2(nx1,2,nz1,5))
      call MPI_RECV(Qp2, 10*nx1*nz1, MPI_REAL8, 2, 0, MPI_COMM_WORLD, status, ierr)
      Qd = Qp2
      !$cuf kernel do(2)<<<*,*>>>
      do k = 1, nz1
        do i = 1, nx1
          ! i odd, k odd
          Qi1(2*i-1,1,2*k-1,:) = 0.5d0 * (Q(i,3,k,:) + Qd(i,2,k,:))
          Qi1(2*i-1,2,2*k-1,:) = Q(i,3,k,:)
      enddo;enddo

      !$cuf kernel do(2)<<<*,*>>>
      do k = 1, nz1-1
        do i = 1, nx1
          ! i odd, k even
          Qi1(2*i-1,1,2*k,:) = 0.25d0 * (Q(i,3,k,:) + Q(i,3,k+1,:) + Qd(i,2,k,:) + Qd(i,2,k+1,:))
          Qi1(2*i-1,2,2*k,:) = 0.5d0 * (Q(i,3,k,:) + Q(i,3,k+1,:))
      enddo;enddo

      !$cuf kernel do(2)<<<*,*>>>
      do k = 1, nz1
        do i = 1, nx1-1
          ! i even, k odd
          Qi1(2*i,1,2*k-1,:) = 0.25d0 * (Q(i,3,k,:) + Q(i+1,3,k,:) + Qd(i,2,k,:) + Qd(i,2,k,:))
          Qi1(2*i,2,2*k-1,:) = 0.5d0 * (Q(i,3,k,:) + Q(i+1,3,k,:))
      enddo;enddo

      !$cuf kernel do(2)<<<*,*>>>
      do k = 1, nz1-1
        do i = 1, nx1-1
          ! i even, k even
          Qi1(2*i,1,2*k,:) = 0.125d0 * (Q(i,3,k,:)  + Q(i,3,k+1,:)  + Q(i+1,3,k,:)  + Q(i+1,3,k+1,:) + &
                                        Qd(i,2,k,:) + Qd(i,2,k+1,:) + Qd(i+1,2,k,:) + Qd(i+1,2,k+1,:))
          Qi1(2*i,2,2*k,:) = 0.25d0 * (Q(i,3,k,:) + Q(i,3,k+1,:) + Q(i+1,3,k,:) + Q(i+1,3,k+1,:))
      enddo;enddo

      Qp1 = Qi1
      call MPI_SEND(Qp1, 10*nx2*nz2, MPI_REAL8, 2, 0, MPI_COMM_WORLD, ierr)
    else
      allocate(Qp1(nx1,2,nz1,5),Qp2(nx2,2,nz2,5))
      !$cuf kernel do(2)<<<*,*>>>
      do k = 1, nz1
        do i = 1, nx1
          Qi2(i,1,k,:) = Q(2*i-1,ny-4,2*k-1,:)
          Qi2(i,2,k,:) = Q(2*i-1,ny-2,2*k-1,:)
      enddo;enddo
      Qp1 = Qi2
      call MPI_SEND(Qp1, 10*nx1*nz1, MPI_REAL8, 0, 0, MPI_COMM_WORLD, ierr)
      call MPI_RECV(Qp2, 10*nx2*nz2, MPI_REAL8, 0, 0, MPI_COMM_WORLD, status, ierr)
      Qd = Qp2
    endif
    deallocate(Qp1,Qp2)
  end subroutine

  subroutine set_bc1(nx,ny,nz,Jacobian,Qre,QJ)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(in), device    :: Qre(2,ny,nz,5)
    real(8), intent(inout), device :: QJ(nx,ny,nz,5)
    integer i, j, k, l, No, Nre
    real(8) :: p_wall
    ! Riemann invariants
    real(8) :: pin, cin, vin, Rp, Rm, rhob, vb, cb, pb
    real(8) :: v0 = 0.d0
    real(8) :: c0 = sqrt(gamma * p0 / rho0)
    real(8), device :: Qd(nx,2,nz,5)
    real(8), device :: Q(nx,ny,nz,5)
    Nre = int(0.3 * nx)
    !$cuf kernel do<<<*,*>>>
    do l = 1, 5
      do k = 3, nz-2
        do j = 3, ny-1
          ! inlet 
          QJ(1,j,k,l)    = Qre(1,j,k,l)
          QJ(2,j,k,l)    = Qre(2,j,k,l)
          ! outlet
          QJ(nx-1,j,k,l) = QJ(nx-3,j,k,l)
          QJ(nx,j,k,l)   = QJ(nx-2,j,k,l)
    enddo;enddo;enddo

    ! Riemann boundary condition
    !$cuf kernel do(2)<<<*,*>>>
    do k = 3, nz-2
      do i = 1, nx
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

    ! cyclic
    !$cuf kernel do(3)<<<*,*>>>
    do l = 1, 5
      do j = 3, ny
        do i = 1, nx
          QJ(i,j,1,l) = QJ(i,j,nz-3,l)
          QJ(i,j,2,l) = QJ(i,j,nz-2,l)
          QJ(i,j,nz-1,l) = QJ(i,j,3,l)
          QJ(i,j,nz,l) = QJ(i,j,4,l)
    enddo;enddo;enddo

    ! parallel
    !$cuf kernel do(4)<<<*,*>>>
    do l = 1, 5
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            Q(i,j,k,l) = QJ(i,j,k,l) * Jacobian(i,j)
    enddo;enddo;enddo;enddo
    call GPU_SENDRECV(0,nx,ny,nz,Q,Qd)
    !$cuf kernel do(3)<<<*,*>>>
    do l = 1, 5
      do k = 1, nz
        do i = 1, nx
          QJ(i,1,k,l) = Qd(i,1,k,l) / Jacobian(i,1)
          QJ(i,2,k,l) = Qd(i,2,k,l) / Jacobian(i,2)
    enddo;enddo;enddo
  end subroutine set_bc1

  subroutine set_bc2(nx,ny,nz,Jacobian,Qre,QJ)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(in), device    :: Qre(2,ny,nz,5)
    real(8), intent(inout), device :: QJ(nx,ny,nz,5)
    integer i, j, k, l, Nre
    real(8) :: p_wall
    real(8), device :: Qd(nx,2,nz,5)
    real(8), device :: Q(nx,ny,nz,5)
    Nre = int(0.3 * nx)
    !$cuf kernel do<<<*,*>>>
    do l = 1, 5
      do k = 3, nz-2
        do j = 2, ny-2
          ! inlet 
          QJ(1,j,k,l)    = Qre(1,j,k,l)
          QJ(2,j,k,l)    = Qre(2,j,k,l)
          ! outlet
          QJ(nx-1,j,k,l) = QJ(nx-3,j,k,l)
          QJ(nx,j,k,l)   = QJ(nx-2,j,k,l)
    enddo;enddo;enddo

    ! bottom
    !$cuf kernel do(2)<<<*,*>>>
    do k = 3, nz-2
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
      do j = 1, ny-2
        do i = 1, nx
          QJ(i,j,1,l) = QJ(i,j,nz-3,l)
          QJ(i,j,2,l) = QJ(i,j,nz-2,l)
          QJ(i,j,nz-1,l) = QJ(i,j,3,l)
          QJ(i,j,nz,l) = QJ(i,j,4,l)
    enddo;enddo;enddo

    ! parallel
    !$cuf kernel do(4)<<<*,*>>>
    do l = 1, 5
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            Q(i,j,k,l) = QJ(i,j,k,l) * Jacobian(i,j)
    enddo;enddo;enddo;enddo
    call GPU_SENDRECV(2,nx,ny,nz,Q,Qd)
    !$cuf kernel do(3)<<<*,*>>>
    do l = 1, 5
      do k = 1, nz
        do i = 1, nx
          QJ(i,ny-1,k,l) = Qd(i,1,k,l) / Jacobian(i,ny-1)
          QJ(i,ny,k,l)   = Qd(i,2,k,l) / Jacobian(i,ny)
    enddo;enddo;enddo
  end subroutine set_bc2

  subroutine set_bc_mut(nx,ny,nz,mut)
    integer, intent(in), value      :: nx, ny, nz
    real(8), intent(inout), device  :: mut(nx,ny,nz)
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

