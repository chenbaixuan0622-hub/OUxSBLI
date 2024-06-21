module calc_time_dev
  use cudafor
  use mpi
  use nvtx
  use mod_globals, only : accuracy, id_hybrid, id_muscl, id_visc, id_turbulence, nt, np, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG
  use calc_physical_quantities
  use calc_steps
  use calc_hybrid
  use calc_visc
  use calc_les
  use set
  use print
  implicit none
  interface calc_EFG
    module procedure calc_EFG_basic, calc_EFG_hybrid
  end interface
  interface RungeKutta
    module procedure RungeKutta_3rd, RungeKutta_4th, RungeKutta_10th
  end interface
contains
  subroutine calc_EFG_basic(id_hybrid,nx,ny,nz,dx,dy,Jacobian,QJ,mut,E,F,G)
    use calc_flux
    integer(kind=2), intent(in)                         :: id_hybrid
    integer, intent(in), value                          :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device        :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1), device        :: dy ! 1 / dy
    real(8), intent(in), dimension(nx,ny), device       :: Jacobian
    real(8), intent(in), dimension(nx,ny,nz,5), device  :: QJ ! Q / Jacobian
    real(8), intent(inout), dimension(nx,ny,nz), device :: mut
    real(8), intent(out), device                        :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device                        :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device                        :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p, T
    integer stat, itr
    call calc_quantities(nx,ny,nz,Jacobian,QJ,rho,u,v,w,p,T)

    call calc_E<<<blocksE,threadsE,1>>>(id_muscl,nx,ny,nz,rho,u,v,w,p,E)
    call calc_F<<<blocksF,threadsF,2>>>(id_muscl,nx,ny,nz,rho,u,v,w,p,F)
    call calc_G<<<blocksG,threadsG,3>>>(id_muscl,nx,ny,nz,rho,u,v,w,p,G)

    if (id_turbulence /= 0) then
      call calc_mut<<<blocks,threads,4>>>(nx,ny,nz,rho,u,v,w,mut)
      stat = cudaDeviceSynchronize()
      call set_bc_mut(nx,ny,nz,mut)
    endif

    stat = cudaDeviceSynchronize()
    if (id_visc == 1 .or. id_turbulence /= 0) then
      call calc_Ev<<<blocksE,threadsE,1>>>(nx,ny,nz,dx,dy,rho,u,v,w,T,p,mut,E)
      call calc_Fv<<<blocksF,threadsF,2>>>(nx,ny,nz,dy,dx,rho,u,v,w,T,p,mut,F)
      call calc_Gv<<<blocksG,threadsG,3>>>(nx,ny,nz,dx,dy,rho,u,v,w,T,p,mut,G)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_basic
  
  subroutine calc_EFG_hybrid(id_hybrid,nx,ny,nz,dx,dy,Jacobian,QJ,mut,E,F,G)
    use calc_flux_hybrid
    integer(kind=4), intent(in)                         :: id_hybrid
    integer, intent(in), value                          :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device        :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1), device        :: dy ! 1 / dy
    real(8), intent(in), dimension(nx,ny), device       :: Jacobian
    real(8), intent(in), dimension(nx,ny,nz,5), device  :: QJ ! Q / Jacobian
    real(8), intent(inout), dimension(nx,ny,nz), device :: mut
    real(8), intent(out), device                        :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device                        :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device                        :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx,ny,nz), device                                :: rho, u, v, w, p, T, fd
    !real(8), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E_upwind
    !real(8), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F_upwind
    !real(8), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G_upwind
    integer stat

    call calc_quantities(nx,ny,nz,Jacobian,QJ,rho,u,v,w,p,T)

    !call calc_E<<<blocksE,threadsE>>>(id_muscl1,nx,ny,nz,rho,u,v,w,p,E)
    !call calc_F<<<blocksF,threadsF>>>(id_muscl1,nx,ny,nz,rho,u,v,w,p,F)
    !call calc_G<<<blocksG,threadsG>>>(id_muscl1,nx,ny,nz,rho,u,v,w,p,G)
    !call calc_E<<<blocksE,threadsE>>>(id_muscl,nx,ny,nz,rho,u,v,w,p,E_upwind)
    !call calc_F<<<blocksF,threadsF>>>(id_muscl,nx,ny,nz,rho,u,v,w,p,F_upwind)
    !call calc_G<<<blocksG,threadsG>>>(id_muscl,nx,ny,nz,rho,u,v,w,p,G_upwind)

    call calc_Ducros<<<blocks,threads>>>(nx,ny,nz,dx,dy,u,v,w,rho,p,fd)
    call calc_E_hybrid<<<blocksE,threadsE,1>>>(nx,ny,nz,rho,u,v,w,p,fd,E)
    call calc_F_hybrid<<<blocksF,threadsF,2>>>(nx,ny,nz,rho,u,v,w,p,fd,F)
    call calc_G_hybrid<<<blocksG,threadsG,3>>>(nx,ny,nz,rho,u,v,w,p,fd,G)
    !call calc_E_hybrid<<<blocksE,threadsE>>>(nx,ny,nz,u,v,w,fd,E_upwind,E)
    !call calc_F_hybrid<<<blocksF,threadsF>>>(nx,ny,nz,u,v,w,fd,F_upwind,F)
    !call calc_G_hybrid<<<blocksG,threadsG>>>(nx,ny,nz,u,v,w,fd,G_upwind,G)
    !print *, trim(cudaGetErrorString(cudaGetLastError()))

    if (id_turbulence /= 0) then
      call calc_mut<<<blocks,threads,4>>>(nx,ny,nz,rho,u,v,w,mut)
      stat = cudaDeviceSynchronize()
      call set_bc_mut(nx,ny,nz,mut)
    endif

    stat = cudaDeviceSynchronize()
    if (id_visc == 1 .or. id_turbulence /= 0) then
      call calc_Ev<<<blocksE,threadsE,1>>>(nx,ny,nz,dx,dy,rho,u,v,w,T,p,mut,E)
      call calc_Fv<<<blocksF,threadsF,2>>>(nx,ny,nz,dy,dx,rho,u,v,w,T,p,mut,F)
      call calc_Gv<<<blocksG,threadsG,3>>>(nx,ny,nz,dx,dy,rho,u,v,w,T,p,mut,G)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_hybrid

  subroutine RungeKutta_3rd(id_RungeKutta,myrank,nx,ny,nz,x,dx_cpu,xix_cpu,y,dy_cpu,etay_cpu,z,Jacobian_cpu,Q)
    integer(kind=2), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1), xix_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1), etay_cpu(ny-1)
    real(8), intent(in)         :: z(nz), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    integer i, j, k, t1, t2, itr, ilen, ierr, stat, request, status(MPI_STATUS_SIZE)
    real(8), dimension(nx,ny,nz)    :: mut_cpu
    real(8), dimension(nx,ny,nz,3)  :: Vmean_cpu
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)          :: prop
    real(8), allocatable, device  :: QJ(:,:,:,:), QJ2(:,:,:,:), QJ3(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device  :: dx(:), xix(:), dy(:), etay(:), Jacobian(:,:), mut(:,:,:), Vmean(:,:,:,:), rhomean(:,:,:), Tmean(:,:,:)
    ! for plot
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    ! check GPU
    stat = cudaSetDevice(0)
    stat = cudaGetDeviceProperties(prop,0)
    ilen = verify(prop%name, ' ', .true.)
    print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", 0, ") is available"

    if (myrank == 0) then
      allocate(QJ(nx,ny,nz,5),QJ2(nx,ny,nz,5),QJ3(nx,ny,nz,5),E(nx-1,ny-2,nz-2,5),F(nx-2,ny-1,nz-2,5),G(nx-2,ny-2,nz-1,5))
      allocate(dx(nx-1),xix(nx-1),dy(ny-1),etay(ny-1),Jacobian(nx,ny),mut(nx,ny,nz),Vmean(nx,ny,nz,3),rhomean(nx,ny,nz),Tmean(nx,ny,nz))

      ! set Q / Jacobian
      do j = 1, ny
        do i = 1, nx
          Q(i,j,:,:) = Q(i,j,:,:) / Jacobian_cpu(i,j)
      enddo;enddo
  
      ! print initial condition
      if (id_turbulence == 0) then
        call print_vtk(0,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),ke0,entropy0)
      else
        call print_vtk(0,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),ke0,entropy0,real(mut_cpu))
      endif

      ! copy on GPU
      QJ = Q
      mut = 0.d0
      xix = xix_cpu
      etay = etay_cpu
      dx = dx_cpu
      dy = dy_cpu
      Jacobian = Jacobian_cpu
    endif

    ! share necessary data
    call MPI_BCAST(ke0, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(entropy0, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    do t2 = 1, np
      if (myrank == 0) then
        do t1 = 1, nt
          call nvtxStartRange("calc 1step",1)
          call nvtxStartRange("calc flux",2)
          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJ,mut,E,F,G)
          call nvtxEndRange
          call nvtxStartRange("calc time dev",3)
          call calc_step(nx,ny,nz,1.d0,0.d0,dx,dy,E,F,G,QJ,QJ2)
          call nvtxEndRange
          call set_bc(nx,ny,nz,Jacobian,QJ2)

          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJ2,mut,E,F,G)
          call calc_step2(nx,ny,nz,0.75d0,0.25d0,0.25d0,1.d0,dx,dy,E,F,G,QJ,QJ2,QJ3)
          call set_bc(nx,ny,nz,Jacobian,QJ3)

          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJ3,mut,E,F,G)
          call calc_step3(nx,ny,nz,dx,dy,E,F,G,QJ3,QJ)
          call set_bc(nx,ny,nz,Jacobian,QJ)
          call calc_mean(t1+(t2-1)*nt,nx,ny,nz,Jacobian,QJ,Vmean,rhomean,Tmean)
          call nvtxEndRange
        enddo
      endif

      ! send and recv device arrays
      if (myrank == 0) then
        Q = QJ
        call nvtxStartRange("MPI_SEND",4)
        if (id_turbulence /= 0) then
          mut_cpu = mut
          call MPI_SEND(mut_cpu, nx*ny*nz, MPI_REAL8, 1, 1, MPI_COMM_WORLD, ierr)
        endif
        call MPI_SEND(Q, nx*ny*nz*5, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
        call nvtxEndRange
      elseif (myrank == 1) then
        call nvtxStartRange("MPI_RECV",5)
        call MPI_RECV(Q, nx*ny*nz*5, MPI_REAL8, 0, 0, MPI_COMM_WORLD, status, ierr)
        call nvtxEndRange
        call nvtxStartRange("print",6)
        if (id_turbulence == 0) then
          call print_vtk(t2,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),ke0,entropy0)
        else
          call MPI_RECV(mut_cpu, nx*ny*nz, MPI_REAL8, 0, 1, MPI_COMM_WORLD, status, ierr)
          call print_vtk(t2,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),ke0,entropy0,real(mut_cpu)) 
        endif
        call nvtxEndRange
      endif
    enddo
    
    if (myrank == 0) then
      deallocate(QJ,QJ2,QJ3,E,F,G,dx,xix,dy,etay,Jacobian,mut,Vmean,rhomean,Tmean)
    endif
  end subroutine RungeKutta_3rd

  subroutine RungeKutta_4th(id_RungeKutta,myrank,nx,ny,nz,x,dx_cpu,xix_cpu,y,dy_cpu,etay_cpu,z,Jacobian_cpu,Q)
    integer(kind=4), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1), xix_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1), etay_cpu(ny-1)
    real(8), intent(in)         :: z(nz), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    integer i, j, k, t1, t2, itr, ierr, ilen, stat, status(MPI_STATUS_SIZE)
    real(8), dimension(nx,ny,nz)    :: mut_cpu
    real(8), dimension(nx,ny,nz,3)  :: Vmean_cpu
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)          :: prop
    real(8), allocatable, device  :: QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device  :: dx(:), xix(:), dy(:), etay(:), Jacobian(:,:), mut(:,:,:), Vmean(:,:,:,:), rhomean(:,:,:), Tmean(:,:,:)
    ! for plot
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    ! check GPU
    stat = cudaSetDevice(0)
    stat = cudaGetDeviceProperties(prop,0)
    ilen = verify(prop%name, ' ', .true.)

    if (myrank == 0) then
      allocate(QJ(nx,ny,nz,5),QJs(nx,ny,nz,5),Rs(nx,ny,nz,5),E(nx-1,ny-2,nz-2,5),F(nx-2,ny-1,nz-2,5),G(nx-2,ny-2,nz-1,5))
      allocate(dx(nx-1),xix(nx-1),dy(ny-1),etay(ny-1),Jacobian(nx,ny),mut(nx,ny,nz),Vmean(nx,ny,nz,3),rhomean(nx,ny,nz),Tmean(nx,ny,nz))

      ! set Q / Jacobian
      do j = 1, ny
        do i = 1, nx
          Q(i,j,:,:) = Q(i,j,:,:) / Jacobian_cpu(i,j)
      enddo;enddo

      ! print initial condition
      if (id_turbulence == 0) then
        call print_vtk(0,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),ke0,entropy0)
      else
        call print_vtk(0,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),ke0,entropy0,real(mut_cpu))
      endif

      ! copy on GPU
      QJ = Q
      mut = 0.d0
      Rs(:,:,:,:) = 0.d0
      xix = xix_cpu
      etay = etay_cpu
      dx = dx_cpu
      dy = dy_cpu
      Jacobian = Jacobian_cpu
    endif

    ! share necessary data
    call MPI_BCAST(ke0, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(entropy0, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    do t2 = 1, np
      if (myrank == 0) then
        do t1 = 1, nt
          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJ,mut,E,F,G)
          call calc_step(nx,ny,nz,0.5d0,1.d0,dx,dy,E,F,G,QJ,QJs,Rs) ! QJs = Q2
          call set_bc(nx,ny,nz,Jacobian,QJs)

          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
          call calc_step(nx,ny,nz,0.5d0,2.d0,dx,dy,E,F,G,QJ,QJs,Rs) ! QJs = Q3
          call set_bc(nx,ny,nz,Jacobian,QJs)

          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
          call calc_step(nx,ny,nz,1.0d0,2.d0,dx,dy,E,F,G,QJ,QJs,Rs) ! QJs = Q4
          call set_bc(nx,ny,nz,Jacobian,QJs)

          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
          call calc_step4(nx,ny,nz,dx,dy,E,F,G,Rs,QJ)
          call set_bc(nx,ny,nz,Jacobian,QJ)
          call calc_mean(t1+(t2-1)*nt,nx,ny,nz,Jacobian,QJ,Vmean,rhomean,Tmean)
        enddo
      endif

      ! send and recv device arrays
      if (myrank == 0) then
        Q = QJ
        if (id_turbulence /= 0) then
          mut_cpu = mut
          call MPI_SEND(mut_cpu, nx*ny*nz, MPI_REAL8, 1, 1, MPI_COMM_WORLD, ierr)
        endif
        call MPI_SEND(Q, nx*ny*nz*5, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
      elseif (myrank == 1) then
        call MPI_RECV(Q, nx*ny*nz*5, MPI_REAL8, 0, 0, MPI_COMM_WORLD, status, ierr)
        if (id_turbulence == 0) then
          call print_vtk(t2,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),ke0,entropy0)
        else
          call MPI_RECV(mut_cpu, nx*ny*nz, MPI_REAL8, 0, 1, MPI_COMM_WORLD, status, ierr)
          call print_vtk(t2,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),ke0,entropy0,real(mut_cpu))
        endif
      endif
    enddo

    if (myrank == 0) then
      deallocate(QJ,QJs,Rs,E,F,G,dx,xix,dy,etay,Jacobian,mut,Vmean,rhomean,Tmean)
    endif
  end subroutine RungeKutta_4th

  subroutine RungeKutta_10th(id_RungeKutta,myrank,nx,ny,nz,x,dx_cpu,xix_cpu,y,dy_cpu,etay_cpu,z,Jacobian_cpu,Q)
    integer(kind=8), intent(in) :: id_RungeKutta  
    integer, intent(in)         :: myrank, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1), xix_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1), etay_cpu(ny-1)
    real(8), intent(in)         :: z(nz), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    integer i, j, k, t1, t2, itr, ierr, ilen, stat, status(MPI_STATUS_SIZE)
    real(8), dimension(nx,ny,nz)    :: mut_cpu
    real(8), dimension(nx,ny,nz,3)  :: Vmean_cpu
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)          :: prop
    real(8), allocatable, device  :: QJ(:,:,:,:), QJs(:,:,:,:), QJ4(:,:,:,:), R4(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device  :: dx(:), xix(:), dy(:), etay(:), Jacobian(:,:), mut(:,:,:), Vmean(:,:,:,:), rhomean(:,:,:), Tmean(:,:,:)
    ! for plot
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    ! check GPU
    stat = cudaSetDevice(0)
    stat = cudaGetDeviceProperties(prop,0)
    ilen = verify(prop%name, ' ', .true.)

    if (myrank == 0) then
      allocate(QJ(nx,ny,nz,5),QJs(nx,ny,nz,5),QJ4(nx,ny,nz,5),R4(nx,ny,nz,5),E(nx-1,ny-2,nz-2,5),F(nx-2,ny-1,nz-2,5),G(nx-2,ny-2,nz-1,5))
      allocate(dx(nx-1),xix(nx-1),dy(ny-1),etay(ny-1),Jacobian(nx,ny),mut(nx,ny,nz),Vmean(nx,ny,nz,3),rhomean(nx,ny,nz),Tmean(nx,ny,nz))

      ! set Q / Jacobian
      do j = 1, ny
        do i = 1, nx
          Q(i,j,:,:) = Q(i,j,:,:) / Jacobian_cpu(i,j)
      enddo;enddo

      ! print initial condition
      if (id_turbulence == 0) then
        call print_vtk(0,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),ke0,entropy0)
      else
        call print_vtk(0,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),ke0,entropy0,real(mut_cpu))
      endif

      ! copy on GPU
      QJ = Q
      mut = 0.d0
      xix = xix_cpu
      etay = etay_cpu
      dx = 1.d0 / dx_cpu
      dy = 1.d0 / dy_cpu
      Jacobian = Jacobian_cpu
    endif

    ! share necessary data
    call MPI_BCAST(ke0, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(entropy0, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    do t2 = 1, np
      if (myrank == 0) then
        do t1 = 1, nt
          QJs = QJ
          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
          call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
          call set_bc(nx,ny,nz,Jacobian,QJs)

          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
          call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
          call set_bc(nx,ny,nz,Jacobian,QJs)

          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
          call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
          call set_bc(nx,ny,nz,Jacobian,QJs)

          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
          call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
          call set_bc(nx,ny,nz,Jacobian,QJs)
          QJ4 = QJs
  
          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
          call calc_step2(nx,ny,nz,9.d0,6.d0,1.d0,15.d0,dx,dy,E,F,G,QJ,QJ4,QJs,R4)
          call set_bc(nx,ny,nz,Jacobian,QJs)

          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
          call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
          call set_bc(nx,ny,nz,Jacobian,QJs)

          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
          call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
          call set_bc(nx,ny,nz,Jacobian,QJs)

          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
          call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
          call set_bc(nx,ny,nz,Jacobian,QJs)

          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
          call calc_step5(nx,ny,nz,1.d0/6.d0,dx,dy,E,F,G,QJs)
          call set_bc(nx,ny,nz,Jacobian,QJs)

          call calc_EFG(id_hybrid,nx,ny,nz,xix,etay,Jacobian,QJs,mut,E,F,G)
          call calc_step10(nx,ny,nz,dx,dy,E,F,G,R4,QJ4,QJs,QJ)
          call set_bc(nx,ny,nz,Jacobian,QJ)
          call calc_mean(t1+(t2-1)*nt,nx,ny,nz,Jacobian,QJ,Vmean,rhomean,Tmean)
        enddo
      endif

      ! send and recv device arrays
      if (myrank == 0) then
        Q = QJ
        if (id_turbulence /= 0) then
          mut_cpu = mut
          call MPI_SEND(mut_cpu, nx*ny*nz, MPI_REAL8, 1, 1, MPI_COMM_WORLD, ierr)
        endif
        call MPI_SEND(Q, nx*ny*nz*5, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
      elseif (myrank == 1) then
        call MPI_RECV(Q, nx*ny*nz*5, MPI_REAL8, 0, 0, MPI_COMM_WORLD, status, ierr)
        if (id_turbulence == 0) then
          call print_vtk(t2,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),ke0,entropy0)
        else
          call MPI_RECV(mut_cpu, nx*ny*nz, MPI_REAL8, 0, 1, MPI_COMM_WORLD, status, ierr)
          call print_vtk(t2,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),ke0,entropy0,real(mut_cpu))
        endif 
      endif
    enddo

    if (myrank == 0) then
      deallocate(QJ,QJs,QJ4,R4,E,F,G,dx,xix,dy,etay,Jacobian,mut,Vmean,rhomean,Tmean)
    endif
  end subroutine RungeKutta_10th
end module calc_time_dev

