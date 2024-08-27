module calc_time_dev
  use cudafor
  use mpi
  use nvtx
  use mod_globals, only : accuracy, id_scheme, id_visc, id_turbulence, id_rescale, nt, np, nre, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG
  use calc_physical_quantities
  use calc_steps
  use calc_hybrid
  use calc_flux
  use calc_visc
  use calc_les
  use calc_rescale
  use set
  use print
  implicit none
  interface RungeKutta
    module procedure RungeKutta_3rd, RungeKutta_4th
  end interface
contains
  subroutine calc_EFG(nx,ny,nz,dx,dy,dz,Jacobian,QJ,mut,E,F,G,sensor)
    integer, intent(in), value                          :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device        :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1), device        :: dy ! 1 / dy
    real(8), intent(in), dimension(nz-1), device        :: dz ! 1 / dz
    real(8), intent(in), dimension(nx,ny,nz), device    :: Jacobian
    real(8), intent(in), dimension(nx,ny,nz,5), device  :: QJ ! Q / Jacobian
    real(8), intent(inout), dimension(nx,ny,nz), device :: mut
    real(8), intent(out), device                        :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device                        :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device                        :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), intent(out), dimension(nx,ny,nz), device   :: sensor
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p, T, fd, qc2
    integer stat
    call calc_quantities(nx,ny,nz,Jacobian,QJ,rho,u,v,w,p,T)
    
    if (id_scheme == 2 .or. 4 <= id_scheme) then
      call calc_Ducros<<<blocks,threads>>>(nx,ny,nz,dx,dy,dz,u,v,w,rho,p,fd)
    endif
    call calc_E<<<blocksE,threadsE,1>>>(nx,ny,nz,rho,u,v,w,p,fd,E)
    call calc_F<<<blocksF,threadsF,2>>>(nx,ny,nz,rho,u,v,w,p,fd,F)
    call calc_G<<<blocksG,threadsG,3>>>(nx,ny,nz,rho,u,v,w,p,fd,G)
  
    qc2 = 0.d0
    if (id_turbulence /= 0) then
      call calc_mut<<<blocks,threads>>>(nx,ny,nz,dx,dy,dz,rho,u,v,w,mut,qc2)
      call set_bc_mut(nx,ny,nz,mut,qc2)
    endif

    stat = cudaDeviceSynchronize()
    if (1 <= id_visc .or. id_turbulence /= 0) then
      call calc_Ev<<<blocksE,threadsE,1>>>(nx,ny,nz,dx,dy,dz,rho,u,v,w,T,p,mut,qc2,E)
      call calc_Fv<<<blocksF,threadsF,2>>>(nx,ny,nz,dy,dx,dz,rho,u,v,w,T,p,mut,qc2,F)
      call calc_Gv<<<blocksG,threadsG,3>>>(nx,ny,nz,dx,dy,dz,rho,u,v,w,T,p,mut,qc2,G)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG

  subroutine RungeKutta_3rd(id_RungeKutta,myrank,nx,ny,nz,x,dx_cpu,xix_cpu,y,dy_cpu,etay_cpu,z,dz_cpu,zetaz_cpu,Jacobian_cpu,Q)
    integer(kind=2), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1), xix_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1), etay_cpu(ny-1)
    real(8), intent(in)         :: z(nz), dz_cpu(nz-1), zetaz_cpu(nz-1), Jacobian_cpu(nx,ny,nz)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    integer i, j, k, t1, t2, itr, ilen, ierr, stat, request, status(MPI_STATUS_SIZE)
    real(8), dimension(nx,ny,nz)    :: mut_cpu, sensor_cpu
    real(8), dimension(nx,ny,nz,3)  :: Vmean_cpu
    ! rescale !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: Qre(:,:,:,:)
    real(8), allocatable :: Qre_cpu(:,:,:,:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)          :: prop
    real(8), allocatable, device  :: QJ(:,:,:,:), QJ2(:,:,:,:), QJ3(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device  :: dx(:), xix(:), dy(:), etay(:), dz(:), zetaz(:), Jacobian(:,:,:), mut(:,:,:), sensor(:,:,:)
    ! for plot
    real(4) :: mass0 = 1.d0, ke0 = 1.d0, entropy0 = 1.d0

    ! check GPU
    stat = cudaSetDevice(0)
    stat = cudaGetDeviceProperties(prop,0)
    ilen = verify(prop%name, ' ', .true.)
    print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", 0, ") is available"

    if (myrank == 0) then
      allocate(QJ(nx,ny,nz,5),QJ2(nx,ny,nz,5),QJ3(nx,ny,nz,5),E(nx-1,ny-2,nz-2,5),F(nx-2,ny-1,nz-2,5),G(nx-2,ny-2,nz-1,5))
      allocate(dx(nx-1),xix(nx-1),dy(ny-1),etay(ny-1),dz(nz-1),zetaz(nz-1),Jacobian(nx,ny,nz),mut(nx,ny,nz),sensor(nx,ny,nz))
      ! rescale
      allocate(Qre(2,ny,nz,5))

      ! set Q / Jacobian
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            Q(i,j,k,:) = Q(i,j,k,:) / Jacobian_cpu(i,j,k)
      enddo;enddo;enddo
      sensor_cpu = 0.d0

      ! print initial condition
      call print_vtk(0,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),real(sensor_cpu),mass0,ke0,entropy0)

      ! copy on GPU
      QJ       = Q
      mut      = 0.d0
      xix      = xix_cpu
      etay     = etay_cpu
      zetaz    = zetaz_cpu
      dx       = dx_cpu
      dy       = dy_cpu
      dz       = dz_cpu
      Jacobian = Jacobian_cpu
      sensor   = sensor_cpu
    endif

    ! share necessary data
    call MPI_BCAST(mass0,    1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ke0,      1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(entropy0, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    ! rescale
    allocate(Qre_cpu(2,ny,nz,5))

    do t2 = 1, np
      if (myrank == 0) then
        do t1 = 1, nt
          call nvtxStartRange("calc 1step",1)
          call nvtxStartRange("calc flux",2)
          call calc_EFG(nx,ny,nz,xix,etay,zetaz,Jacobian,QJ,mut,E,F,G,sensor)
          call nvtxEndRange
          call nvtxStartRange("calc time dev",3)
          call calc_step(nx,ny,nz,1.d0,0.d0,dx,dy,dz,E,F,G,QJ,QJ2)
          call nvtxEndRange
          call set_bc(nx,ny,nz,Jacobian,QJ2,Qre)

          call calc_EFG(nx,ny,nz,xix,etay,zetaz,Jacobian,QJ2,mut,E,F,G,sensor)
          call calc_step2(nx,ny,nz,0.75d0,0.25d0,0.25d0,1.d0,dx,dy,dz,E,F,G,QJ,QJ2,QJ3)
          call set_bc(nx,ny,nz,Jacobian,QJ3,Qre)

          call calc_EFG(nx,ny,nz,xix,etay,zetaz,Jacobian,QJ3,mut,E,F,G,sensor)
          call calc_step3(nx,ny,nz,dx,dy,dz,E,F,G,QJ3,QJ)
          call set_bc(nx,ny,nz,Jacobian,QJ,Qre)
          call nvtxEndRange
        enddo
      endif

      ! send and recv device arrays
      if (myrank == 0) then
        Q = QJ
        sensor_cpu = sensor
        call nvtxStartRange("MPI_SEND",4)
        call MPI_SEND(Q,          nx*ny*nz*5, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
        call MPI_SEND(sensor_cpu, nx*ny*nz,   MPI_REAL8, 1, 1, MPI_COMM_WORLD, ierr) 
        call nvtxEndRange
      elseif (myrank == 1) then
        call nvtxStartRange("MPI_RECV",5)
        call MPI_RECV(Q,          nx*ny*nz*5, MPI_REAL8, 0, 0, MPI_COMM_WORLD, status, ierr)
        call MPI_RECV(sensor_cpu, nx*ny*nz,   MPI_REAL8, 0, 1, MPI_COMM_WORLD, status, ierr)
        call nvtxEndRange
        call nvtxStartRange("print",6)
        call print_vtk(t2,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),real(sensor_cpu),mass0,ke0,entropy0)
        call nvtxEndRange
      endif
    enddo
    
    if (myrank == 0) then
      deallocate(QJ,QJ2,QJ3,E,F,G,dx,xix,dy,etay,dz,zetaz,Jacobian,mut,sensor,Qre)
    endif
  end subroutine RungeKutta_3rd

  subroutine RungeKutta_4th(id_RungeKutta,myrank,nx,ny,nz,x,dx_cpu,xix_cpu,y,dy_cpu,etay_cpu,z,dz_cpu,zetaz_cpu,Jacobian_cpu,Q)
    integer(kind=4), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1), xix_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1), etay_cpu(ny-1)
    real(8), intent(in)         :: z(nz), dz_cpu(nz-1), zetaz_cpu(nz-1), Jacobian_cpu(nx,ny,nz)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    integer i, j, k, t1, t2, itr, ierr, ilen, stat, ireq, ireq2, istat(MPI_STATUS_SIZE), status(MPI_STATUS_SIZE)
    real(8), dimension(nx,ny,nz)    :: mut_cpu, sensor_cpu
    real(8), dimension(nx,ny,nz,3)  :: Vmean_cpu
    ! rescale !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: Qre(:,:,:,:)
    real(8), allocatable :: Qre_cpu(:,:,:,:), Um(:,:), Vm(:,:), Wm(:,:), pm(:,:), Tm(:,:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)          :: prop
    real(8), allocatable, device  :: QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device  :: dx(:), xix(:), dy(:), etay(:), dz(:), zetaz(:), Jacobian(:,:,:), mut(:,:,:), sensor(:,:,:)
    ! for plot
    real(4) :: mass0 = 1.d0, ke0 = 1.d0, entropy0 = 1.d0

    ! check GPU
    stat = cudaSetDevice(0)
    stat = cudaGetDeviceProperties(prop,0)
    ilen = verify(prop%name, ' ', .true.)

    if (myrank == 0) then
      allocate(QJ(nx,ny,nz,5),QJs(nx,ny,nz,5),Rs(nx,ny,nz,5),E(nx-1,ny-2,nz-2,5),F(nx-2,ny-1,nz-2,5),G(nx-2,ny-2,nz-1,5))
      allocate(dx(nx-1),xix(nx-1),dy(ny-1),etay(ny-1),dz(nz-1),zetaz(nz-1),Jacobian(nx,ny,nz),mut(nx,ny,nz),sensor(nx,ny,nz))
      ! rescale
      allocate(Qre(2,ny,nz,5))

      ! set Q / Jacobian
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            Q(i,j,k,:) = Q(i,j,k,:) / Jacobian_cpu(i,j,k)
      enddo;enddo;enddo
      sensor_cpu = 0.d0

      ! print initial condition
      call print_vtk(0,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),real(sensor_cpu),mass0,ke0,entropy0)

      ! copy on GPU
      QJ          = Q
      mut         = 0.d0
      Rs(:,:,:,:) = 0.d0
      xix         = xix_cpu
      etay        = etay_cpu
      zetaz       = zetaz_cpu
      dx          = dx_cpu
      dy          = dy_cpu
      dz          = dz_cpu
      Jacobian    = Jacobian_cpu
      sensor      = sensor_cpu
    endif

    ! share necessary data
    call MPI_BCAST(mass0,    1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ke0,      1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(entropy0, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    ! rescale
    allocate(Qre_cpu(2,ny,nz,5),Um(2,ny),Vm(2,ny),Wm(2,ny),pm(2,ny),Tm(2,ny))

    do t2 = 1, np
      do t1 = 1, nt
        if (myrank == 0) then
          call nvtxStartRange("Send Qre", 1)
          if (kind(id_rescale) == 4) then
            Qre_cpu(:,:,:,:) = QJ(nre+1:nre+2,:,:,:)
            call MPI_ISEND(Qre_cpu, 10*ny*nz, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ireq, ierr)
          endif
          call nvtxEndRange
          call nvtxStartRange("calc flux", 2)
          call calc_EFG(nx,ny,nz,xix,etay,zetaz,Jacobian,QJ,mut,E,F,G,sensor)
          call calc_step(nx,ny,nz,0.5d0,1.d0,dx,dy,dz,E,F,G,QJ,QJs,Rs) ! QJs = Q2
          call nvtxEndRange
          call nvtxStartRange("Recv Qre", 3)
          if (kind(id_rescale) == 4) then
            call MPI_RECV(Qre_cpu, 10*ny*nz, MPI_REAL8, 1, 1, MPI_COMM_WORLD, istat, ierr)
            Qre = Qre_cpu
          endif
          call nvtxEndRange
          call nvtxStartRange("set bc", 4)
          call set_bc(nx,ny,nz,Jacobian,QJs,Qre)
          call nvtxEndRange
        elseif (myrank == 1 .and. kind(id_rescale) == 4) then
          call MPI_IRECV(Qre_cpu, 10*ny*nz, MPI_REAL8, 0, 0, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAIT(ireq, istat, ierr)
          call nvtxStartRange("rescale", 5)
          call set_rescale(t1+(t2-1)*nt,nx,ny,nz,nre,2.d-3,y,Jacobian_cpu,Um,Vm,Wm,pm,Tm,Qre_cpu)
          call nvtxEndRange
          call MPI_SEND(Qre_cpu, 10*ny*nz, MPI_REAL8, 0, 1, MPI_COMM_WORLD, ierr)
        endif

        if (myrank == 0) then
          if (kind(id_rescale) == 4) then
            Qre_cpu(:,:,:,:) = QJs(nre+1:nre+2,:,:,:)
            call MPI_ISEND(Qre_cpu, 10*ny*nz, MPI_REAL8, 1, 2, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(nx,ny,nz,xix,etay,zetaz,Jacobian,QJs,mut,E,F,G,sensor)
          call calc_step(nx,ny,nz,0.5d0,2.d0,dx,dy,dz,E,F,G,QJ,QJs,Rs) ! QJs = Q3
          if (kind(id_rescale) == 4) then
            call MPI_RECV(Qre_cpu, 10*ny*nz, MPI_REAL8, 1, 3, MPI_COMM_WORLD, istat, ierr)
            Qre = Qre_cpu
          endif
          call set_bc(nx,ny,nz,Jacobian,QJs,Qre)
        elseif (myrank == 1 .and. kind(id_rescale) == 4) then
          call MPI_IRECV(Qre_cpu, 10*ny*nz, MPI_REAL8, 0, 2, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAit(ireq, istat, ierr)
          call set_rescale(t1+(t2-1)*nt,nx,ny,nz,nre,2.d-3,y,Jacobian_cpu,Um,Vm,Wm,pm,Tm,Qre_cpu)
          call MPI_SEND(Qre_cpu, 10*ny*nz, MPI_REAL8, 0, 3, MPI_COMM_WORLD, ierr)
        endif

        if (myrank == 0) then
          if (kind(id_rescale) == 4) then
            Qre_cpu(:,:,:,:) = QJs(nre+1:nre+2,:,:,:)
            call MPI_ISEND(Qre_cpu, 10*ny*nz, MPI_REAL8, 1, 4, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(nx,ny,nz,xix,etay,zetaz,Jacobian,QJs,mut,E,F,G,sensor)
          call calc_step(nx,ny,nz,1.0d0,2.d0,dx,dy,dz,E,F,G,QJ,QJs,Rs) ! QJs = Q4
          if (kind(id_rescale) == 4) then
            call MPI_RECV(Qre_cpu, 10*ny*nz, MPI_REAL8, 1, 5, MPI_COMM_WORLD, istat, ierr)
            Qre = Qre_cpu
          endif
          call set_bc(nx,ny,nz,Jacobian,QJs,Qre)
        elseif (myrank == 1 .and. kind(id_rescale) == 4) then
          call MPI_IRECV(Qre_cpu, 10*ny*nz, MPI_REAL8, 0, 4, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAIT(ireq, istat, ierr)
          call set_rescale(t1+(t2-1)*nt,nx,ny,nz,nre,2.d-3,y,Jacobian_cpu,Um,Vm,Wm,pm,Tm,Qre_cpu)
          call MPI_SEND(Qre_cpu, 10*ny*nz, MPI_REAL8, 0, 5, MPI_COMM_WORLD, ierr)
        endif

        if (myrank == 0) then
          if (kind(id_rescale) == 4) then
            Qre_cpu(:,:,:,:) = QJs(nre+1:nre+2,:,:,:)
            call MPI_ISEND(Qre_cpu, 10*ny*nz, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(nx,ny,nz,xix,etay,zetaz,Jacobian,QJs,mut,E,F,G,sensor)
          call calc_step4(nx,ny,nz,dx,dy,dz,E,F,G,Rs,QJ)
          if (kind(id_rescale) == 4) then
            call MPI_RECV(Qre_cpu, 10*ny*nz, MPI_REAL8, 1, 1, MPI_COMM_WORLD, istat, ierr)
            Qre = Qre_cpu
          endif
          call set_bc(nx,ny,nz,Jacobian,QJ,Qre)
        elseif (myrank == 1 .and. kind(id_rescale) == 4) then
          call MPI_IRECV(Qre_cpu, 10*ny*nz, MPI_REAL8, 0, 0, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAIT(ireq, istat, ierr)
          call set_rescale(t1+(t2-1)*nt,nx,ny,nz,nre,2.d-3,y,Jacobian_cpu,Um,Vm,Wm,pm,Tm,Qre_cpu)
          call MPI_SEND(Qre_cpu, 10*ny*nz, MPI_REAL8, 0, 1, MPI_COMM_WORLD, ierr)
        endif
      enddo

      ! send and recv device arrays
      if (myrank == 0) then
        Q          = QJ
        sensor_cpu = sensor
        call MPI_SEND(Q,          nx*ny*nz*5, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
        call MPI_SEND(sensor_cpu, nx*ny*nz,   MPI_REAL8, 1, 1, MPI_COMM_WORLD, ierr) 
      elseif (myrank == 1) then
        call MPI_RECV(Q,          nx*ny*nz*5, MPI_REAL8, 0, 0, MPI_COMM_WORLD, istat, ierr)
        call MPI_RECV(sensor_cpu, nx*ny*nz,   MPI_REAL8, 0, 1, MPI_COMM_WORLD, istat, ierr)
        call print_vtk(t2,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),real(sensor_cpu),mass0,ke0,entropy0)
      endif
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (myrank == 0) then
      deallocate(QJ,QJs,Rs,E,F,G,dx,xix,dy,etay,dz,zetaz,Jacobian,mut,sensor,Qre)
    endif
    deallocate(Qre_cpu,Um,Vm,Wm,pm,Tm)
  end subroutine RungeKutta_4th
end module calc_time_dev

