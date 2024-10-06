module calc_time_dev2
  use cudafor
  use mpi
  use nvtx
  use mod_globals, only : accuracy, id_rescale, nt, np, nre
  use calc_physical_quantities
  use calc_steps
  use calc_hybrid
  use calc_flux
  use calc_visc
  use calc_rescale
  use set
  use print
  implicit none
contains
  subroutine calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,dx,dy,dz,Jacobian,QJ,E,F,G,sensor)
    integer, intent(in), value                         :: nx, ny, nz
    type(dim3), intent(in)                             :: blocksE, blocksF, blocksG, blocks
    type(dim3), intent(in)                             :: threadsE, threadsF, threadsG, threads
    real(8), intent(in), dimension(nx-1), device       :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1), device       :: dy ! 1 / dy
    real(8), intent(in), dimension(nz-1), device       :: dz ! 1 / dz
    real(8), intent(in), dimension(nx,ny,nz), device   :: Jacobian
    real(8), intent(in), dimension(nx,ny,nz,5), device :: QJ ! Q / Jacobian
    real(8), intent(out), device                       :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device                       :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device                       :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), intent(out), dimension(nx,ny,nz), device  :: sensor
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p, T
    integer stat
    call calc_quantities(nx,ny,nz,Jacobian,QJ,rho,u,v,w,p,T)

    call calc_Ducros<<<blocks,threads>>>(nx,ny,nz,dx,dy,dz,u,v,w,rho,p,sensor)
    call calc_E<<<blocksE,threadsE,1>>>(nx,ny,nz,rho,u,v,w,p,sensor,E)
    call calc_F<<<blocksF,threadsF,2>>>(nx,ny,nz,rho,u,v,w,p,sensor,F)
    call calc_G<<<blocksG,threadsG,3>>>(nx,ny,nz,rho,u,v,w,p,sensor,G)
  
    stat = cudaDeviceSynchronize()
    call calc_Ev<<<blocksE,threadsE,1>>>(nx,ny,nz,dx,dy,dz,rho,u,v,w,T,p,E)
    call calc_Fv<<<blocksF,threadsF,2>>>(nx,ny,nz,dy,dx,dz,rho,u,v,w,T,p,F)
    call calc_Gv<<<blocksG,threadsG,3>>>(nx,ny,nz,dx,dy,dz,rho,u,v,w,T,p,G)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG

  subroutine RungeKutta(myrank,nx,ny,nz,x,dx_cpu,xix_cpu,y,dy_cpu,etay_cpu,z,dz_cpu,zetaz_cpu,Jacobian_cpu,Q)
    integer, intent(in)    :: myrank, nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1), xix_cpu(nx-1)
    real(8), intent(in)    :: y(ny), dy_cpu(ny-1), etay_cpu(ny-1)
    real(8), intent(in)    :: z(nz), dz_cpu(nz-1), zetaz_cpu(nz-1), Jacobian_cpu(nx,ny,nz)
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer i, j, k, t1, t2, itr, ilen, ierr, stat, ndevices, ireq, ireq2, istat(MPI_STATUS_SIZE)
    real(8), dimension(nx,ny,nz) :: sensor_cpu
    ! rescale
    real(8), allocatable, pinned :: Qre_cpu(:,:,:,:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)         :: prop
    type(dim3)                   :: blocksE, blocksF, blocksG, blocks
    type(dim3)                   :: threadsE, threadsF, threadsG, threads
    real(8), allocatable, device :: QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: dx(:), xix(:), dy(:), etay(:), dz(:), zetaz(:), Jacobian(:,:,:), sensor(:,:,:)
    ! for plot
    real(8) :: mass0 = 1.d0, ke0 = 1.d0, entropy0 = 1.d0
    real(8), allocatable :: rhomt(:), pmt(:), Tmt(:), Mmt(:), vmt(:)

    ! check GPU
    if (mod(myrank,2) == 0) then
      if (myrank == 0) then
        stat = cudaGetDeviceCount(ndevices)
        print '(2x, i2, a)', ndevices, " GPU devices are found"
      endif
      stat = cudaSetDevice(myrank/2)
      stat = cudaGetDeviceProperties(prop,myrank/2)
      ilen = verify(prop%name, ' ', .true.)
      print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", myrank/2, ") is available"
    endif

    allocate(rhomt(nx*ny*nz), pmt(nx*ny*nz), Tmt(nx*ny*nz), Mmt(nx*ny*nz), vmt(3*nx*ny*nz))
    rhomt(:) = 0.d0
    pmt(:)   = 0.d0
    Tmt(:)   = 0.d0
    Mmt(:)   = 0.d0
    vmt(:)   = 0.d0

    if (mod(myrank,2) == 0) then
      allocate(QJ(nx,ny,nz,5),QJs(nx,ny,nz,5),Rs(nx,ny,nz,5),E(nx-1,ny-2,nz-2,5),F(nx-2,ny-1,nz-2,5),G(nx-2,ny-2,nz-1,5))
      allocate(dx(nx-1),xix(nx-1),dy(ny-1),etay(ny-1),dz(nz-1),zetaz(nz-1),Jacobian(nx,ny,nz),sensor(nx,ny,nz))
      call set_blocks_threads(myrank,nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads)

      ! set Q / Jacobian
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            Q(i,j,k,:) = Q(i,j,k,:) / Jacobian_cpu(i,j,k)
      enddo;enddo;enddo
      sensor_cpu = 0.d0

      ! print initial condition
      call print_vtk(0,nx,ny,nz,x,y,z,Jacobian_cpu,Q,rhomt,pmt,Tmt,Mmt,vmt,mass0,ke0,entropy0,myrank+1)

      ! copy on GPU
      QJ       = Q
      Rs       = 0.d0
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
    if (mod(myrank,2) == 0) then
      call MPI_SEND(mass0,    1, MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
      call MPI_SEND(ke0,      1, MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
      call MPI_SEND(entropy0, 1, MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
    else
      call MPI_RECV(mass0,    1, MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(ke0,      1, MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
    endif

    ! rescale
    if (myrank == 0 .or. myrank == 1) then
      allocate(Qre_cpu(10,ny,nz,5))
    endif

    do t2 = 1, np
      do t1 = 1, nt
        if (mod(myrank,2) == 0) then
          if (myrank == 0) then
            Qre_cpu(:,:,:,:) = QJ(nre:nre+9,:,:,:)
            call MPI_ISEND(Qre_cpu, 50*ny*nz, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,xix,etay,zetaz,Jacobian,QJ,E,F,G,sensor)
          call calc_step(nx,ny,nz,0.5d0,1.d0,dx,dy,dz,E,F,G,QJ,QJs,Rs) ! QJs = Q2
          if (myrank == 0) then
            call MPI_RECV(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 1, 1, MPI_COMM_WORLD, istat, ierr)
            QJs(1,:,:,:) = Qre_cpu(1,:,:,:)
          endif
          call set_bc(myrank,nx,ny,nz,Jacobian,QJs)
        elseif (myrank == 1) then
          call MPI_IRECV(Qre_cpu, 50*ny*nz, MPI_REAL8, 0, 0, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAIT(ireq, istat, ierr)
          call set_rescale(t1+(t2-1)*nt,nx,ny,nz,nre,2.d-3,y,Jacobian_cpu,Qre_cpu)
          call MPI_SEND(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 0, 1, MPI_COMM_WORLD, ierr)
        endif

        if (mod(myrank,2) == 0) then
          if (myrank == 0) then
            Qre_cpu(:,:,:,:) = QJ(nre:nre+9,:,:,:)
            call MPI_ISEND(Qre_cpu, 50*ny*nz, MPI_REAL8, 1, 2, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,xix,etay,zetaz,Jacobian,QJs,E,F,G,sensor)
          call calc_step(nx,ny,nz,0.5d0,2.d0,dx,dy,dz,E,F,G,QJ,QJs,Rs) ! QJs = Q3
          if (myrank == 0) then
            call MPI_RECV(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 1, 3, MPI_COMM_WORLD, istat, ierr)
            QJs(1,:,:,:) = Qre_cpu(1,:,:,:)
          endif
          call set_bc(myrank,nx,ny,nz,Jacobian,QJs)
        elseif (myrank == 1) then
          call MPI_IRECV(Qre_cpu, 50*ny*nz, MPI_REAL8, 0, 2, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAIT(ireq, istat, ierr)
          call set_rescale(t1+(t2-1)*nt,nx,ny,nz,nre,2.d-3,y,Jacobian_cpu,Qre_cpu)
          call MPI_SEND(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 0, 3, MPI_COMM_WORLD, ierr)
        endif

        if (mod(myrank,2) == 0) then
          if (myrank == 0) then
            Qre_cpu(:,:,:,:) = QJ(nre:nre+9,:,:,:)
            call MPI_ISEND(Qre_cpu, 50*ny*nz, MPI_REAL8, 1, 4, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,xix,etay,zetaz,Jacobian,QJs,E,F,G,sensor)
          call calc_step(nx,ny,nz,1.d0,2.d0,dx,dy,dz,E,F,G,QJ,QJs,Rs) ! QJs = Q4
          if (myrank == 0) then
            call MPI_RECV(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 1, 5, MPI_COMM_WORLD, istat, ierr)
            QJs(1,:,:,:) = Qre_cpu(1,:,:,:)
          endif
          call set_bc(myrank,nx,ny,nz,Jacobian,QJs)
        elseif (myrank == 1) then
          call MPI_IRECV(Qre_cpu, 50*ny*nz, MPI_REAL8, 0, 4, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAIT(ireq, istat, ierr)
          call set_rescale(t1+(t2-1)*nt,nx,ny,nz,nre,2.d-3,y,Jacobian_cpu,Qre_cpu)
          call MPI_SEND(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 0, 5, MPI_COMM_WORLD, ierr)
        endif

        if (mod(myrank,2) == 0) then
          if (myrank == 0) then
            Qre_cpu(:,:,:,:) = QJ(nre:nre+9,:,:,:)
            call MPI_ISEND(Qre_cpu, 50*ny*nz, MPI_REAL8, 1, 6, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,xix,etay,zetaz,Jacobian,QJs,E,F,G,sensor)
          call calc_step4(nx,ny,nz,dx,dy,dz,E,F,G,Rs,QJ)
          if (myrank == 0) then
            call MPI_RECV(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 1, 7, MPI_COMM_WORLD, istat, ierr)
            QJ(1,:,:,:) = Qre_cpu(1,:,:,:)
          endif
          call set_bc(myrank,nx,ny,nz,Jacobian,QJ)
        elseif (myrank == 1) then
          call MPI_IRECV(Qre_cpu, 50*ny*nz, MPI_REAL8, 0, 6, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAIT(ireq, istat, ierr)
          call set_rescale(t1+(t2-1)*nt,nx,ny,nz,nre,2.d-3,y,Jacobian_cpu,Qre_cpu)
          call MPI_SEND(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 0, 7, MPI_COMM_WORLD, ierr)
        endif
      enddo

      ! send and recv device arrays
      if (mod(myrank,2) == 0) then
        Q = QJ
        call MPI_SEND(Q, nx*ny*nz*5, MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr) 
      else
        call MPI_RECV(Q, nx*ny*nz*5, MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
        call print_vtk(t2,nx,ny,nz,x,y,z,Jacobian_cpu,Q,rhomt,pmt,Tmt,Mmt,vmt,mass0,ke0,entropy0,myrank)
      endif
    enddo
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (mod(myrank,2) == 0) then
      deallocate(QJ,QJs,Rs,E,F,G,dx,xix,dy,etay,dz,zetaz,Jacobian,sensor)
    endif
    if (myrank == 0 .or. myrank == 1) then
      deallocate(Qre_cpu)
    endif
    deallocate(rhomt,pmt,Tmt,Mmt,vmt)
  end subroutine RungeKutta
end module calc_time_dev2

