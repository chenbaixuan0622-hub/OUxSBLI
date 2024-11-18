module calc_time_dev2
  use, intrinsic :: iso_fortran_env
  use cudafor
  use mpi
  use nvtx
  use mod_globals, only : accuracy, id_rescale, nt, np, nre1, rerank, overlap, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use calc_physical_quantities
  use calc_steps
  use calc_hybrid
  use calc_flux
  use calc_visc
  use calc_rescale
  use calc_para
  use set
  use print
  implicit none
contains
  subroutine calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,dx,dy,dz,Jacobian,QJ,E,F,G)
    integer, intent(in), value   :: nx, ny, nz
    type(dim3), intent(in)       :: blocksE, blocksF, blocksG, blocks
    type(dim3), intent(in)       :: threadsE, threadsF, threadsG, threads
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(ny)
    real(8), intent(in), device  :: QJ(nx,ny,nz,5) ! Q / Jacobian
    real(8), intent(out), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p, sensor
    integer stat
    call calc_quantities_3D(nx,ny,nz,Jacobian,QJ,rho,u,v,w,p)
    
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, u, v, w, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(nx, ny, nz, rho, u, v, w, p, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(nx, ny, nz, rho, u, v, w, p, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(nx, ny, nz, rho, u, v, w, p, sensor, G)
  
    stat = cudaDeviceSynchronize()
    call calc_Ev<<<blocksE,threadsE,1>>>(nx,ny,nz,dx,dy,dz,rho,u,v,w,p,E)
    call calc_Fv<<<blocksF,threadsF,2>>>(nx,ny,nz,dy,dx,dz,rho,u,v,w,p,F)
    call calc_Gv<<<blocksG,threadsG,3>>>(nx,ny,nz,dx,dy,dz,rho,u,v,w,p,G)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG

  subroutine RungeKutta(myrank,nx,ny,nz,x,dx_cpu,y,dy_cpu,z,dz_cpu,Jacobian_cpu,Q)
    integer, intent(in)    :: myrank, nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)    :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)    :: z(nz), dz_cpu(nz-1), Jacobian_cpu(ny)
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer i, j, k, t1, t2, itr, ilen, ierr, stat, nranks, ndevices, ireq, ireqs(2), istat(MPI_STATUS_SIZE), istats(MPI_STATUS_SIZE,2)
    integer :: canaccess1, canaccess2, canaccess = 1
    ! rescale
    real(8), allocatable, device :: Qre(:,:,:), Qm(:,:)
    real(8), allocatable         :: Qre_cpu(:,:,:), Qm_cpu(:,:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)         :: prop
    real(8), allocatable, device :: QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), zetaz(:), Jacobian(:)
    ! for plot
    real(8) :: mass0 = 1.d0, ke0 = 1.d0, entropy0 = 1.d0
    real(8) t_start, t_end

    ! check GPU
    if (mod(myrank,2) == 0) then
      if (myrank == 0) then
        stat = cudaGetDeviceCount(ndevices)
        print '(2x, i2, a)', ndevices, " GPU devices are found"
      endif
     
      call MPI_BCAST(ndevices, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (ndevices /= 0) then
        stat = cudaSetDevice(myrank/2)
        stat = cudaGetDeviceProperties(prop,myrank/2)
        ilen = verify(prop%name, ' ', .true.)
        print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", myrank/2, ") is available"
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
        if (nranks == 8*2) then
          ! check peer to peer access
          if (2 <= myrank .and. myrank <= 12) then
            istat = cudaDeviceCanAccessPeer(canaccess1, myrank/2, myrank/2-1)
            istat = cudaDeviceCanAccessPeer(canaccess2, myrank/2, myrank/2+1)
          elseif (myrank == 0) then
            istat = cudaDeviceCanAccessPeer(canaccess1, 0, 1)
            istat = cudaDeviceCanAccessPeer(canaccess2, 0, 7)
          else
            istat = cudaDeviceCanAccessPeer(canaccess1, 7, 0)
            istat = cudaDeviceCanAccessPeer(canaccess2, 7, 6)
          endif

          ! enable peer to peer access
          if (canaccess1 == 1 .and. canaccess2 == 1) then
            if (2 <= myrank .and. myrank <= 12) then
              istat = cudaDeviceEnablePeerAccess(myrank/2-1,0)
              istat = cudaDeviceEnablePeerAccess(myrank/2+1,0)
            elseif (myrank == 0) then
              istat = cudaDeviceEnablePeerAccess(1,0)
              istat = cudaDeviceEnablePeerAccess(7,0)
            else
              istat = cudaDeviceEnablePeerAccess(0,0)
              istat = cudaDeviceEnablePeerAccess(6,0)
            endif
            print *, "GPU:", myrank/2, " Device Peer to Peer Access is available"
          else
            print *, "GPU:", myrank/2, " Device Peer to Peer Access is unavailable"
          endif
            
          if ((myrank == rerank .or. myrank == 0) .and. kind(id_rescale) == 4) then
            ! check peer to peer access for rescaling
            istat = cudaDeviceCanAccessPeer(canaccess1, 0, rerank/2)
            istat = cudaDeviceCanAccessPeer(canaccess2, rerank/2, 0)
              
            ! enable peer to peer access
            if (canaccess1 == 1 .and. canaccess2 == 1) then
              if (myrank == 0) then
                istat = cudaDeviceEnablePeerAccess(rerank/2,0)
              elseif (myrank == rerank) then
                istat = cudaDeviceEnablePeerAccess(0,0)
              endif
            else
              print *, "For rescaling GPU:", myrank/2, " Device Peer to Peer Access is available"
            endif
              print *, "For rescaling GPU:", myrank/2, " Device Peer to Peer Access is available"
          endif
        elseif (nranks == 2*2) then
          ! check peer to peer access
          istat = cudaDeviceCanAccessPeer(canaccess1, 0, 1)
          istat = cudaDeviceCanAccessPeer(canaccess2, 1, 0)

          ! enable peer to peer access
          if (canaccess1 == 1 .and. canaccess2 == 1) then
            if (myrank == 0) then
              istat = cudaDeviceEnablePeerAccess(1,0)
            else
              istat = cudaDeviceEnablePeerAccess(0,0)
            endif
            print *, "GPU:", myrank/2, " Device Peer to Peer Access is available"
          else
            print *, "GPU:", myrank/2, " Device Peer to Peer Access is unavailable"
          endif
        else
          print *, "wrong MPI_COMM_SIZE"
        endif
      else
        stat = cudaSetDevice(0)
        stat = cudaGetDeviceProperties(prop,0)
        ilen = verify(prop%name, ' ', .true.)
        print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", 0, ") is available"
      endif
    endif

    ! check peer to perr access
    if (mod(myrank,2) == 0) then
      if (canaccess1 == 0 .or. canaccess2 == 0) then
        canaccess = canaccess1 * canaccess2
        call MPI_BCAST(canaccess, 1, MPI_INTEGER, myrank, MPI_COMM_WORLD, ierr)
      endif
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (mod(myrank,2) == 0) then
      allocate(QJ(nx,ny,nz,5),QJs(nx,ny,nz,5),Rs(nx-2,ny-2,nz-2,5),E(nx-1,ny-2,nz-2,5),F(nx-2,ny-1,nz-2,5),G(nx-2,ny-2,nz-1,5))
      allocate(xix(nx-1),etay(ny-1),zetaz(nz-1),Jacobian(ny))

      print *, "myrank is ", myrank, " memory allocation has completed"

      ! set Q / Jacobian
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            Q(i,j,k,:) = Q(i,j,k,:) / Jacobian_cpu(j)
      enddo;enddo;enddo

      ! print initial condition
      call print_vtk(0,nx,ny,nz,x,y,z,Jacobian_cpu,Q,mass0,ke0,entropy0,myrank+1)

      ! copy on GPU
      QJ       = Q
      Rs       = 0.d0
      xix      = 1.d0 / dx_cpu
      etay     = 1.d0 / dy_cpu
      zetaz    = 1.d0 / dz_cpu
      Jacobian = Jacobian_cpu
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
    if (myrank == rerank) then
      allocate(Qre_cpu(ny,nz,5),Qm_cpu(ny,5),Qm(ny,5))
    elseif (myrank == rerank+1) then
      allocate(Qre_cpu(ny,nz,5),Qm_cpu(ny,5))
    endif
    if (myrank == 0) then
      allocate(Qre(ny,nz,5),Qre_cpu(ny,nz,5))
    endif

    call cpu_time(t_start)
    do t2 = 1, np
      do t1 = 1, nt
        if (mod(myrank,2) == 0) then
          call nvtxStartRange("Send Qre", 1)
          if (myrank == rerank .and. kind(id_rescale) == 4) then
            Qre_cpu = QJ(nre1,:,:,:)
            call MPI_ISEND(Qre_cpu, 5*ny*nz, MPI_REAL8, rerank+1, 0, MPI_COMM_WORLD, ireq, ierr)
            call calc_mean(nx,ny,nz,Jacobian,QJ,Qm)
            Qm_cpu  = Qm
            call MPI_ISEND(Qm_cpu,  5*ny,    MPI_REAL8, rerank+1, 1, MPI_COMM_WORLD, ireq, ierr)
          endif
          call nvtxEndRange
          call nvtxStartRange("cal flux", 2)
          call calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,xix,etay,zetaz,Jacobian,QJ,E,F,G)
          call calc_step(nx,ny,nz,0.5d0,1.d0,xix,etay,zetaz,E,F,G,QJ,QJs,Rs) ! QJs = Q2
          call nvtxEndRange
          call nvtxStartRange("Recv Qre", 3)
          if (myrank == 0 .and. kind(id_rescale) == 4) then
            call MPI_RECV(Qre_cpu, 5*ny*nz, MPI_REAL8, rerank+1, 2, MPI_COMM_WORLD, istat, ierr)
            Qre = Qre_cpu
          endif
          call nvtxEndRange
          call nvtxStartRange("set bc", 4)
          call exchange(myrank,nranks,overlap,nx,ny,nz,QJs,Jacobian_cpu)
          call set_bc(myrank,nx,ny,nz,Jacobian,QJs,Qre)
          call nvtxEndRange
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call nvtxStartRange("recv", 5)
          call MPI_IRECV(Qre_cpu, 5*ny*nz, MPI_REAL8, rerank, 0, MPI_COMM_WORLD, ireqs(1), ierr)
          call MPI_IRECV(Qm_cpu,  5*ny,    MPI_REAL8, rerank, 1, MPI_COMM_WORLD, ireqs(2), ierr)
          call MPI_WAITALL(2,ireqs, istats, ierr)
          call nvtxEndRange
          call nvtxStartRange("rescale", 6)
          call set_rescale(t2,nx,ny,nz,y,Jacobian_cpu,Qm_cpu,Qre_cpu)
          call nvtxEndRange
          call MPI_SEND(Qre_cpu, 5*ny*nz, MPI_REAL8, 0, 2, MPI_COMM_WORLD, ierr)
        endif

        if (mod(myrank,2) == 0) then
          if (myrank == rerank .and. kind(id_rescale) == 4) then
            Qre_cpu = QJs(nre1,:,:,:)
            call MPI_ISEND(Qre_cpu, 5*ny*nz, MPI_REAL8, rerank+1, 3, MPI_COMM_WORLD, ireq, ierr)
            call calc_mean(nx,ny,nz,Jacobian,QJs,Qm)
            Qm_cpu  = Qm
            call MPI_ISEND(Qm_cpu,  5*ny,    MPI_REAL8, rerank+1, 4, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,xix,etay,zetaz,Jacobian,QJs,E,F,G)
          call calc_step(nx,ny,nz,0.5d0,2.d0,xix,etay,zetaz,E,F,G,QJ,QJs,Rs) ! QJs = Q3
          if (myrank == 0 .and. kind(id_rescale) == 4) then
            call MPI_RECV(Qre_cpu, 5*ny*nz, MPI_REAL8, rerank+1, 5, MPI_COMM_WORLD, istat, ierr)
            Qre = Qre_cpu
          endif
          call exchange(myrank,nranks,overlap,nx,ny,nz,QJs,Jacobian_cpu)
          call set_bc(myrank,nx,ny,nz,Jacobian,QJs,Qre)
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call MPI_IRECV(Qre_cpu, 5*ny*nz, MPI_REAL8, rerank, 3, MPI_COMM_WORLD, ireqs(1), ierr)
          call MPI_IRECV(Qm_cpu,  5*ny,    MPI_REAL8, rerank, 4, MPI_COMM_WORLD, ireqs(2), ierr)
          call MPI_WAITALL(2,ireqs, istats, ierr)
          call set_rescale(t2,nx,ny,nz,y,Jacobian_cpu,Qm_cpu,Qre_cpu)
          call MPI_SEND(Qre_cpu, 5*ny*nz, MPI_REAL8, 0, 5, MPI_COMM_WORLD, ierr)
        endif

        if (mod(myrank,2) == 0) then
          if (myrank == rerank .and. kind(id_rescale) == 4) then
            Qre_cpu = QJs(nre1,:,:,:)
            call MPI_ISEND(Qre_cpu, 5*ny*nz, MPI_REAL8, rerank+1, 6, MPI_COMM_WORLD, ireq, ierr)
            call calc_mean(nx,ny,nz,Jacobian,QJs,Qm)
            Qm_cpu  = Qm
            call MPI_ISEND(Qm_cpu,  5*ny,    MPI_REAL8, rerank+1, 7, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,xix,etay,zetaz,Jacobian,QJs,E,F,G)
          call calc_step(nx,ny,nz,1.d0,2.d0,xix,etay,zetaz,E,F,G,QJ,QJs,Rs) ! QJs = Q4
          if (myrank == 0 .and. kind(id_rescale) == 4) then
            call MPI_RECV(Qre_cpu, 5*ny*nz, MPI_REAL8, rerank+1, 8, MPI_COMM_WORLD, istat, ierr)
            Qre = Qre_cpu
          endif
          call exchange(myrank,nranks,overlap,nx,ny,nz,QJs,Jacobian_cpu)
          call set_bc(myrank,nx,ny,nz,Jacobian,QJs,Qre)
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call MPI_IRECV(Qre_cpu, 5*ny*nz, MPI_REAL8, rerank, 6, MPI_COMM_WORLD, ireqs(1), ierr)
          call MPI_IRECV(Qm_cpu,  5*ny,    MPI_REAL8, rerank, 7, MPI_COMM_WORLD, ireqs(2), ierr)
          call MPI_WAITAll(2, ireqs, istats, ierr)
          call set_rescale(t2,nx,ny,nz,y,Jacobian_cpu,Qm_cpu,Qre_cpu)
          call MPI_SEND(Qre_cpu, 5*ny*nz, MPI_REAL8, 0, 8, MPI_COMM_WORLD, ierr)
        endif

        if (mod(myrank,2) == 0) then
          if (myrank == rerank .and. kind(id_rescale) == 4) then
            Qre_cpu = QJs(nre1,:,:,:)
            call MPI_ISEND(Qre_cpu, 5*ny*nz, MPI_REAL8, rerank+1, 9, MPI_COMM_WORLD, ireq, ierr)
            call calc_mean(nx,ny,nz,Jacobian,QJs,Qm)
            Qm_cpu  = Qm
            call MPI_ISEND(Qm_cpu,  5*ny,    MPI_REAL8, rerank+1,10, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,xix,etay,zetaz,Jacobian,QJs,E,F,G)
          call calc_step4(nx,ny,nz,xix,etay,zetaz,E,F,G,Rs,QJ)
          if (myrank == 0 .and. kind(id_rescale) == 4) then
            call MPI_RECV(Qre_cpu, 5*ny*nz, MPI_REAL8, rerank+1,11, MPI_COMM_WORLD, istat, ierr)
            Qre = Qre_cpu
          endif
          call exchange(myrank,nranks,overlap,nx,ny,nz,QJ,Jacobian_cpu)
          call set_bc(myrank,nx,ny,nz,Jacobian,QJ,Qre)
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call MPI_IRECV(Qre_cpu, 5*ny*nz, MPI_REAL8, rerank, 9, MPI_COMM_WORLD, ireqs(1), ierr)
          call MPI_IRECV(Qm_cpu,  5*ny,    MPI_REAL8, rerank,10, MPI_COMM_WORLD, ireqs(2), ierr)
          call MPI_WAITALL(2, ireqs, istats, ierr)
          call set_rescale(t2,nx,ny,nz,y,Jacobian_cpu,Qm_cpu,Qre_cpu)
          call MPI_SEND(Qre_cpu, 5*ny*nz, MPI_REAL8, 0,11, MPI_COMM_WORLD, ierr)
        endif
      enddo

      ! send and recv device arrays
      if (mod(myrank,2) == 0) then
        Q = QJ
        call MPI_SEND(Q, nx*ny*nz*5, MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr) 
      else
        call MPI_RECV(Q, nx*ny*nz*5, MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
        call print_vtk(t2,nx,ny,nz,x,y,z,Jacobian_cpu,Q,mass0,ke0,entropy0,myrank)
      endif
    enddo
    call cpu_time(t_end)
    print *, "calculation time:", t_end - t_start

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (mod(myrank,2) == 0) then
      deallocate(QJ,QJs,Rs,E,F,G,xix,etay,zetaz,Jacobian)
    endif
    if (myrank == rerank) then
      deallocate(Qre_cpu,Qm_cpu,Qm)
    elseif (myrank == rerank+1) then
      deallocate(Qre_cpu,Qm_cpu)
    endif
    if (myrank == 0) then
      deallocate(Qre,Qre_cpu)
    endif
  end subroutine RungeKutta
end module calc_time_dev2

