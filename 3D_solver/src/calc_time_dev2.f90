module calc_time_dev2
  use cudafor
  use mpi
  use nvtx
  use mod_globals, only : accuracy, id_rescale, nt, np, nre, rerank
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
    integer stat
    call calc_E<<<blocksE,threadsE,1>>>(nx,ny,nz,Jacobian,QJ,E)
    call calc_F<<<blocksF,threadsF,2>>>(nx,ny,nz,Jacobian,QJ,F)
    call calc_G<<<blocksG,threadsG,3>>>(nx,ny,nz,Jacobian,QJ,G)
  
    stat = cudaDeviceSynchronize()
    call calc_Ev<<<blocksE,threadsE,1>>>(nx,ny,nz,dx,dy,dz,Jacobian,QJ,E)
    call calc_Fv<<<blocksF,threadsF,2>>>(nx,ny,nz,dy,dx,dz,Jacobian,QJ,F)
    call calc_Gv<<<blocksG,threadsG,3>>>(nx,ny,nz,dx,dy,dz,Jacobian,QJ,G)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG

  subroutine RungeKutta(myrank,nx,ny,nz,x,dx_cpu,y,dy_cpu,z,dz_cpu,Jacobian_cpu,Q)
    integer, intent(in)    :: myrank, nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)    :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)    :: z(nz), dz_cpu(nz-1), Jacobian_cpu(ny)
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer i, j, k, t1, t2, itr, ilen, ierr, stat, ndevices, ireq, ireq2, istat(MPI_STATUS_SIZE)
    ! rescale
    real(8), allocatable         :: Qre_cpu(:,:,:,:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)         :: prop
    type(dim3)                   :: blocksE, blocksF, blocksG, blocks
    type(dim3)                   :: threadsE, threadsF, threadsG, threads
    real(8), allocatable, device :: QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), zetaz(:), Jacobian(:)
    ! for plot
    real(8) :: mass0 = 1.d0, ke0 = 1.d0, entropy0 = 1.d0

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

    if (mod(myrank,2) == 0) then
      allocate(QJ(nx,ny,nz,5),QJs(nx,ny,nz,5),Rs(nx-2,ny-2,nz-2,5),E(nx-1,ny-2,nz-2,5),F(nx-2,ny-1,nz-2,5),G(nx-2,ny-2,nz-1,5))
      allocate(xix(nx-1),etay(ny-1),zetaz(nz-1),Jacobian(ny))

      print *, "myrank is ", myrank, " memory allocation has completed"
      call set_blocks_threads(myrank,nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads)

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
    if (myrank == 0 .or. myrank == rerank .or. myrank == rerank+1) then
      allocate(Qre_cpu(10,ny,nz,5))
    endif

    do t2 = 1, np
      do t1 = 1, nt
        if (mod(myrank,2) == 0) then
          if (myrank == rerank) then
            Qre_cpu(:,:,:,:) = QJ(nre:nre+9,:,:,:)
            call MPI_ISEND(Qre_cpu, 50*ny*nz, MPI_REAL8, rerank+1, 0, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,xix,etay,zetaz,Jacobian,QJ,E,F,G)
          call calc_step(nx,ny,nz,0.5d0,1.d0,xix,etay,zetaz,E,F,G,QJ,QJs,Rs) ! QJs = Q2
          if (myrank == 0) then
            call MPI_RECV(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, rerank+1, 1, MPI_COMM_WORLD, istat, ierr)
            QJs(1,:,:,:) = Qre_cpu(1,:,:,:)
          endif
          call set_bc(myrank,nx,ny,nz,Jacobian,QJs)
        elseif (myrank == rerank+1) then
          call MPI_IRECV(Qre_cpu, 50*ny*nz, MPI_REAL8, rerank, 0, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAIT(ireq, istat, ierr)
          call set_rescale(t2,nx,ny,nz,nre,y,Jacobian_cpu,Qre_cpu)
          call MPI_SEND(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 0, 1, MPI_COMM_WORLD, ierr)
        endif

        if (mod(myrank,2) == 0) then
          if (myrank == rerank) then
            Qre_cpu(:,:,:,:) = QJ(nre:nre+9,:,:,:)
            call MPI_ISEND(Qre_cpu, 50*ny*nz, MPI_REAL8, rerank+1, 2, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,xix,etay,zetaz,Jacobian,QJs,E,F,G)
          call calc_step(nx,ny,nz,0.5d0,2.d0,xix,etay,zetaz,E,F,G,QJ,QJs,Rs) ! QJs = Q3
          if (myrank == 0) then
            call MPI_RECV(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, rerank+1, 3, MPI_COMM_WORLD, istat, ierr)
            QJs(1,:,:,:) = Qre_cpu(1,:,:,:)
          endif
          call set_bc(myrank,nx,ny,nz,Jacobian,QJs)
        elseif (myrank == rerank+1) then
          call MPI_IRECV(Qre_cpu, 50*ny*nz, MPI_REAL8, rerank, 2, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAIT(ireq, istat, ierr)
          call set_rescale(t2,nx,ny,nz,nre,y,Jacobian_cpu,Qre_cpu)
          call MPI_SEND(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 0, 3, MPI_COMM_WORLD, ierr)
        endif

        if (mod(myrank,2) == 0) then
          if (myrank == rerank) then
            Qre_cpu(:,:,:,:) = QJ(nre:nre+9,:,:,:)
            call MPI_ISEND(Qre_cpu, 50*ny*nz, MPI_REAL8, rerank+1, 4, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,xix,etay,zetaz,Jacobian,QJs,E,F,G)
          call calc_step(nx,ny,nz,1.d0,2.d0,xix,etay,zetaz,E,F,G,QJ,QJs,Rs) ! QJs = Q4
          if (myrank == 0) then
            call MPI_RECV(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, rerank+1, 5, MPI_COMM_WORLD, istat, ierr)
            QJs(1,:,:,:) = Qre_cpu(1,:,:,:)
          endif
          call set_bc(myrank,nx,ny,nz,Jacobian,QJs)
        elseif (myrank == rerank+1) then
          call MPI_IRECV(Qre_cpu, 50*ny*nz, MPI_REAL8, rerank, 4, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAIT(ireq, istat, ierr)
          call set_rescale(t2,nx,ny,nz,nre,y,Jacobian_cpu,Qre_cpu)
          call MPI_SEND(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 0, 5, MPI_COMM_WORLD, ierr)
        endif

        if (mod(myrank,2) == 0) then
          if (myrank == rerank) then
            Qre_cpu(:,:,:,:) = QJ(nre:nre+9,:,:,:)
            call MPI_ISEND(Qre_cpu, 50*ny*nz, MPI_REAL8, rerank+1, 6, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,xix,etay,zetaz,Jacobian,QJs,E,F,G)
          call calc_step4(nx,ny,nz,xix,etay,zetaz,E,F,G,Rs,QJ)
          if (myrank == 0) then
            call MPI_RECV(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, rerank+1, 7, MPI_COMM_WORLD, istat, ierr)
            QJ(1,:,:,:) = Qre_cpu(1,:,:,:)
          endif
          call set_bc(myrank,nx,ny,nz,Jacobian,QJ)
        elseif (myrank == rerank+1) then
          call MPI_IRECV(Qre_cpu, 50*ny*nz, MPI_REAL8, rerank, 6, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAIT(ireq, istat, ierr)
          call set_rescale(t2,nx,ny,nz,nre,y,Jacobian_cpu,Qre_cpu)
          call MPI_SEND(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 0, 7, MPI_COMM_WORLD, ierr)
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
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (mod(myrank,2) == 0) then
      deallocate(QJ,QJs,Rs,E,F,G,xix,etay,zetaz,Jacobian)
    endif
    if (myrank == 0 .or. myrank == rerank .or. myrank == rerank+1) then
      deallocate(Qre_cpu)
    endif
  end subroutine RungeKutta
end module calc_time_dev2

