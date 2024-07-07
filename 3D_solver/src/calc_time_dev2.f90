module calc_time_dev2
  use cudafor
  use mpi
  use nvtx
  use mod_globals, only : accuracy, id_hybrid, id_muscl, id_visc, id_turbulence, nt, np, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG
  use calc_physical_quantities
  use calc_steps
  use calc_hybrid
  use calc_flux
  use calc_flux_hybrid
  use calc_visc
  use calc_les
  use calc_rescale
  use set
  use print
  implicit none
contains
  subroutine calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,dx,dy,dz,Jacobian,QJ,mut,E,F,G)
    integer, intent(in), value                          :: nx, ny, nz
    type(dim3), intent(in)                              :: blocksE, blocksF, blocksG, blocks
    type(dim3), intent(in)                              :: threadsE, threadsF, threadsG, threads
    real(8), intent(in), dimension(nx-1), device        :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1), device        :: dy ! 1 / dy
    real(8), intent(in), dimension(nz-1), device        :: dz ! 1 / dz
    real(8), intent(in), dimension(nx,ny,nz), device    :: Jacobian
    real(8), intent(in), dimension(nx,ny,nz,5), device  :: QJ ! Q / Jacobian
    real(8), intent(inout), dimension(nx,ny,nz), device :: mut
    real(8), intent(out), device                        :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device                        :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device                        :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p, T, fd
    integer stat
    call calc_quantities(nx,ny,nz,Jacobian,QJ,rho,u,v,w,p,T)

    if (kind(id_hybrid) == 2) then
      call calc_E<<<blocksE,threadsE,1>>>(id_muscl,nx,ny,nz,rho,u,v,w,p,E)
      call calc_F<<<blocksF,threadsF,2>>>(id_muscl,nx,ny,nz,rho,u,v,w,p,F)
      call calc_G<<<blocksG,threadsG,3>>>(id_muscl,nx,ny,nz,rho,u,v,w,p,G)
    elseif (kind(id_hybrid) == 4) then
      call calc_Ducros<<<blocks,threads>>>(nx,ny,nz,dx,dy,dz,u,v,w,rho,p,fd)
      call calc_E_hybrid<<<blocksE,threadsE,1>>>(nx,ny,nz,rho,u,v,w,p,fd,E)
      call calc_F_hybrid<<<blocksF,threadsF,2>>>(nx,ny,nz,rho,u,v,w,p,fd,F)
      call calc_G_hybrid<<<blocksG,threadsG,3>>>(nx,ny,nz,rho,u,v,w,p,fd,G)
    endif
  
    if (id_turbulence /= 0) then
      call calc_mut<<<blocks,threads,4>>>(nx,ny,nz,rho,u,v,w,mut)
      stat = cudaDeviceSynchronize()
      call set_bc_mut(nx,ny,nz,mut)
    endif

    stat = cudaDeviceSynchronize()
    if (id_visc == 1 .or. id_turbulence /= 0) then
      call calc_Ev<<<blocksE,threadsE,1>>>(nx,ny,nz,dx,dy,dz,rho,u,v,w,T,p,mut,E)
      call calc_Fv<<<blocksF,threadsF,2>>>(nx,ny,nz,dy,dx,dz,rho,u,v,w,T,p,mut,F)
      call calc_Gv<<<blocksG,threadsG,3>>>(nx,ny,nz,dx,dy,dz,rho,u,v,w,T,p,mut,G)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG

  subroutine RungeKutta(myrank,nx,ny,nz,x,dx_cpu,xix_cpu,y,dy_cpu,etay_cpu,z,dz_cpu,zetaz_cpu,Jacobian_cpu,Q)
    integer, intent(in)    :: myrank, nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1), xix_cpu(nx-1)
    real(8), intent(in)    :: y(ny), dy_cpu(ny-1), etay_cpu(ny-1)
    real(8), intent(in)    :: z(nz), dz_cpu(nz-1), zetaz_cpu(nz-1), Jacobian_cpu(nx,ny,nz)
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer i, j, k, nre, t1, t2, itr, ilen, ierr, stat, request, status(MPI_STATUS_SIZE)
    real(8), allocatable   :: Qre_CPU(:,:,:,:), Um(:,:), Vm(:,:), Wm(:,:), pm(:,:), Tm(:,:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)         :: prop
    type(dim3)                   :: blocksE, blocksF, blocksG, blocks
    type(dim3)                   :: threadsE, threadsF, threadsG, threads
    real(8), allocatable, device :: QJ(:,:,:,:), QJ2(:,:,:,:), QJ3(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: dx(:), xix(:), dy(:), etay(:), dz(:), zetaz(:), Jacobian(:,:,:), mut(:,:,:)
    ! for plot
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0
    ! for rescaling
    real(8), allocatable, device :: Qre(:,:,:,:)
    ! rescaling plane
    nre = nx-4

    ! check GPU
    stat = cudaSetDevice(0)
    stat = cudaGetDeviceProperties(prop,0)
    ilen = verify(prop%name, ' ', .true.)
    print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", 0, ") is available"

    allocate(QJ(nx,ny,nz,5),Qre_cpu(2,ny,nz,5))
    if (mod(myrank,2) == 0) then
      allocate(QJ2(nx,ny,nz,5),QJ3(nx,ny,nz,5),E(nx-1,ny-2,nz-2,5),F(nx-2,ny-1,nz-2,5),G(nx-2,ny-2,nz-1,5))
      allocate(dx(nx-1),xix(nx-1),dy(ny-1),etay(ny-1),dz(nz-1),zetaz(nz-1),Jacobian(nx,ny,nz),mut(nx,ny,nz),Qre(2,ny,nz,5))
      call set_blocks_threads(myrank,nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads)

      ! set Q / Jacobian
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            Q(i,j,k,:) = Q(i,j,k,:) / Jacobian_cpu(i,j,k)
      enddo;enddo;enddo
  
      ! print initial condition
      call print_vtk(0,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),ke0,entropy0,myrank+1)

      ! copy on GPU
      QJ    = Q
      mut   = 0.d0
      xix   = xix_cpu
      etay  = etay_cpu
      zetaz = zetaz_cpu
      dx    = dx_cpu
      dy    = dy_cpu
      dz    = dz_cpu
      Jacobian = Jacobian_cpu
    else
      allocate(Um(2,ny),Vm(2,ny),Wm(2,ny),pm(2,ny),Tm(2,ny))
    endif

    ! share necessary data
    if (mod(myrank,2) == 0) then
      call MPI_SEND(ke0,      1, MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
      call MPI_SEND(entropy0, 1, MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
    else
      call MPI_RECV(ke0,      1, MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, status, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, status, ierr)
    endif

    do t2 = 1, np
      do t1 = 1, nt
        if (mod(myrank,2) == 0) then
          Qre_cpu(:,:,:,:) = QJ(nre+1:nre+2,:,:,:)
          call MPI_SEND(Qre_cpu, 10*ny*nz, MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
          call nvtxStartRange("calc 1step",1)
          call nvtxStartRange("calc flux",2)
          call calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,xix,etay,dz,Jacobian,QJ,mut,E,F,G)
          call nvtxEndRange
          call nvtxStartRange("calc time dev",3)
          call calc_step(nx,ny,nz,1.d0,0.d0,dx,dy,dz,E,F,G,QJ,QJ2)
          call nvtxEndRange
          call MPI_RECV(Qre_cpu, 10*ny*nz, MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, status, ierr)
          Qre = Qre_cpu
          if (myrank == 0) then 
            call set_bc1(nx,ny,nz,Jacobian,Qre,QJ2)
          elseif (myrank == 2) then
            call set_bc2(nx,ny,nz,Jacobian,Qre,QJ2)
          endif
        else
          call MPI_RECV(Qre_cpu, 10*ny*nz, MPI_REAL8, myrank-1, myrank, MPI_COMM_WORLD, status, ierr)
          call set_rescale(t1+(t2-1)*nt,myrank,nx,ny,nz,nre,2.d-3,y,Jacobian_cpu,Um,Vm,Wm,pm,Tm,Qre_cpu)
          call MPI_SEND(Qre_cpu, 10*ny*nz, MPI_REAL8, myrank-1, myrank, MPI_COMM_WORLD, ierr)
        endif

        if (mod(myrank,2) == 0) then
          Qre_cpu(:,:,:,:) = QJ(nre+1:nre+2,:,:,:)
          call MPI_SEND(Qre_cpu, 10*ny*nz, MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
          call calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,xix,etay,dz,Jacobian,QJ2,mut,E,F,G)
          call calc_step2(nx,ny,nz,0.75d0,0.25d0,0.25d0,1.d0,dx,dy,dz,E,F,G,QJ,QJ2,QJ3)
          call MPI_RECV(Qre_cpu, 10*ny*nz, MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, status, ierr)
          Qre = Qre_cpu
          if (myrank == 0) then 
            call set_bc1(nx,ny,nz,Jacobian,Qre,QJ3)
          elseif (myrank == 2) then
            call set_bc2(nx,ny,nz,Jacobian,Qre,QJ3)
          endif
        else
          call MPI_RECV(Qre_cpu, 10*ny*nz, MPI_REAL8, myrank-1, myrank, MPI_COMM_WORLD, status, ierr)
          call set_rescale(t1+(t2-1)*nt,myrank,nx,ny,nz,nre,2.d-3,y,Jacobian_cpu,Um,Vm,Wm,pm,Tm,Qre_cpu)
          call MPI_SEND(Qre_cpu, 10*ny*nz, MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, ierr)
        endif

        if (mod(myrank,2) == 0) then
          Qre_cpu(:,:,:,:) = QJ(nre+1:nre+2,:,:,:)
          call MPI_SEND(Qre_cpu, 10*ny*nz, MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
          call calc_EFG(nx,ny,nz,blocksE,blocksF,blocksG,blocks,threadsE,threadsF,threadsG,threads,xix,etay,dz,Jacobian,QJ3,mut,E,F,G)
          call calc_step3(nx,ny,nz,dx,dy,dz,E,F,G,QJ3,QJ)
          call MPI_RECV(Qre_cpu, 10*ny*nz, MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, status, ierr)
          Qre = Qre_cpu
          if (myrank == 0) then 
            call set_bc1(nx,ny,nz,Jacobian,Qre,QJ)
          elseif (myrank == 2) then
            call set_bc2(nx,ny,nz,Jacobian,Qre,QJ)
          endif
          call nvtxEndRange
        else
          call MPI_RECV(Qre_cpu, 10*ny*nz, MPI_REAL8, myrank-1, myrank, MPI_COMM_WORLD, status, ierr)
          call set_rescale(t1+(t2-1)*nt,myrank,nx,ny,nz,nre,2.d-3,y,Jacobian_cpu,Um,Vm,Wm,pm,Tm,Qre_cpu)
          call MPI_SEND(Qre_cpu, 10*ny*nz, MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, ierr)
        endif
      enddo

      ! send and recv device arrays
      if (mod(myrank,2) == 0) then
        Q = QJ
        call nvtxStartRange("MPI_SEND",4)
        call MPI_SEND(Q, nx*ny*nz*5, MPI_REAL8, myrank+1, myrank+1, MPI_COMM_WORLD, ierr) 
        call nvtxEndRange
      else
        call nvtxStartRange("MPI_RECV",5)
        call MPI_RECV(Q, nx*ny*nz*5, MPI_REAL8, myrank-1, myrank,   MPI_COMM_WORLD, status, ierr)
        call nvtxEndRange
        call nvtxStartRange("print",6)
        call print_vtk(t2,nx,ny,nz,real(x),real(y),real(z),real(Jacobian_cpu),real(Q),ke0,entropy0,myrank)
        call nvtxEndRange
      endif
    enddo
    
    if (mod(myrank,2) == 0) then
      deallocate(QJ2,QJ3,E,F,G,dx,xix,dy,etay,dz,zetaz,Jacobian,mut,Qre)
    else
      deallocate(Um,Vm,Wm,pm,Tm)
    endif
    deallocate(QJ,Qre_cpu)
  end subroutine RungeKutta
end module calc_time_dev2

