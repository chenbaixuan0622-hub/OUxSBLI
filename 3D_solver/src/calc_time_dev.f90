module calc_time_dev
  use cudafor
  use mpi
  use nvtx
  use mod_globals, only : accuracy, id_scheme, id_turbulence, id_rescale, nt, np, nre, &
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
  interface calc_EFG
    module procedure calc_EFG_Euler, calc_EFG_visc, calc_EFG_LES
  end interface calc_EFG

  interface RungeKutta
    module procedure RungeKutta_3rd, RungeKutta_4th
  end interface
contains
  subroutine calc_EFG_Euler(id_visc,nx,ny,nz,dx,dy,dz,Jacobian,QJ,E,F,G)
    integer(kind=2), intent(in), value :: id_visc
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(in), device        :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device        :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device        :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device        :: Jacobian(ny)
    real(8), intent(in), device        :: QJ(nx,ny,nz,5) ! Q / Jacobian
    real(8), intent(out), device       :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device       :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device       :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    integer stat
    call calc_E<<<blocksE,threadsE,1>>>(nx,ny,nz,Jacobian,QJ,E)
    call calc_F<<<blocksF,threadsF,2>>>(nx,ny,nz,Jacobian,QJ,F)
    call calc_G<<<blocksG,threadsG,3>>>(nx,ny,nz,Jacobian,QJ,G)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_Euler

  subroutine calc_EFG_visc(id_visc,nx,ny,nz,dx,dy,dz,Jacobian,QJ,E,F,G)
    integer(kind=4), intent(in), value :: id_visc
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(in), device        :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device        :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device        :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device        :: Jacobian(ny)
    real(8), intent(in), device        :: QJ(nx,ny,nz,5) ! Q / Jacobian
    real(8), intent(out), device       :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device       :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device       :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    integer stat
    call calc_E<<<blocksE,threadsE,1>>>(nx,ny,nz,Jacobian,QJ,E)
    call calc_F<<<blocksF,threadsF,2>>>(nx,ny,nz,Jacobian,QJ,F)
    call calc_G<<<blocksG,threadsG,3>>>(nx,ny,nz,Jacobian,QJ,G)

    stat = cudaDeviceSynchronize()
    call calc_Ev<<<blocksE,threadsE,1>>>(nx,ny,nz,dx,dy,dz,Jacobian,QJ,E)
    call calc_Fv<<<blocksF,threadsF,2>>>(nx,ny,nz,dy,dx,dz,Jacobian,QJ,F)
    call calc_Gv<<<blocksG,threadsG,3>>>(nx,ny,nz,dx,dy,dz,Jacobian,QJ,G)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_visc
  
  subroutine calc_EFG_LES(id_visc,nx,ny,nz,dx,dy,dz,Jacobian,QJ,E,F,G)
    integer(kind=8), intent(in), value :: id_visc
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(in), device        :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device        :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device        :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device        :: Jacobian(ny)
    real(8), intent(in), device        :: QJ(nx,ny,nz,5) ! Q / Jacobian
    real(8), intent(out), device       :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device       :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device       :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), allocatable, device       :: mut(:,:,:), qc2(:,:,:)
    integer stat
    allocate(mut(nx,ny,nz), qc2(nx,ny,nz))
    mut = 0.d0
    qc2 = 0.d0
    
    call calc_E<<<blocksE,threadsE,1>>>(nx,ny,nz,Jacobian,QJ,E)
    call calc_F<<<blocksF,threadsF,2>>>(nx,ny,nz,Jacobian,QJ,F)
    call calc_G<<<blocksG,threadsG,3>>>(nx,ny,nz,Jacobian,QJ,G)
    call calc_mut<<<blocks,threads,4>>>(nx,ny,nz,dx,dy,dz,Jacobian,QJ,mut,qc2)
    stat = cudaDeviceSynchronize()
    call set_bc_mut(nx,ny,nz,mut,qc2)

    call calc_Ev_LES<<<blocksE,threadsE,1>>>(nx,ny,nz,dx,dy,dz,Jacobian,QJ,mut,qc2,E)
    call calc_Fv_LES<<<blocksF,threadsF,2>>>(nx,ny,nz,dy,dx,dz,Jacobian,QJ,mut,qc2,F)
    call calc_Gv_LES<<<blocksG,threadsG,3>>>(nx,ny,nz,dx,dy,dz,Jacobian,QJ,mut,qc2,G)
    stat = cudaDeviceSynchronize()
    deallocate(mut, qc2)
  end subroutine calc_EFG_LES

  subroutine RungeKutta_3rd(id_RungeKutta,myrank,nx,ny,nz,x,dx_cpu,y,dy_cpu,z,dz_cpu,Jacobian_cpu,Q)
    use mod_globals, only : id_visc
    integer(kind=2), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z(nz), dz_cpu(nz-1), Jacobian_cpu(ny)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    integer i, j, k, t1, t2, itr, ilen, ierr, stat, request, status(MPI_STATUS_SIZE)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)         :: prop
    real(8), allocatable, device :: QJ(:,:,:,:), QJ2(:,:,:,:), QJ3(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), zetaz(:), Jacobian(:)
    ! for plot
    real(8) :: mass0 = 1.d0, ke0 = 1.d0, entropy0 = 1.d0

    ! check GPU
    stat = cudaSetDevice(0)
    stat = cudaGetDeviceProperties(prop,0)
    ilen = verify(prop%name, ' ', .true.)
    print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", 0, ") is available"

    if (myrank == 0) then
      allocate(QJ(nx,ny,nz,5),QJ2(nx,ny,nz,5),QJ3(nx,ny,nz,5),E(nx-1,ny-2,nz-2,5),F(nx-2,ny-1,nz-2,5),G(nx-2,ny-2,nz-1,5))
      allocate(xix(nx-1),etay(ny-1),zetaz(nz-1),Jacobian(ny))

      ! set Q / Jacobian
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            Q(i,j,k,:) = Q(i,j,k,:) / Jacobian_cpu(j)
      enddo;enddo;enddo

      ! print initial condition
      call print_vtk(0,nx,ny,nz,x,y,z,Jacobian_cpu,Q,mass0,ke0,entropy0)

      ! copy on GPU
      QJ       = Q
      xix      = 1.d0 / dx_cpu
      etay     = 1.d0 / dy_cpu
      zetaz    = 1.d0 / dz_cpu
      Jacobian = Jacobian_cpu
    endif

    ! share necessary data
    call MPI_BCAST(mass0,    1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ke0,      1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(entropy0, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    do t2 = 1, np
      if (myrank == 0) then
        do t1 = 1, nt
          call nvtxStartRange("calc 1step",1)
          call nvtxStartRange("calc flux",2)
          call calc_EFG(id_visc,nx,ny,nz,xix,etay,zetaz,Jacobian,QJ,E,F,G)
          call nvtxEndRange
          call nvtxStartRange("calc time dev",3)
          call calc_step(nx,ny,nz,1.d0,0.d0,xix,etay,zetaz,E,F,G,QJ,QJ2)
          call nvtxEndRange
          call set_bc(nx,ny,nz,Jacobian,QJ2)

          call calc_EFG(id_visc,nx,ny,nz,xix,etay,zetaz,Jacobian,QJ2,E,F,G)
          call calc_step2(nx,ny,nz,0.75d0,0.25d0,0.25d0,1.d0,xix,etay,zetaz,E,F,G,QJ,QJ2,QJ3)
          call set_bc(nx,ny,nz,Jacobian,QJ3)

          call calc_EFG(id_visc,nx,ny,nz,xix,etay,zetaz,Jacobian,QJ3,E,F,G)
          call calc_step3(nx,ny,nz,xix,etay,zetaz,E,F,G,QJ3,QJ)
          call set_bc(nx,ny,nz,Jacobian,QJ)
          call nvtxEndRange
        enddo
      endif

      ! send and recv device arrays
      if (myrank == 0) then
        Q = QJ
        call nvtxStartRange("MPI_SEND",4)
        call MPI_SEND(Q, nx*ny*nz*5, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
        call nvtxEndRange
      elseif (myrank == 1) then
        call nvtxStartRange("MPI_RECV",5)
        call MPI_RECV(Q, nx*ny*nz*5, MPI_REAL8, 0, 0, MPI_COMM_WORLD, status, ierr)
        call nvtxEndRange
        call nvtxStartRange("print",6)
        call print_vtk(t2,nx,ny,nz,x,y,z,Jacobian_cpu,Q,mass0,ke0,entropy0)
        call nvtxEndRange
      endif
    enddo
    
    if (myrank == 0) then
      deallocate(QJ,QJ2,QJ3,E,F,G,xix,etay,zetaz,Jacobian)
    endif
  end subroutine RungeKutta_3rd

  subroutine RungeKutta_4th(id_RungeKutta,myrank,nx,ny,nz,x,dx_cpu,y,dy_cpu,z,dz_cpu,Jacobian_cpu,Q)
    use mod_globals, only : id_visc
    integer(kind=4), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z(nz), dz_cpu(nz-1), Jacobian_cpu(ny)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    integer i, j, k, t1, t2, itr, ierr, ilen, stat, ireq, ireq2, istat(MPI_STATUS_SIZE), status(MPI_STATUS_SIZE)
    ! rescale !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable         :: Qre_cpu(:,:,:,:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)         :: prop
    real(8), allocatable, device :: QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), zetaz(:), Jacobian(:)
    ! for plot
    real(8) :: mass0 = 1.d0, ke0 = 1.d0, entropy0 = 1.d0

    ! check GPU
    stat = cudaSetDevice(0)
    stat = cudaGetDeviceProperties(prop,0)
    ilen = verify(prop%name, ' ', .true.)

    if (myrank == 0) then
      allocate(QJ(nx,ny,nz,5),QJs(nx,ny,nz,5),Rs(nx-2,ny-2,nz-2,5),E(nx-1,ny-2,nz-2,5),F(nx-2,ny-1,nz-2,5),G(nx-2,ny-2,nz-1,5))
      allocate(xix(nx-1),etay(ny-1),zetaz(nz-1),Jacobian(ny))

      ! set Q / Jacobian
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            Q(i,j,k,:) = Q(i,j,k,:) / Jacobian_cpu(j)
      enddo;enddo;enddo

      ! print initial condition
      call print_vtk(0,nx,ny,nz,x,y,z,Jacobian_cpu,Q,mass0,ke0,entropy0)

      ! copy on GPU
      QJ       = Q
      Rs       = 0.d0
      xix      = 1.d0 / dx_cpu
      etay     = 1.d0 / dy_cpu
      zetaz    = 1.d0 / dz_cpu
      Jacobian = Jacobian_cpu
    endif

    ! share necessary data
    call MPI_BCAST(mass0,    1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ke0,      1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(entropy0, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    ! rescale
    allocate(Qre_cpu(10,ny,nz,5))

    do t2 = 1, np
      do t1 = 1, nt
        if (myrank == 0) then
          call nvtxStartRange("Send Qre", 1)
          if (kind(id_rescale) == 4) then
            Qre_cpu(:,:,:,:) = QJ(nre:nre+45:5,:,:,:)
            call MPI_ISEND(Qre_cpu, 50*ny*nz, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ireq, ierr)
          endif
          call nvtxEndRange
          call nvtxStartRange("calc flux", 2)
          call calc_EFG(id_visc,nx,ny,nz,xix,etay,zetaz,Jacobian,QJ,E,F,G)
          call calc_step(nx,ny,nz,0.5d0,1.d0,xix,etay,zetaz,E,F,G,QJ,QJs,Rs) ! QJs = Q2
          call nvtxEndRange
          call nvtxStartRange("Recv Qre", 3)
          if (kind(id_rescale) == 4) then
            call MPI_RECV(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 1, 1, MPI_COMM_WORLD, istat, ierr)
            QJs(1,:,:,:) = Qre_cpu(1,:,:,:)
          endif
          call nvtxEndRange
          call nvtxStartRange("set bc", 4)
          call set_bc(nx,ny,nz,Jacobian,QJs)
          call nvtxEndRange
        elseif (myrank == 1 .and. kind(id_rescale) == 4) then
          call MPI_IRECV(Qre_cpu, 50*ny*nz, MPI_REAL8, 0, 0, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAIT(ireq, istat, ierr)
          call nvtxStartRange("rescale", 5)
          call set_rescale(t2,nx,ny,nz,nre,y,Jacobian_cpu,Qre_cpu)
          call nvtxEndRange
          call MPI_SEND(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 0, 1, MPI_COMM_WORLD, ierr)
        endif

        if (myrank == 0) then
          if (kind(id_rescale) == 4) then
            Qre_cpu(:,:,:,:) = QJs(nre:nre+45:5,:,:,:)
            call MPI_ISEND(Qre_cpu, 50*ny*nz, MPI_REAL8, 1, 2, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(id_visc,nx,ny,nz,xix,etay,zetaz,Jacobian,QJs,E,F,G)
          call calc_step(nx,ny,nz,0.5d0,2.d0,xix,etay,zetaz,E,F,G,QJ,QJs,Rs) ! QJs = Q3
          if (kind(id_rescale) == 4) then
            call MPI_RECV(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 1, 3, MPI_COMM_WORLD, istat, ierr)
            QJs(1,:,:,:) = Qre_cpu(1,:,:,:)
          endif
          call set_bc(nx,ny,nz,Jacobian,QJs)
        elseif (myrank == 1 .and. kind(id_rescale) == 4) then
          call MPI_IRECV(Qre_cpu, 50*ny*nz, MPI_REAL8, 0, 2, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAit(ireq, istat, ierr)
          call set_rescale(t2,nx,ny,nz,nre,y,Jacobian_cpu,Qre_cpu)
          call MPI_SEND(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 0, 3, MPI_COMM_WORLD, ierr)
        endif

        if (myrank == 0) then
          if (kind(id_rescale) == 4) then
            Qre_cpu(:,:,:,:) = QJs(nre:nre+45:5,:,:,:)
            call MPI_ISEND(Qre_cpu, 50*ny*nz, MPI_REAL8, 1, 4, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(id_visc,nx,ny,nz,xix,etay,zetaz,Jacobian,QJs,E,F,G)
          call calc_step(nx,ny,nz,1.0d0,2.d0,xix,etay,zetaz,E,F,G,QJ,QJs,Rs) ! QJs = Q4
          if (kind(id_rescale) == 4) then
            call MPI_RECV(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 1, 5, MPI_COMM_WORLD, istat, ierr)
            QJs(1,:,:,:) = Qre_cpu(1,:,:,:)
          endif
          call set_bc(nx,ny,nz,Jacobian,QJs)
        elseif (myrank == 1 .and. kind(id_rescale) == 4) then
          call MPI_IRECV(Qre_cpu, 50*ny*nz, MPI_REAL8, 0, 4, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAIT(ireq, istat, ierr)
          call set_rescale(t2,nx,ny,nz,nre,y,Jacobian_cpu,Qre_cpu)
          call MPI_SEND(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 0, 5, MPI_COMM_WORLD, ierr)
        endif

        if (myrank == 0) then
          if (kind(id_rescale) == 4) then
            Qre_cpu(:,:,:,:) = QJs(nre:nre+45:5,:,:,:)
            call MPI_ISEND(Qre_cpu, 50*ny*nz, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ireq, ierr)
          endif
          call calc_EFG(id_visc,nx,ny,nz,xix,etay,zetaz,Jacobian,QJs,E,F,G)
          call calc_step4(nx,ny,nz,xix,etay,zetaz,E,F,G,Rs,QJ)
          if (kind(id_rescale) == 4) then
            call MPI_RECV(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 1, 1, MPI_COMM_WORLD, istat, ierr)
            QJ(1,:,:,:) = Qre_cpu(1,:,:,:)
          endif
          call set_bc(nx,ny,nz,Jacobian,QJ)
        elseif (myrank == 1 .and. kind(id_rescale) == 4) then
          call MPI_IRECV(Qre_cpu, 50*ny*nz, MPI_REAL8, 0, 0, MPI_COMM_WORLD, ireq, ierr)
          call MPI_WAIT(ireq, istat, ierr)
          call set_rescale(t2,nx,ny,nz,nre,y,Jacobian_cpu,Qre_cpu)
          call MPI_SEND(Qre_cpu(1,:,:,:), 5*ny*nz, MPI_REAL8, 0, 1, MPI_COMM_WORLD, ierr)
        endif
      enddo

      ! send and recv device arrays
      if (myrank == 0) then
        Q = QJ
        call MPI_SEND(Q, nx*ny*nz*5, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
      elseif (myrank == 1) then
        call MPI_RECV(Q, nx*ny*nz*5, MPI_REAL8, 0, 0, MPI_COMM_WORLD, istat, ierr)
        call print_vtk(t2,nx,ny,nz,x,y,z,Jacobian_cpu,Q,mass0,ke0,entropy0)
      endif
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (myrank == 0) then
      deallocate(QJ,QJs,Rs,E,F,G,xix,etay,zetaz,Jacobian)
    endif
    deallocate(Qre_cpu)
  end subroutine RungeKutta_4th
end module calc_time_dev

