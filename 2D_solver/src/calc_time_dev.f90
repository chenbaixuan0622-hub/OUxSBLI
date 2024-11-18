module calc_time_dev
  use cudafor
  use mpi
  use mod_globals, only : accuracy, id_scheme, nt, np, &
  & blocks, threads, blocksE, blocksF, threadsE, threadsF, &
  & blocksEv, blocksFv, threadsEv, threadsFv
  use calc_physical_quantities
  use calc_steps
  use calc_hybrid
  use calc_flux
  use calc_visc
  use set
  use print
  implicit none
  interface calc_EFG
    module procedure calc_EFG_Euler, calc_EFG_visc
  end interface calc_EFG
contains
  subroutine calc_EFG_Euler(id_visc,nx,ny,dx,dy,Jacobian,QJ,E,F)
    integer(kind=2), intent(in), value :: id_visc
    integer, intent(in), value         :: nx, ny
    real(8), intent(in), device        :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device        :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device        :: Jacobian(ny)
    real(8), intent(in), device        :: QJ(nx,ny,4) ! Q / Jacobian
    real(8), intent(out), device       :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(out), device       :: F(nx-accuracy,ny-accuracy+1,4)
    real(8), dimension(nx,ny), device  :: rho, u, v, p
    integer stat
    call calc_quantities_2D(nx,ny,Jacobian,QJ,rho,u,v,p)
    call calc_E<<<blocksE,threadsE,1>>>(nx,ny,rho,u,v,p,E)
    call calc_F<<<blocksF,threadsF,2>>>(nx,ny,rho,u,v,p,F)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_Euler

  subroutine calc_EFG_visc(id_visc,nx,ny,dx,dy,Jacobian,QJ,E,F)
    integer(kind=4), intent(in), value :: id_visc
    integer, intent(in), value         :: nx, ny
    real(8), intent(in), device        :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device        :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device        :: Jacobian(ny)
    real(8), intent(in), device        :: QJ(nx,ny,4) ! Q / Jacobian
    real(8), intent(out), device       :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), intent(out), device       :: F(nx-accuracy,ny-accuracy+1,4)
    real(8), dimension(nx,ny), device  :: rho, u, v, p
    integer stat
    call calc_quantities_2D(nx,ny,Jacobian,QJ,rho,u,v,p)
    call calc_E<<<blocksE,threadsE,1>>>(nx,ny,rho,u,v,p,E)
    call calc_F<<<blocksF,threadsF,2>>>(nx,ny,rho,u,v,p,F)

    stat = cudaDeviceSynchronize()
    call calc_Ev<<<blocksEv,threadsEv,1>>>(nx,ny,dx,dy,rho,u,v,p,E)
    call calc_Fv<<<blocksFv,threadsFv,2>>>(nx,ny,dy,dx,rho,u,v,p,F)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_visc
  
  subroutine RungeKutta(id_RungeKutta,myrank,nx,ny,nz,x,dx_cpu,y,dy_cpu,z,dz_cpu,Jacobian_cpu,Q)
    use mod_globals, only : id_visc
    integer(kind=4), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z(nz), dz_cpu(nz-1), Jacobian_cpu(ny)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    integer i, j, k, t1, t2, itr, ierr, ilen, ndevices, stat, ireq, ireqs(2), istat(MPI_STATUS_SIZE), istats(MPI_STATUS_SIZE,2)
    real(8) Q2d(nx,ny,4)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type(cudaDeviceProp)         :: prop
    real(8), allocatable, device :: QJ(:,:,:), QJs(:,:,:), Rs(:,:,:), E(:,:,:), F(:,:,:)
    real(8), allocatable, device :: xix(:), etay(:), Jacobian(:)

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

    if (myrank == 0) then
      allocate(QJ(nx,ny,4),QJs(nx,ny,4),Rs(nx-2,ny-2,4),E(nx-1,ny-2,4),F(nx-2,ny-1,4))
      allocate(xix(nx-1),etay(ny-1),Jacobian(ny))

      ! set Q / Jacobian
      do j = 1, ny
        do i = 1, nx
          Q2d(i,j,1) = Q(i,j,1,1) / Jacobian_cpu(j)
          Q2d(i,j,2) = Q(i,j,1,2) / Jacobian_cpu(j)
          Q2d(i,j,3) = Q(i,j,1,3) / Jacobian_cpu(j)
          Q2d(i,j,4) = Q(i,j,1,5) / Jacobian_cpu(j)
      enddo;enddo

      ! print initial condition
      call print_vtk(0,nx,ny,x,y,Jacobian_cpu,Q2d)

      ! copy on GPU
      QJ       = Q2d
      Rs       = 0.d0
      xix      = 1.d0 / dx_cpu
      etay     = 1.d0 / dy_cpu
      Jacobian = Jacobian_cpu
    endif

    do t2 = 1, np
      do t1 = 1, nt
        if (myrank == 0) then
          call calc_EFG(id_visc,nx,ny,xix,etay,Jacobian,QJ,E,F)
          call calc_step(nx,ny,0.5d0,1.d0,xix,etay,E,F,QJ,QJs,Rs) ! QJs = Q2
          call set_bc(nx,ny,Jacobian,QJs)

          call calc_EFG(id_visc,nx,ny,xix,etay,Jacobian,QJs,E,F)
          call calc_step(nx,ny,0.5d0,2.d0,xix,etay,E,F,QJ,QJs,Rs) ! QJs = Q3
          call set_bc(nx,ny,Jacobian,QJs)

          call calc_EFG(id_visc,nx,ny,xix,etay,Jacobian,QJs,E,F)
          call calc_step(nx,ny,1.0d0,2.d0,xix,etay,E,F,QJ,QJs,Rs) ! QJs = Q4
          call set_bc(nx,ny,Jacobian,QJs)

          call calc_EFG(id_visc,nx,ny,xix,etay,Jacobian,QJs,E,F)
          call calc_step4(nx,ny,xix,etay,E,F,Rs,QJ)
          call set_bc(nx,ny,Jacobian,QJ)
        endif
      enddo

      ! send and recv device arrays
      if (myrank == 0) then
        Q2d = QJ
        call MPI_SEND(Q2d, nx*ny*4, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
      elseif (myrank == 1) then
        call MPI_RECV(Q2d, nx*ny*4, MPI_REAL8, 0, 0, MPI_COMM_WORLD, istat, ierr)
        call print_vtk(t2,nx,ny,x,y,Jacobian_cpu,Q2d)
      endif
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (myrank == 0) then
      deallocate(QJ,QJs,Rs,E,F,xix,etay,Jacobian)
    endif
  end subroutine RungeKutta
end module calc_time_dev

