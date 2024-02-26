module calc_time_dev
  implicit none
contains
  attributes(global) subroutine calc_quantities(nx,ny,nz,gamma,Q,rho,u,v,w,p)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: gamma
    real(8), intent(in), dimension(nx,ny,nz,5), device :: Q
    real(8), intent(out), dimension(nx,ny,nz), device :: rho, u, v, w, p
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    rho(i,j,k) = Q(i,j,k,1)
    u(i,j,k) = Q(i,j,k,2) / rho(i,j,k)
    v(i,j,k) = Q(i,j,k,3) / rho(i,j,k)
    w(i,j,k) = Q(i,j,k,4) / rho(i,j,k)
    p(i,j,k) = (gamma-1.d0)*(Q(i,j,k,5)-0.5d0*rho(i,j,k)*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2))
  end subroutine

  subroutine RungeKutta(nx,ny,nz,nt,np,dx,dy,dz,dt,gamma,mu,kappa,Cp,Q)
    use iso_fortran_env
    use cudafor
    use mod_globals, only : accuracy
    use calc_steps
    use calc_flux
    use set_bc
    use print
    integer, intent(in) :: nx, ny, nz, nt, np
    real(8), intent(in) :: dx, dy, dz, dt, gamma, mu, kappa, Cp
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer t1, t2
    integer(kind=2**(accuracy/2)) :: id
    real(8) dxi, dyi, dzi
    ! GPU
    integer stat, len
    type(cudaDeviceProp) :: prop
    type(dim3) :: blocks, threads
    type(dim3) :: blocks_offset, threads_offset
    type(dim3) :: blocksE, blocksF, blocksG, threadsE, threadsF, threadsG
    real(8), dimension(nx,ny,nz,5), device :: Q_d, Q2, Q3
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p
    real(8), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    dxi = 1.0d0 / dx
    dyi = 1.0d0 / dy
    dzi = 1.0d0 / dz

    ! check active device
    print *,"\nChecking for GPU"
    stat = cudaSetDevice(0)
    stat = cudaGetDeviceProperties(prop,0)
    len = verify(prop%name, ' ', .true.)
    print '(1x, a, a, i1,a)', prop%name(1:len), " (GPU) is available"

    ! thread num must be less than 1024
    blocks = dim3(nx/10,ny/10,nz/10)
    threads = dim3(10,10,10)
    blocks_offset = dim3((nx-accuracy)/7,(ny-accuracy)/7,(nz-accuracy)/7)
    threads_offset = dim3(7,7,7)
    blocksE = dim3((nx-accuracy+1)/3,(ny-accuracy)/7,(nz-accuracy)/7)
    blocksF = dim3((nx-accuracy)/7,(ny-accuracy+1)/3,(nz-accuracy)/7)
    blocksG = dim3((nx-accuracy)/7,(ny-accuracy)/7,(nz-accuracy+1)/3)
    threadsE = dim3(3,7,7)
    threadsF = dim3(7,3,7)
    threadsG = dim3(7,7,3)

    Q_d = Q
    do t2 = 1, np
      do t1 = 1, nt
        call calc_quantities<<<blocks,threads>>>(nx,ny,nz,gamma,Q_d,rho,u,v,w,p)
        !print *, trim(cudaGetErrorString(cudaGetLastError()))
        stat = cudaDeviceSynchronize()
        call calc_E<<<blocksE,threadsE>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,E)
        !print *, trim(cudaGetErrorString(cudaGetLastError()))
        call calc_F<<<blocksF,threadsF>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,F)
        !print *, trim(cudaGetErrorString(cudaGetLastError()))
        call calc_G<<<blocksG,threadsG>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,G)
        !print *, trim(cudaGetErrorString(cudaGetLastError()))
        stat = cudaDeviceSynchronize()
        call calc_step1<<<blocks_offset,threads_offset>>>(nx,ny,nz,dxi,dyi,dzi,dt,E,F,G,Q_d,Q2)
        !print *, trim(cudaGetErrorString(cudaGetLastError()))
        stat = cudaDeviceSynchronize()
        call set_cyclic_bc_d(id,nx,ny,nz,Q2)
        stat = cudaDeviceSynchronize()

        call calc_quantities<<<blocks,threads>>>(nx,ny,nz,gamma,Q2,rho,u,v,w,p)
        stat = cudaDeviceSynchronize()
        call calc_E<<<blocksE,threadsE>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,E)
        call calc_F<<<blocksF,threadsF>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,F)
        call calc_G<<<blocksG,threadsG>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,G)
        stat = cudaDeviceSynchronize()
        call calc_step2<<<blocks_offset,threads_offset>>>(nx,ny,nz,dxi,dyi,dzi,dt,E,F,G,Q_d,Q2,Q3)
        stat = cudaDeviceSynchronize()
        call set_cyclic_bc_d(id,nx,ny,nz,Q3)
        stat = cudaDeviceSynchronize()

        call calc_quantities<<<blocks,threads>>>(nx,ny,nz,gamma,Q3,rho,u,v,w,p)
        stat = cudaDeviceSynchronize()
        call calc_E<<<blocksE,threadsE>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,E)
        call calc_F<<<blocksF,threadsF>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,F)
        call calc_G<<<blocksG,threadsG>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,G)
        stat = cudaDeviceSynchronize()
        call calc_step3<<<blocks_offset,threads_offset>>>(nx,ny,nz,dxi,dyi,dzi,dt,E,F,G,Q3,Q_d)
        stat = cudaDeviceSynchronize()
        call set_cyclic_bc_d(id,nx,ny,nz,Q_d)
        stat = cudaDeviceSynchronize()
      enddo
      Q = Q_d
      call print_vtk(t2,nx,ny,nz,dx,dy,dz,gamma,Q)
    enddo
  end subroutine RungeKutta
end module calc_time_dev
