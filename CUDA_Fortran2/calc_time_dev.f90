module calc_time_dev
  implicit none
contains
  subroutine RungeKutta(nx,ny,nz,nyp,nzp,ni,nj,nk,blockDim,nt,np,dX,dY,dZ,dt,gamma,mu,kappa,Cp,Q)
    use cudafor
    use calc_steps
    use set_bc
    use print
    ! CPU
    integer, intent(in) :: nx, ny, nz, nyp, nzp, ni, nj, nk, nt, np
    integer, intent(in) :: blockDim(3)
    real(8), intent(in) :: dX, dY, dZ, dt, gamma, mu, kappa, Cp
    real(8), intent(inout) :: Q(nx,nyp,nzp,5,4)
    integer tc, tp, threadNum, stat
    ! GPU
    real(8), allocatable, device :: Q_d(:,:,:,:,:), Q2(:,:,:,:,:), Q3(:,:,:,:,:)
    real(8), allocatable, device :: E(:,:,:,:,:)
    real(8), allocatable, device :: F(:,:,:,:,:)
    real(8), allocatable, device :: G(:,:,:,:,:)
    real(8), allocatable, device :: rho(:,:,:,:), u(:,:,:,:), v(:,:,:,:), w(:,:,:,:), p(:,:,:,:)
    allocate(Q_d(nx,nyp,nzp,5,4),Q2(nx,nyp,nzp,5,4),Q3(nx,nyp,nzp,5,4))
    allocate(E(nx-3,nyp-4,nzp-4,5,4),F(nx-4,nyp-3,nzp-4,5,4),G(nx-4,nyp-4,nzp-3,5,4))
    allocate(rho(nx,nyp,nzp,4),u(nx,nyp,nzp,4),v(nx,nyp,nzp,4),w(nx,nyp,nzp,4),p(nx,nyp,nzp,4))
    ! define block size
    threadNum = blockDim(1) * blockDim(2) * blockDim(3)
    
    ! send variable to GPU
    Q_d = Q
    do tp = 1, np
      do tc = 1, nt
        ! GPU
        call calc_quantity<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,gamma,Q_d,rho,u,v,w,p)
        stat = cudaDeviceSynchronize()
        !call calc_E<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,gamma,u,v,w,p,E)
        !call calc_F<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,gamma,u,v,w,p,F)
        !call calc_G<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,gamma,u,v,w,p,G)
        stat = cudaDeviceSynchronize()
        call calc_step1<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,dX,dY,dZ,dt,Q_d,Q2,E,F,G)
        stat = cudaDeviceSynchronize()
        call set_cyclic_bc_d(nx,nyp,nzp,Q2)

        call calc_quantity<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,gamma,Q2,rho,u,v,w,p)
        stat = cudaDeviceSynchronize()
        !call calc_E<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,gamma,u,v,w,p,E)
        !call calc_F<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,gamma,u,v,w,p,F)
        !call calc_G<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,gamma,u,v,w,p,G)
        stat = cudaDeviceSynchronize()
        call calc_step2<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,dX,dY,dZ,dt,Q_d,Q2,Q3,E,F,G)
        stat = cudaDeviceSynchronize()
        call set_cyclic_bc_d(nx,nyp,nzp,Q3)
        
        call calc_quantity<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,gamma,Q3,rho,u,v,w,p)
        stat = cudaDeviceSynchronize()
        !call calc_E<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,gamma,u,v,w,p,E)
        !call calc_F<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,gamma,u,v,w,p,F)
        !call calc_G<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,gamma,u,v,w,p,G)
        stat = cudaDeviceSynchronize()
        call calc_step3<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,dX,dY,dZ,dt,Q3,Q_d,E,F,G)
        stat = cudaDeviceSynchronize()
        call set_cyclic_bc_d(nx,nyp,nzp,Q_d)
      enddo
      ! CPU
      Q = Q_d
      call print_vtk(tp,nx,ny,nz,nyp,nzp,dX,dY,dZ,gamma,Q)
    enddo
    deallocate(Q_d,Q2,Q3,E,F,G,rho,u,v,w,p)
  end subroutine RungeKutta
end module calc_time_dev
