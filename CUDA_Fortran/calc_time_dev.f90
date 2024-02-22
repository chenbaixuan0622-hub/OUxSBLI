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
    real(8), dimension(nx,nyp,nzp,5,4) :: Q2, Q3
    integer tc, tp, threadNum, stat
    ! GPU
    real(8), dimension(nx,nyp,nzp,5,4), device :: Q_d, Q2_d, Q3_d
    real(8), dimension(ni+4,nj+4,nk+4,4), device :: Q_local
    real(8), dimension(ni,nj,nk,4), device :: R
    ! define block size
    threadNum = blockDim(1) * blockDim(2) * blockDim(3)
    
    ! send variable to GPU
    Q_d = Q
    do tp = 1, np
      do tc = 1, nt
        ! GPU
        call calc_step1<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,dX,dY,dZ,dt,gamma,mu,kappa,Cp,Q_d,Q2_d,Q_local,R)
        stat = cudaDeviceSynchronize()
        call set_cyclic_bc_d(nx,nyp,nzp,Q2_d)
        
        call calc_step2<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,dX,dY,dZ,dt,gamma,mu,kappa,Cp,Q_d,Q2_d,Q3_d,Q_local,R)
        stat = cudaDeviceSynchronize()
        call set_cyclic_bc_d(nx,nyp,nzp,Q3_d)
        
        call calc_step3<<<4, threadNum>>>(nx,nyp,nzp,ni,nj,nk,dX,dY,dZ,dt,gamma,mu,kappa,Cp,Q3_d,Q_d,Q_local,R)
        stat = cudaDeviceSynchronize()
        call set_cyclic_bc_d(nx,nyp,nzp,Q_d)
      enddo
      ! CPU
      Q = Q_d
      call print_vtk(tp,nx,ny,nz,nyp,nzp,dX,dY,dZ,gamma,Q)
    enddo
  end subroutine RungeKutta
end module calc_time_dev
