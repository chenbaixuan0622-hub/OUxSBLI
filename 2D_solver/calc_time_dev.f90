module calc_time_dev
  implicit none
contains
  subroutine RungeKutta(nx,ny,nt,np,dx,dy,dt,gamma,T0,Q)
    use iso_fortran_env
    use cudafor
    use mod_globals, only : accuracy, id_visc, id_scheme
    use calc_physical_quantities
    use calc_steps
    use calc_KEEP, calc_E_KEEP => calc_E, calc_F_KEEP => calc_F
    use calc_SLAU, calc_E_SLAU => calc_E, calc_F_SLAU => calc_F
    use calc_visc
    use set_bc
    use print
    integer, intent(in) :: nx, ny, nt, np
    real(8), intent(in) :: dx, dy, dt, gamma
    real(8), intent(inout) :: Q(nx,ny,4)
    real(8), intent(in) :: T0(nx,ny)
    integer t1, t2, itr
    integer(kind=2**(accuracy/2)) :: id
    real(8) k, b, dxi, dyi, dtdx, dtdy
    ! GPU
    integer stat, len
    type(cudaDeviceProp) :: prop
    type(dim3) :: blocksE, blocksF, threadsE, threadsF, blocks, threads
    real(8), dimension(nx,ny,4), device :: Q_d, Q2, Q3
    real(8), dimension(nx,ny), device :: rho, u, v, p, T
    real(8), device :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), device :: F(nx-accuracy,ny-accuracy+1,4)
    real(8), device :: Ev(nx-accuracy+1,ny-accuracy,4)
    real(8), device :: Fv(nx-accuracy,ny-accuracy+1,4)
    real(8), device :: ke(nt*np)
    real(8), device :: s(nt*np)
    dxi = 1.0d0 / dx
    dyi = 1.0d0 / dy
    dtdx = dt * dxi
    dtdy = dt * dyi

    k = 0.d0
    b = (3.d0 - k) / (1.d0 - k)

    ! check active device
    stat = cudaSetDevice(0)
    stat = cudaGetDeviceProperties(prop,0)
    len = verify(prop%name, ' ', .true.) 

    ! thread num must be less than 1024
    if (accuracy == 2) then
      blocksE = dim3((nx-accuracy+1)/32,(ny-accuracy)/5,1)
      blocksF = dim3((nx-accuracy)/5,(ny-accuracy+1)/32,1)
      threadsE = dim3(32,5,1)
      threadsF = dim3(5,32,1)
      blocks = dim3((nx-accuracy)/5,(ny-accuracy)/5,1)
      threads = dim3(5,5,1)
    else if (accuracy == 4) then
      blocksE = dim3((nx-accuracy+1),(ny-accuracy)/16,1)
      blocksF = dim3((nx-accuracy)/16,(ny-accuracy+1),1)
      threadsE = dim3(1,16,1)
      threadsF = dim3(16,1,1)
    endif

    ! copy on GPU
    Q_d = Q
    T = T0
    do t2 = 1, np
      do t1 = 1, nt
        if (id_visc == 1) then
          call calc_quantities_T(nx,ny,gamma,Q_d,rho,u,v,p,T)
        else
          call calc_quantities(nx,ny,gamma,Q_d,rho,u,v,p)
        endif

        if (id_scheme == 1) then
          call calc_E_KEEP<<<blocksE,threadsE>>>(id,nx,ny,gamma,rho,u,v,p,E)
          call calc_F_KEEP<<<blocksF,threadsF>>>(id,nx,ny,gamma,rho,u,v,p,F)
          !print *, trim(cudaGetErrorString(cudaGetLastError()))
        else
          call calc_E_SLAU<<<blocksE,threadsE>>>(nx,ny,gamma,k,b,rho,u,v,p,E)
          call calc_F_SLAU<<<blocksF,threadsF>>>(nx,ny,gamma,k,b,rho,u,v,p,F)
        endif

        if (id_visc == 1) then 
          call calc_Ev<<<blocksE,threadsE>>>(nx,ny,dxi,dyi,u,v,T,Ev)
          call calc_Fv<<<blocksF,threadsF>>>(nx,ny,dxi,dyi,u,v,T,Fv)
        endif

        stat = cudaDeviceSynchronize()

        if (id_visc == 1) then
          call vecadd(nx-accuracy+1,ny-accuracy,Ev,E)
          call vecadd(nx-accuracy,ny-accuracy+1,Fv,F)
        endif

        call calc_step1<<<blocks,threads>>>(nx,ny,dtdx,dtdy,E,F,Q_d,Q2)
        !call set_tube_bc(id,nx,ny,gamma,Q2)
        call wind_tunnel_with_a_step(nx,ny,gamma,Q2)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if (id_visc == 1) then
          call calc_quantities_T(nx,ny,gamma,Q2,rho,u,v,p,T)
        else
          call calc_quantities(nx,ny,gamma,Q2,rho,u,v,p)
        endif
        
        if (id_scheme == 1) then
          call calc_E_KEEP<<<blocksE,threadsE>>>(id,nx,ny,gamma,rho,u,v,p,E)
          call calc_F_KEEP<<<blocksF,threadsF>>>(id,nx,ny,gamma,rho,u,v,p,F)
        else
          call calc_E_SLAU<<<blocksE,threadsE>>>(nx,ny,gamma,k,b,rho,u,v,p,E)
          call calc_F_SLAU<<<blocksF,threadsF>>>(nx,ny,gamma,k,b,rho,u,v,p,F)
        endif
        
        if (id_visc == 1) then
          call calc_Ev<<<blocksE,threadsE>>>(nx,ny,dxi,dyi,u,v,T,Ev)
          call calc_Fv<<<blocksF,threadsF>>>(nx,ny,dxi,dyi,u,v,T,Fv)
        endif      
        
        stat = cudaDeviceSynchronize()

        if (id_visc == 1) then
          call vecadd(nx-accuracy+1,ny-accuracy,Ev,E)
          call vecadd(nx-accuracy,ny-accuracy+1,Fv,F)
        endif

        call calc_step2<<<blocks,threads>>>(nx,ny,dtdx,dtdy,E,F,Q_d,Q2,Q3)
        !call set_tube_bc(id,nx,ny,gamma,Q3)
        call wind_tunnel_with_a_step(nx,ny,gamma,Q3)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if (id_visc == 1) then
          call calc_quantities_T(nx,ny,gamma,Q3,rho,u,v,p,T)
        else
          call calc_quantities(nx,ny,gamma,Q3,rho,u,v,p)
        endif

        if (id_scheme == 1) then
          call calc_E_KEEP<<<blocksE,threadsE>>>(id,nx,ny,gamma,rho,u,v,p,E)
          call calc_F_KEEP<<<blocksF,threadsF>>>(id,nx,ny,gamma,rho,u,v,p,F)
        else
          call calc_E_SLAU<<<blocksE,threadsE>>>(nx,ny,gamma,k,b,rho,u,v,p,E)
          call calc_F_SLAU<<<blocksF,threadsF>>>(nx,ny,gamma,k,b,rho,u,v,p,F)
        endif
        
        if (id_visc == 1) then
          call calc_Ev<<<blocksE,threadsE>>>(nx,ny,dxi,dyi,u,v,T,Ev)
          call calc_Fv<<<blocksF,threadsF>>>(nx,ny,dxi,dyi,u,v,T,Fv)
        endif
        
        stat = cudaDeviceSynchronize()

        if (id_visc == 1) then
          call vecadd(nx-accuracy+1,ny-accuracy,Ev,E)
          call vecadd(nx-accuracy,ny-accuracy+1,Fv,F)
        endif

        call calc_step3<<<blocks,threads>>>(nx,ny,dtdx,dtdy,E,F,Q3,Q_d)
        !call set_tube_bc(id,nx,ny,gamma,Q_d)
        call wind_tunnel_with_a_step(nx,ny,gamma,Q_d)
      enddo
      Q = Q_d
      call print_vtk(t2,nx,ny,dx,dy,gamma,Q)
    enddo
  end subroutine RungeKutta
end module calc_time_dev

