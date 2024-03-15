module calc_time_dev
  implicit none
contains
  subroutine RungeKutta(nx,ny,nt,np,dx,dy,dt,gamma,T0,Q)
    use cudafor
    use mod_globals, only : accuracy, id_visc, id_scheme, &
            & blocksE, blocksF, blocks, threadsE, threadsF, threads
    use calc_physical_quantities
    use calc_steps
    use calc_KEEP, calc_E_KEEP => calc_E, calc_F_KEEP => calc_F
    use calc_SLAU, calc_E_SLAU => calc_E, calc_F_SLAU => calc_F
    use calc_visc
    use set
    use print
    integer, intent(in) :: nx, ny, nt, np
    real(8), intent(in) :: dx, dy, dt, gamma
    real(8), intent(inout) :: Q(nx,ny,4)
    real(8), intent(in) :: T0(nx,ny)
    integer t1, t2, itr
    integer(kind=2**(accuracy/2)) :: id
    real(8) k, b, dxi, dyi, dtdx, dtdy
    ! GPU
    integer stat, i, j
    type(cudaDeviceProp) :: prop
    real(8), dimension(nx,ny,4), device :: Q_d, Q2, Q3
    real(8), dimension(nx,ny), device :: rho, u, v, p, T
    real(8), device :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), device :: F(nx-accuracy,ny-accuracy+1,4)
    real(8), device :: Ev(nx-accuracy+1,ny-accuracy,4)
    real(8), device :: Fv(nx-accuracy,ny-accuracy+1,4)
    dxi = 1.0d0 / dx
    dyi = 1.0d0 / dy
    dtdx = dt * dxi
    dtdy = dt * dyi

    k = 0.d0
    b = (3.d0 - k) / (1.d0 - k)

    ! copy on GPU
    Q_d = Q
    T = T0
    Ev = 0.d0
    Fv = 0.d0
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

        call calc_step1<<<blocks,threads>>>(nx,ny,dtdx,dtdy,E,F,Ev,Fv,Q_d,Q2)
        call set_bc(nx,ny,gamma,Q2)
        
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

        call calc_step2<<<blocks,threads>>>(nx,ny,dtdx,dtdy,E,F,Ev,Fv,Q_d,Q2,Q3)
        call set_bc(nx,ny,gamma,Q3)

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

        call calc_step3<<<blocks,threads>>>(nx,ny,dtdx,dtdy,E,F,Ev,Fv,Q3,Q_d)
        call set_bc(nx,ny,gamma,Q_d)
      enddo
      Q = Q_d
      call print_vtk(t2,nx,ny,dx,dy,gamma,Q)
    enddo
  end subroutine RungeKutta
end module calc_time_dev

