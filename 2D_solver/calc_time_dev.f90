module calc_time_dev
  implicit none
contains
  subroutine RungeKutta(nx,ny,nt,np,dx,dy,dt,gamma,mu,kappa,Cv,Cp,Q)
    use iso_fortran_env
    use cudafor
    use mod_globals, only : accuracy, id_visc
    use calc_physical_quantities
    use calc_steps
    use calc_flux
    use calc_visc
    use set_bc
    use print
    integer, intent(in) :: nx, ny, nt, np
    real(8), intent(in) :: dx, dy, dt, gamma, kappa, Cv, Cp
    real(8), intent(inout) :: mu
    real(8), intent(inout) :: Q(nx,ny,4)
    integer t1, t2, itr
    integer(kind=2**(accuracy/2)) :: id
    real(8) dxi, dyi
    ! GPU
    integer stat, len
    type(cudaDeviceProp) :: prop
    type(dim3) :: blocksE, blocksF, threadsE, threadsF
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

    ! check active device
    stat = cudaSetDevice(0)
    stat = cudaGetDeviceProperties(prop,0)
    len = verify(prop%name, ' ', .true.) 

    ! thread num must be less than 1024
    if (accuracy == 2) then
      blocksE = dim3((nx-accuracy+1)/32,(ny-accuracy)/5,(nz-accuracy)/5)
      blocksF = dim3((nx-accuracy)/5,(ny-accuracy+1)/32,(nz-accuracy)/5)
      threadsE = dim3(32,5,5)
      threadsF = dim3(5,32,5)
    else if (accuracy == 4) then
      blocksE = dim3((nx-accuracy+1),(ny-accuracy)/16,(nz-accuracy)/16)
      blocksF = dim3((nx-accuracy)/16,(ny-accuracy+1),(nz-accuracy)/16)
      threadsE = dim3(1,16,16)
      threadsF = dim3(16,1,16)
    endif

    ! copy on GPU
    Q_d = Q
    do t2 = 1, np
      do t1 = 1, nt
        call calc_quantities(nx,ny,gamma,Cp,Q_d,rho,u,v,p,T)
        
        ! Euler
        call calc_E<<<blocksE,threadsE>>>(id,nx,ny,gamma,rho,u,v,p,E)
        call calc_F<<<blocksF,threadsF>>>(id,nx,ny,gamma,rho,u,v,p,F)
        !print *, trim(cudaGetErrorString(cudaGetLastError()))

        if (id_visc == 1) then 
          call calc_Ev<<<blocksE,threadsE>>>(id,nx,ny,dxi,dyi,mu,kappa,u,v,T,Ev)
          call calc_Fv<<<blocksF,threadsF>>>(id,nx,ny,dxi,dyi,mu,kappa,u,v,T,Fv)
        endif

        stat = cudaDeviceSynchronize()

        if (id_visc == 1) then
          call calc_step1(nx,ny,dxi,dyi,dt,E,F,Ev,Fv,Q_d,Q2)
        else
          call calc_step1(nx,ny,dxi,dyi,dt,E,F,Q_d,Q2)
        endif
        
        call set_tube_bc(id,nx,ny,Q2)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call calc_quantities(nx,ny,gamma,Cp,Q2,rho,u,v,p,T)
        
        ! Euler
        call calc_E<<<blocksE,threadsE>>>(id,nx,ny,gamma,rho,u,v,p,E)
        call calc_F<<<blocksF,threadsF>>>(id,nx,ny,gamma,rho,u,v,p,F)
        
        if (id_visc == 1) then
          call calc_Ev<<<blocksE,threadsE>>>(id,nx,ny,dxi,dyi,mu,kappa,u,v,T,Ev)
          call calc_Fv<<<blocksF,threadsF>>>(id,nx,ny,dxi,dyi,mu,kappa,u,v,T,Fv)
        endif      
        
        stat = cudaDeviceSynchronize()

        if (id_visc == 1) then
          call calc_step2(nx,ny,dxi,dyi,dt,E,F,Ev,Fv,Q_d,Q2,Q3)
        else
          call calc_step2(nx,ny,dxi,dyi,dt,E,F,Q_d,Q2,Q3)
        endif

        call set_tube_bc(id,nx,ny,Q3)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        call calc_quantities(nx,ny,gamma,Cp,Q3,rho,u,v,p,T)
        
        ! Euler
        call calc_E<<<blocksE,threadsE>>>(id,nx,ny,gamma,rho,u,v,p,E)
        call calc_F<<<blocksF,threadsF>>>(id,nx,ny,gamma,rho,u,v,p,F)
        
        if (id_visc == 1) then
          call calc_Ev<<<blocksE,threadsE>>>(id,nx,ny,dxi,dyi,mu,kappa,u,v,T,Ev)
          call calc_Fv<<<blocksF,threadsF>>>(id,nx,ny,dxi,dyi,mu,kappa,u,v,T,Fv)
        endif
        
        stat = cudaDeviceSynchronize()

        if (id_visc == 1) then
          call calc_step3(nx,ny,dxi,dyi,dt,E,F,Ev,Fv,Q3,Q_d)
        else
          call calc_step3(nx,ny,dxi,dyi,dt,E,F,Q3,Q_d)
        endif
        
        call set_tube_bc(id,nx,ny,Q_d)
      enddo
      Q = Q_d
      call print_vtk(t2,nx,ny,dx,dy,gamma,Q)
    enddo
  end subroutine RungeKutta
end module calc_time_dev
