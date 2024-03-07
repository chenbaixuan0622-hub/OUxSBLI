module calc_time_dev
  implicit none
contains
  subroutine RungeKutta(nx,ny,nz,nt,np,dx,dy,dz,dt,gamma,T0,Cs,Q)
    use iso_fortran_env
    use cudafor
    use mod_globals, only : accuracy, id_visc, id_turbulence
    use calc_physical_quantities
    use calc_steps
    use calc_flux
    use calc_visc
    use set_bc
    use print
    integer, intent(in) :: nx, ny, nz, nt, np
    real(8), intent(in) :: dx, dy, dz, dt, gamma, Cs
    real(8), intent(in) :: T0(nx,ny,nz)
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer t1, t2, itr
    integer(kind=2**(accuracy/2)) :: id
    real(8) dxi, dyi, dzi, delta
    ! GPU
    integer stat, len
    type(cudaDeviceProp) :: prop
    type(dim3) :: blocksE, blocksF, blocksG, threadsE, threadsF, threadsG
    real(8), dimension(nx,ny,nz,5), device :: Q_d, Q2, Q3
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p, T
    real(8), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), device :: Ev(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), device :: Fv(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), device :: Gv(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), device :: ke(nt*np)
    real(8), device :: s(nt*np)
    dxi = 1.0d0 / dx
    dyi = 1.0d0 / dy
    dzi = 1.0d0 / dz
    delta = (dx * dy * dz) ** (1.d0 / 3.d0)

    ! check active device
    stat = cudaSetDevice(0)
    stat = cudaGetDeviceProperties(prop,0)
    len = verify(prop%name, ' ', .true.) 

    ! thread num must be less than 1024
    if (accuracy == 2) then
      blocksE = dim3((nx-accuracy+1)/32,(ny-accuracy)/5,(nz-accuracy)/5)
      blocksF = dim3((nx-accuracy)/5,(ny-accuracy+1)/32,(nz-accuracy)/5)
      blocksG = dim3((nx-accuracy)/5,(ny-accuracy)/5,(nz-accuracy+1)/32)
      threadsE = dim3(32,5,5)
      threadsF = dim3(5,32,5)
      threadsG = dim3(5,5,32)
    else if (accuracy == 4) then
      blocksE = dim3((nx-accuracy+1)/2,(ny-accuracy)/11,(nz-accuracy)/11)
      blocksF = dim3((nx-accuracy)/11,(ny-accuracy+1)/2,(nz-accuracy)/11)
      blocksG = dim3((nx-accuracy)/11,(ny-accuracy)/11,(nz-accuracy+1)/2)
      threadsE = dim3(2,11,11)
      threadsF = dim3(11,2,11)
      threadsG = dim3(11,11,2)
    endif

    ! copy on GPU
    Q_d = Q
    T = T0
    do t2 = 1, np
      do t1 = 1, nt
        if (id_visc == 1 .or. id_turbulence == 1) then
          call calc_quantities_T(nx,ny,nz,gamma,Q_d,rho,u,v,w,p,T)
        else
          call calc_quantities(nx,ny,nz,gamma,Q_d,rho,u,v,w,p)
        endif

        call calc_E<<<blocksE,threadsE>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,E)
        call calc_F<<<blocksF,threadsF>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,F)
        call calc_G<<<blocksG,threadsG>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,G)
        !print *, trim(cudaGetErrorString(cudaGetLastError()))

        if (id_visc == 1 .or. id_turbulence == 1) then
          call calc_Ev<<<blocksE,threadsE>>>(nx,ny,nz,dxi,dyi,dzi,delta,Cs,rho,u,v,w,T,Ev)
          call calc_Fv<<<blocksF,threadsF>>>(nx,ny,nz,dxi,dyi,dzi,delta,Cs,rho,u,v,w,T,Fv)
          call calc_Gv<<<blocksG,threadsG>>>(nx,ny,nz,dxi,dyi,dzi,delta,Cs,rho,u,v,w,T,Gv)
        endif

        stat = cudaDeviceSynchronize()

        if (id_visc == 1 .or. id_turbulence == 1) then
          call calc_step1(nx,ny,nz,dxi,dyi,dzi,dt,E,F,G,Ev,Fv,Gv,Q_d,Q2)
        else
          call calc_step1(nx,ny,nz,dxi,dyi,dzi,dt,E,F,G,Q_d,Q2)
        endif
        
        call set_cyclic_bc_d(id,nx,ny,nz,Q2)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (id_visc == 1 .or. id_turbulence == 1) then
          call calc_quantities_T(nx,ny,nz,gamma,Q2,rho,u,v,w,p,T)
        else
          call calc_quantities(nx,ny,nz,gamma,Q2,rho,u,v,w,p)
        endif

        call calc_E<<<blocksE,threadsE>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,E)
        call calc_F<<<blocksF,threadsF>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,F)
        call calc_G<<<blocksG,threadsG>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,G)
        
        if (id_visc == 1 .and. id_turbulence == 1) then
          call calc_Ev<<<blocksE,threadsE>>>(nx,ny,nz,dxi,dyi,dzi,delta,Cs,rho,u,v,w,T,Ev)
          call calc_Fv<<<blocksF,threadsF>>>(nx,ny,nz,dxi,dyi,dzi,delta,Cs,rho,u,v,w,T,Fv)
          call calc_Gv<<<blocksG,threadsG>>>(nx,ny,nz,dxi,dyi,dzi,delta,Cs,rho,u,v,w,T,Gv)
        endif      
        
        stat = cudaDeviceSynchronize()

        if (id_visc == 1 .or. id_turbulence == 1) then
          call calc_step2(nx,ny,nz,dxi,dyi,dzi,dt,E,F,G,Ev,Fv,Gv,Q_d,Q2,Q3)
        else
          call calc_step2(nx,ny,nz,dxi,dyi,dzi,dt,E,F,G,Q_d,Q2,Q3)
        endif

        call set_cyclic_bc_d(id,nx,ny,nz,Q3)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if (id_visc == 1 .or. id_turbulence == 1) then
          call calc_quantities_T(nx,ny,nz,gamma,Q3,rho,u,v,w,p,T)
        else
          call calc_quantities(nx,ny,nz,gamma,Q3,rho,u,v,w,p)
        endif
        
        call calc_E<<<blocksE,threadsE>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,E)
        call calc_F<<<blocksF,threadsF>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,F)
        call calc_G<<<blocksG,threadsG>>>(id,nx,ny,nz,gamma,rho,u,v,w,p,G)
        
        if (id_visc == 1 .or. id_turbulence == 1) then
          call calc_Ev<<<blocksE,threadsE>>>(nx,ny,nz,dxi,dyi,dzi,delta,Cs,rho,u,v,w,T,Ev)
          call calc_Fv<<<blocksF,threadsF>>>(nx,ny,nz,dxi,dyi,dzi,delta,Cs,rho,u,v,w,T,Fv)
          call calc_Gv<<<blocksG,threadsG>>>(nx,ny,nz,dxi,dyi,dzi,delta,Cs,rho,u,v,w,T,Gv)
        endif
        
        stat = cudaDeviceSynchronize()

        if (id_visc == 1 .or. id_turbulence == 1) then
          call calc_step3(nx,ny,nz,dxi,dyi,dzi,dt,E,F,G,Ev,Fv,Gv,Q3,Q_d)
        else
          call calc_step3(nx,ny,nz,dxi,dyi,dzi,dt,E,F,G,Q3,Q_d)
        endif
        
        call set_cyclic_bc_d(id,nx,ny,nz,Q_d)
      enddo
      Q = Q_d
      call print_vtk(t2,nx,ny,nz,dx,dy,dz,gamma,Q)
    enddo
  end subroutine RungeKutta
end module calc_time_dev
