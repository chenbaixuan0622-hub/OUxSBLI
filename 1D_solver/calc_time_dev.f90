module calc_time_dev
  implicit none
contains
  subroutine RungeKutta(nx,nt,np,dx,dt,gamma,T0,Q)
    use iso_fortran_env
    use mod_globals, only : accuracy, id_visc, id_scheme
    use calc_physical_quantities
    use calc_steps
    !use calc_KEEP
    use calc_SLAU
    use calc_visc
    use set_bc
    use print
    integer, intent(in) :: nx, nt, np
    real(8), intent(in) :: dx, dt, gamma
    real(8), intent(inout) :: Q(nx,3)
    real(8), intent(in) :: T0(nx)
    integer t1, t2, itr
    integer(kind=2**(accuracy/2)) :: id
    real(8) dxi, dtdx
    real(8) k, b
    real(8), dimension(nx,3) :: Q2, Q3, Qs
    real(8), dimension(nx) :: rho, u, p, T
    real(8) :: E(nx-accuracy+1,3)
    real(8) :: Ev(nx-accuracy+1,3)
    dxi = 1.0d0 / dx
    dtdx = dt * dxi

    k = 0.d0
    b = (3.d0 - k) / (1.d0 - k)

    T = T0
    do t2 = 1, np
      do t1 = 1, nt
        if (id_visc == 1) then
          call calc_quantities_T(nx,gamma,Q,rho,u,p,T)
        else
          call calc_quantities(nx,gamma,Q,rho,u,p)
        endif

        Qs(:,1) = rho(:)
        Qs(:,2) = u(:)
        Qs(:,3) = p(:)
        call calc_E(nx,gamma,k,b,Qs,E)
        !call calc_E(id,nx,gamma,rho,u,p,E)

        if (id_visc == 1) then 
          call calc_Ev(nx,dxi,u,T,Ev)
        endif

        if (id_visc == 1) then
          call calc_step1(nx,dtdx,E,Ev,Q,Q2)
        else
          call calc_step1(nx,dtdx,E,Q,Q2)
        endif
        
        call set_tube_bc(id,nx,gamma,Q2)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if (id_visc == 1) then
          call calc_quantities_T(nx,gamma,Q2,rho,u,p,T)
        else
          call calc_quantities(nx,gamma,Q2,rho,u,p)
        endif
        
        Qs(:,1) = rho(:)
        Qs(:,2) = u(:)
        Qs(:,3) = p(:)
        call calc_E(nx,gamma,k,b,Qs,E)
        !call calc_E(id,nx,gamma,rho,u,p,E)
        
        if (id_visc == 1) then
          call calc_Ev(nx,dxi,u,T,Ev)
        endif      

        if (id_visc == 1) then
          call calc_step2(nx,dtdx,E,Ev,Q,Q2,Q3)
        else
          call calc_step2(nx,dtdx,E,Q,Q2,Q3)
        endif

        call set_tube_bc(id,nx,gamma,Q3)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if (id_visc == 1) then
          call calc_quantities_T(nx,gamma,Q3,rho,u,p,T)
        else
          call calc_quantities(nx,gamma,Q3,rho,u,p)
        endif

        Qs(:,1) = rho(:)
        Qs(:,2) = u(:)
        Qs(:,3) = p(:)
        call calc_E(nx,gamma,k,b,Qs,E)
        !call calc_E(id,nx,gamma,rho,u,p,E)
        
        if (id_visc == 1) then
          call calc_Ev(nx,dxi,u,T,Ev)
        endif

        if (id_visc == 1) then
          call calc_step3(nx,dtdx,E,Ev,Q3,Q)
        else
          call calc_step3(nx,dtdx,E,Q3,Q)
        endif
        
        call set_tube_bc(id,nx,gamma,Q)
      enddo
      call print_vtk(t2,nx,dx,gamma,Q)
    enddo
  end subroutine RungeKutta
end module calc_time_dev

