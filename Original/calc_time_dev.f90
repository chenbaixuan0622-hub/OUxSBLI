module calc_time_dev
  implicit none
contains
  subroutine RungeKutta(nx,ny,nz,nt,np,dx,dy,dz,dt,gamma,mu,kappa,Cp,Cs,Q)
    use iso_fortran_env
    use mod_globals, only : accuracy
    use calc_physical_quantities
    use calc_steps
    use calc_KEEP
    use calc_SLAU
    use calc_visc
    use set_bc
    use print
    integer, intent(in) :: nx, ny, nz, nt, np
    real(8), intent(in) :: dx, dy, dz, dt, gamma, kappa, Cp, Cs
    real(8), intent(inout) :: mu
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    integer t1, t2
    integer(kind=2**(accuracy/2)) :: id
    real(8) dxi, dyi, dzi, delta
    real(8), dimension(nx,ny,nz,5) :: Q2, Q3
    real(8) E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8) F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8) G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8) Ev(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8) Fv(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8) Gv(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    ! KEEP
    real(8), dimension(nx,ny,nz) :: rho, u, v, w, p, T
    dxi = 1.0d0 / dx
    dyi = 1.0d0 / dy
    dzi = 1.0d0 / dz
    delta = (dx * dy * dz) ** (1.d0 / 3.d0)

    Ev = 0.d0
    Fv = 0.d0
    Gv = 0.d0
    mu = 0.d0
    do t2 = 1, np
      do t1 = 1, nt
        call calc_quantities(nx,ny,nz,gamma,Cp,Q,rho,u,v,w,p,T)
        ! Euler
        call calc_E(id,nx,ny,nz,gamma,rho,u,v,w,p,E)
        call calc_F(id,nx,ny,nz,gamma,rho,u,v,w,p,F)
        call calc_G(id,nx,ny,nz,gamma,rho,u,v,w,p,G)
        ! visc 
        !call calc_Ev(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,u,v,w,T,Ev)
        !call calc_Fv(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,u,v,w,T,Fv)
        !call calc_Gv(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,u,v,w,T,Gv)
        ! visc + LES
        call calc_Ev(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,delta,Cs,rho,u,v,w,T,Ev)
        call calc_Fv(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,delta,Cs,rho,u,v,w,T,Fv)
        call calc_Gv(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,delta,Cs,rho,u,v,w,T,Gv)
        call calc_step1(nx,ny,nz,dxi,dyi,dzi,dt,(E-Ev),(F-Fv),(G-Gv),Q,Q2)
        call set_cyclic_bc(id,nx,ny,nz,Q2)

        call calc_quantities(nx,ny,nz,gamma,Cp,Q2,rho,u,v,w,p,T)
        ! Euler
        call calc_E(id,nx,ny,nz,gamma,rho,u,v,w,p,E)
        call calc_F(id,nx,ny,nz,gamma,rho,u,v,w,p,F)
        call calc_G(id,nx,ny,nz,gamma,rho,u,v,w,p,G)
        ! visc
        !call calc_Ev(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,u,v,w,T,Ev)
        !call calc_Fv(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,u,v,w,T,Fv)
        !call calc_Gv(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,u,v,w,T,Gv)
        ! visc + LES
        call calc_Ev(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,delta,Cs,rho,u,v,w,T,Ev)
        call calc_Fv(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,delta,Cs,rho,u,v,w,T,Fv)
        call calc_Gv(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,delta,Cs,rho,u,v,w,T,Gv)
        call calc_step2(nx,ny,nz,dxi,dyi,dzi,dt,(E-Ev),(F-Fv),(G-Gv),Q,Q2,Q3)
        call set_cyclic_bc(id,nx,ny,nz,Q3)

        call calc_quantities(nx,ny,nz,gamma,Cp,Q3,rho,u,v,w,p,T)
        ! Euler
        call calc_E(id,nx,ny,nz,gamma,rho,u,v,w,p,E)
        call calc_F(id,nx,ny,nz,gamma,rho,u,v,w,p,F)
        call calc_G(id,nx,ny,nz,gamma,rho,u,v,w,p,G)
        ! visc
        !call calc_Ev(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,u,v,w,T,Ev)
        !call calc_Fv(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,u,v,w,T,Fv)
        !call calc_Gv(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,u,v,w,T,Gv)
        ! visc + LES
        call calc_Ev(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,delta,Cs,rho,u,v,w,T,Ev)
        call calc_Fv(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,delta,Cs,rho,u,v,w,T,Fv)
        call calc_Gv(id,nx,ny,nz,dxi,dyi,dzi,mu,kappa,delta,Cs,rho,u,v,w,T,Gv)
        call calc_step3(nx,ny,nz,dxi,dyi,dzi,dt,(E-Ev),(F-Fv),(G-Gv),Q3,Q)
        call set_cyclic_bc(id,nx,ny,nz,Q)
      enddo
      call print_vtk(t2,nx,ny,nz,dx,dy,dz,gamma,Q)
    enddo
  end subroutine RungeKutta
end module calc_time_dev
