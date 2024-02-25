module calc_time_dev
  implicit none
contains
  subroutine calc_quantities(nx,ny,nz,gamma,Q,rho,u,v,w,p)
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in) :: gamma
    real(8), intent(in), dimension(nx,ny,nz,5) :: Q
    real(8), intent(out), dimension(nx,ny,nz) :: rho, u, v, w, p
    integer i, j, k
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rho(i,j,k) = Q(i,j,k,1)
          u(i,j,k) = Q(i,j,k,2) / rho(i,j,k)
          v(i,j,k) = Q(i,j,k,3) / rho(i,j,k)
          w(i,j,k) = Q(i,j,k,4) / rho(i,j,k)
          p(i,j,k) = (gamma-1.d0)*(Q(i,j,k,5)-0.5d0*rho(i,j,k)*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2))
        enddo
      enddo
    enddo
  end subroutine

  subroutine RungeKutta(nx,ny,nz,nt,np,dx,dy,dz,dt,gamma,mu,kappa,Cp,Q)
    use iso_fortran_env
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
    real(8), dimension(nx,ny,nz,5) :: Q2, Q3
    real(8), dimension(nx,ny,nz) :: rho, u, v, w, p
    real(8) E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8) F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8) G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    dxi = 1.0d0 / dx
    dyi = 1.0d0 / dy
    dzi = 1.0d0 / dz

    do t2 = 1, np
      do t1 = 1, nt
        call calc_quantities(nx,ny,nz,gamma,Q,rho,u,v,w,p)
        call calc_E(id,nx,ny,nz,gamma,rho,u,v,w,p,E)
        call calc_F(id,nx,ny,nz,gamma,rho,u,v,w,p,F)
        call calc_G(id,nx,ny,nz,gamma,rho,u,v,w,p,G)
        call calc_step1(accuracy,nx,ny,nz,dxi,dyi,dzi,dt,E,F,G,Q,Q2)
        call set_cyclic_bc(id,nx,ny,nz,Q2)

        call calc_quantities(nx,ny,nz,gamma,Q2,rho,u,v,w,p)
        call calc_E(id,nx,ny,nz,gamma,rho,u,v,w,p,E)
        call calc_F(id,nx,ny,nz,gamma,rho,u,v,w,p,F)
        call calc_G(id,nx,ny,nz,gamma,rho,u,v,w,p,G)
        call calc_step2(accuracy,nx,ny,nz,dxi,dyi,dzi,dt,E,F,G,Q,Q2,Q3)
        call set_cyclic_bc(id,nx,ny,nz,Q3)

        call calc_quantities(nx,ny,nz,gamma,Q3,rho,u,v,w,p)
        call calc_E(id,nx,ny,nz,gamma,rho,u,v,w,p,E)
        call calc_F(id,nx,ny,nz,gamma,rho,u,v,w,p,F)
        call calc_G(id,nx,ny,nz,gamma,rho,u,v,w,p,G)
        call calc_step3(accuracy,nx,ny,nz,dxi,dyi,dzi,dt,E,F,G,Q3,Q)
        call set_cyclic_bc(id,nx,ny,nz,Q)
      enddo
      call print_vtk(t2,nx,ny,nz,dx,dy,dz,gamma,Q)
    enddo
  end subroutine RungeKutta
end module calc_time_dev
