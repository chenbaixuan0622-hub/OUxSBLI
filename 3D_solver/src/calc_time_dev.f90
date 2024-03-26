module calc_time_dev
  use cudafor
  use mod_globals, only : accuracy, id_hybrid, id_muscl, id_visc, nx, ny, nz, nt, np, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG
  use calc_physical_quantities
  use calc_steps
  use calc_flux
  use calc_hybrid
  use calc_visc
  use set
  use print
  implicit none
  interface calc_EFG
    module procedure calc_EFG_basic, calc_EFG_hybrid
  end interface
contains
  subroutine calc_EFG_basic(id_hybrid,Q,T,E,F,G)
    integer(kind=2), intent(in) :: id_hybrid
    real(8), intent(in), device :: Q(nx,ny,nz,5)
    real(8), intent(inout), device :: T(nx,ny,nz)
    real(8), intent(out), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p
    integer stat
    call calc_quantities(Q,rho,u,v,w,p,T)

    call calc_E<<<blocksE,threadsE>>>(id_muscl,rho,u,v,w,p,E)
    call calc_F<<<blocksF,threadsF>>>(id_muscl,rho,u,v,w,p,F)
    call calc_G<<<blocksG,threadsG>>>(id_muscl,rho,u,v,w,p,G)
    !print *, trim(cudaGetErrorString(cudaGetLastError()))
    stat = cudaDeviceSynchronize()

    if (id_visc == 1) then
      call calc_Ev<<<blocksE,threadsE>>>(u,v,w,T,E)
      call calc_Fv<<<blocksF,threadsF>>>(u,v,w,T,F)
      call calc_Gv<<<blocksG,threadsG>>>(u,v,w,T,G)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_basic
  
  subroutine calc_EFG_hybrid(id_hybrid,Q,T,E_hybrid,F_hybrid,G_hybrid)
    integer(kind=4), intent(in) :: id_hybrid
    real(8), intent(in), device :: Q(nx,ny,nz,5)
    real(8), intent(inout), device :: T(nx,ny,nz)
    real(8), intent(out), device :: E_hybrid(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), intent(out), device :: F_hybrid(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), intent(out), device :: G_hybrid(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p
    real(8), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E_keep, E_roe, E_tvd
    real(8), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F_keep, F_roe, F_tvd
    real(8), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G_keep, G_roe, G_tvd
    integer stat
    integer(kind=2) :: id_muscl1
    integer(kind=4) :: id_muscl2
    call calc_quantities(Q,rho,u,v,w,p,T)
    call calc_E<<<blocksE,threadsE>>>(id_muscl1,rho,u,v,w,p,E_keep)
    call calc_F<<<blocksF,threadsF>>>(id_muscl1,rho,u,v,w,p,F_keep)
    call calc_G<<<blocksG,threadsG>>>(id_muscl1,rho,u,v,w,p,G_keep)
    call calc_E<<<blocksE,threadsE>>>(id_muscl2,rho,u,v,w,p,E_roe)
    call calc_F<<<blocksF,threadsF>>>(id_muscl2,rho,u,v,w,p,F_roe)
    call calc_G<<<blocksG,threadsG>>>(id_muscl2,rho,u,v,w,p,G_roe)
    stat = cudaDeviceSynchronize()

    call calc_E_tvd<<<blocksE,threadsE>>>(Q,E_keep,E_roe,E_tvd)
    call calc_F_tvd<<<blocksF,threadsF>>>(Q,F_keep,F_roe,F_tvd)
    call calc_G_tvd<<<blocksG,threadsG>>>(Q,G_keep,G_roe,G_tvd)
    call calc_E_hybrid<<<blocksE,threadsE>>>(u,v,w,E_keep,E_tvd,E_hybrid)
    call calc_F_hybrid<<<blocksF,threadsF>>>(u,v,w,F_keep,F_tvd,F_hybrid)
    call calc_G_hybrid<<<blocksG,threadsG>>>(u,v,w,G_keep,G_tvd,G_hybrid)
    if (id_visc == 1) then
      call calc_Ev<<<blocksE,threadsE>>>(u,v,w,T,E_hybrid)
      call calc_Fv<<<blocksF,threadsF>>>(u,v,w,T,F_hybrid)
      call calc_Gv<<<blocksG,threadsG>>>(u,v,w,T,G_hybrid)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_hybrid

  subroutine RungeKutta(T0,Q)
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    real(8), intent(in) :: T0(nx,ny,nz)
    integer t1, t2, itr, stat
    integer(kind=2) :: id_muscl
    real(8), dimension(nx,ny,nz,5), device :: Q_d, Q2, Q3
    real(8), dimension(nx,ny,nz), device :: T
    real(8), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)

    ! copy on GPU
    Q_d = Q
    T = T0
    do t2 = 1, np
      do t1 = 1, nt
        call calc_EFG(id_hybrid,Q_d,T,E,F,G)
        call calc_step1(E,F,G,Q_d,Q2)
        call set_bc(id_accuracy,Q2)

        call calc_EFG(id_hybrid,Q2,T,E,F,G)
        call calc_step2(E,F,G,Q_d,Q2,Q3)
        call set_bc(id_accuracy,Q3)

        call calc_EFG(id_hybrid,Q3,T,E,F,G)
        call calc_step3(E,F,G,Q3,Q_d)
        call set_bc(id_accuracy,Q_d)
      enddo
      Q = Q_d
      call print_vtk(t2,Q)
    enddo
  end subroutine RungeKutta
end module calc_time_dev

