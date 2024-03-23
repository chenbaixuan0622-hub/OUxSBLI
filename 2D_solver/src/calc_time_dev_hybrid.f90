module calc_time_dev_hybrid
  implicit none
contains
  subroutine RungeKutta(T0,Q)
    use cudafor
    use mod_globals, only : accuracy, id_visc, nx, ny, nt, np, blocksE, blocksF, threadsE, threadsF
    use calc_physical_quantities
    use calc_steps
    use calc_flux
    use calc_flux_hybrid
    use calc_visc
    use set
    use print
    real(8), intent(inout) :: Q(nx,ny,4)
    real(8), intent(in) :: T0(nx,ny)
    integer t1, t2, itr, stat
    integer(kind=2) :: id_muscl1
    integer(kind=4) :: id_muscl2
    real(8), dimension(nx,ny,4), device :: Q_d, Q2, Q3
    real(8), dimension(nx,ny), device :: rho, u, v, p, T
    real(8), dimension(nx-accuracy+1,ny-accuracy,4), device :: E_hybrid, E_tvd, E_keep, E_roe
    real(8), dimension(nx-accuracy,ny-accuracy+1,4), device :: F_hybrid, F_tvd, F_keep, F_roe

    ! copy on GPU
    Q_d = Q
    T = T0
    do t2 = 1, np
      do t1 = 1, nt
        call calc_quantities(Q_d,rho,u,v,p,T)
        !print *, trim(cudaGetErrorString(cudaGetLastError()))

        call calc_E<<<blocksE,threadsE>>>(id_muscl1,rho,u,v,p,E_keep)
        call calc_F<<<blocksF,threadsF>>>(id_muscl1,rho,u,v,p,F_keep)
        call calc_E<<<blocksE,threadsE>>>(id_muscl2,rho,u,v,p,E_roe)
        call calc_F<<<blocksF,threadsF>>>(id_muscl2,rho,u,v,p,F_roe)
        stat = cudaDeviceSynchronize()
        
        call calc_E_tvd<<<blocksE,threadsE>>>(Q_d,E_keep,E_roe,E_tvd)
        call calc_E_hybrid<<<blocksE,threadsE>>>(u,v,E_keep,E_tvd,E_hybrid)
        call calc_F_tvd<<<blocksF,threadsF>>>(Q_d,F_keep,F_roe,F_tvd)
        call calc_F_hybrid<<<blocksF,threadsF>>>(u,v,F_keep,F_tvd,F_hybrid)

        if (id_visc == 1) then 
          call calc_Ev<<<blocksE,threadsE>>>(u,v,T,E_hybrid)
          call calc_Fv<<<blocksF,threadsF>>>(u,v,T,F_hybrid)
        endif
        stat = cudaDeviceSynchronize()

        call calc_step1(E_hybrid,F_hybrid,Q_d,Q2)
        call set_bc(Q2)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        call calc_quantities(Q2,rho,u,v,p,T)
        
        call calc_E<<<blocksE,threadsE>>>(id_muscl1,rho,u,v,p,E_keep)
        call calc_F<<<blocksF,threadsF>>>(id_muscl1,rho,u,v,p,F_keep)
        call calc_E<<<blocksE,threadsE>>>(id_muscl2,rho,u,v,p,E_roe)
        call calc_F<<<blocksF,threadsF>>>(id_muscl2,rho,u,v,p,F_roe)
        stat = cudaDeviceSynchronize()
        
        call calc_E_tvd<<<blocksE,threadsE>>>(Q2,E_keep,E_roe,E_tvd)
        call calc_E_hybrid<<<blocksE,threadsE>>>(u,v,E_keep,E_tvd,E_hybrid)
        call calc_F_tvd<<<blocksF,threadsF>>>(Q2,F_keep,F_roe,F_tvd)
        call calc_F_hybrid<<<blocksF,threadsF>>>(u,v,F_keep,F_tvd,F_hybrid)
        
        if (id_visc == 1) then
          call calc_Ev<<<blocksE,threadsE>>>(u,v,T,E_hybrid)
          call calc_Fv<<<blocksF,threadsF>>>(u,v,T,F_hybrid)
        endif      
        stat = cudaDeviceSynchronize()

        call calc_step2(E_hybrid,F_hybrid,Q_d,Q2,Q3)
        call set_bc(Q3)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        call calc_quantities(Q3,rho,u,v,p,T)

        call calc_E<<<blocksE,threadsE>>>(id_muscl1,rho,u,v,p,E_keep)
        call calc_F<<<blocksF,threadsF>>>(id_muscl1,rho,u,v,p,F_keep)
        call calc_E<<<blocksE,threadsE>>>(id_muscl2,rho,u,v,p,E_roe)
        call calc_F<<<blocksF,threadsF>>>(id_muscl2,rho,u,v,p,F_roe)
        stat = cudaDeviceSynchronize()

        call calc_E_tvd<<<blocksE,threadsE>>>(Q3,E_keep,E_roe,E_tvd)
        call calc_E_hybrid<<<blocksE,threadsE>>>(u,v,E_keep,E_tvd,E_hybrid)
        call calc_F_tvd<<<blocksF,threadsF>>>(Q3,F_keep,F_roe,F_tvd)
        call calc_F_hybrid<<<blocksF,threadsF>>>(u,v,F_keep,F_tvd,F_hybrid)
        
        if (id_visc == 1) then
          call calc_Ev<<<blocksE,threadsE>>>(u,v,T,E_hybrid)
          call calc_Fv<<<blocksF,threadsF>>>(u,v,T,F_hybrid)
        endif
        stat = cudaDeviceSynchronize()

        call calc_step3(E_hybrid,F_hybrid,Q3,Q_d)
        call set_bc(Q_d)
      enddo
      Q = Q_d
      call print_vtk(t2,Q)
    enddo
  end subroutine RungeKutta
end module calc_time_dev_hybrid

