module calc_time_dev
  implicit none
contains
  subroutine RungeKutta(T0,Q)
    use cudafor
    use mod_globals, only : accuracy, id_muscl, id_visc, nx, ny, nt, np, blocksE, blocksF, threadsE, threadsF
    use calc_physical_quantities
    use calc_steps
    use calc_flux, only : calc_E, calc_F
    use calc_visc
    use set
    use print
    real(8), intent(inout) :: Q(nx,ny,4)
    real(8), intent(in) :: T0(nx,ny)
    integer t1, t2, itr, stat
    real(8), dimension(nx,ny,4), device :: Q_d, Q2, Q3
    real(8), dimension(nx,ny), device :: rho, u, v, p, T
    real(8), device :: E(nx-accuracy+1,ny-accuracy,4)
    real(8), device :: F(nx-accuracy,ny-accuracy+1,4)

    ! copy on GPU
    Q_d = Q
    T = T0
    do t2 = 1, np
      do t1 = 1, nt
        call calc_quantities(Q_d,rho,u,v,p,T)
        !print *, trim(cudaGetErrorString(cudaGetLastError()))

        call calc_E<<<blocksE,threadsE>>>(id_muscl,rho,u,v,p,E)
        call calc_F<<<blocksF,threadsF>>>(id_muscl,rho,u,v,p,F)
        stat = cudaDeviceSynchronize()

        if (id_visc == 1) then 
          call calc_Ev<<<blocksE,threadsE>>>(u,v,T,E)
          call calc_Fv<<<blocksF,threadsF>>>(u,v,T,F)
        endif
        stat = cudaDeviceSynchronize()

        call calc_step1(E,F,Q_d,Q2)
        call set_bc(Q2)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        call calc_quantities(Q2,rho,u,v,p,T)
        
        call calc_E<<<blocksE,threadsE>>>(id_muscl,rho,u,v,p,E)
        call calc_F<<<blocksF,threadsF>>>(id_muscl,rho,u,v,p,F)
        stat = cudaDeviceSynchronize()
        
        if (id_visc == 1) then
          call calc_Ev<<<blocksE,threadsE>>>(u,v,T,E)
          call calc_Fv<<<blocksF,threadsF>>>(u,v,T,F)
        endif      
        stat = cudaDeviceSynchronize()

        call calc_step2(E,F,Q_d,Q2,Q3)
        call set_bc(Q3)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        call calc_quantities(Q3,rho,u,v,p,T)

        call calc_E<<<blocksE,threadsE>>>(id_muscl,rho,u,v,p,E)
        call calc_F<<<blocksF,threadsF>>>(id_muscl,rho,u,v,p,F)
        stat = cudaDeviceSynchronize()

        if (id_visc == 1) then
          call calc_Ev<<<blocksE,threadsE>>>(u,v,T,E)
          call calc_Fv<<<blocksF,threadsF>>>(u,v,T,F)
        endif
        stat = cudaDeviceSynchronize()

        call calc_step3(E,F,Q3,Q_d)
        call set_bc(Q_d)
      enddo
      Q = Q_d
      call print_vtk(t2,Q)
    enddo
  end subroutine RungeKutta
end module calc_time_dev

