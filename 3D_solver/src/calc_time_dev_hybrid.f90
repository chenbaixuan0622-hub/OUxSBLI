module calc_time_dev_hybrid
  implicit none
contains
  subroutine RungeKutta(T0,Q)
    use cudafor
    use mod_globals, only : accuracy, id_visc, nx, ny, nz, nt, np, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG
    use calc_physical_quantities
    use calc_steps
    use calc_KEEP
    !use calc_riemann_solver
    use calc_visc
    use set
    use print
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    real(8), intent(in) :: T0(nx,ny,nz)
    integer t1, t2, itr, stat
    integer(kind=2) :: id_muscl
    real(8), dimension(nx,ny,nz,5), device :: Q_d, Q2, Q3
    real(8), dimension(nx,ny,nz), device :: rho, u, v, w, p, T
    real(8), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    real(8), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    real(8), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)

    ! copy on GPU
    Q_d = Q
    T = T0
    do t2 = 1, np
      do t1 = 1, nt
        call calc_quantities(Q_d,rho,u,v,w,p,T)

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

        call calc_step1(E,F,G,Q_d,Q2)
        call set_bc(id_accuracy,Q2)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call calc_quantities(Q2,rho,u,v,w,p,T)

        call calc_E<<<blocksE,threadsE>>>(id_muscl,rho,u,v,w,p,E)
        call calc_F<<<blocksF,threadsF>>>(id_muscl,rho,u,v,w,p,F)
        call calc_G<<<blocksG,threadsG>>>(id_muscl,rho,u,v,w,p,G)
        stat = cudaDeviceSynchronize()
        
        if (id_visc == 1) then
          call calc_Ev<<<blocksE,threadsE>>>(u,v,w,T,E)
          call calc_Fv<<<blocksF,threadsF>>>(u,v,w,T,F)
          call calc_Gv<<<blocksG,threadsG>>>(u,v,w,T,G)
        endif      
        stat = cudaDeviceSynchronize()

        call calc_step2(E,F,G,Q_d,Q2,Q3)
        call set_bc(id_accuracy,Q3)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        call calc_quantities(Q3,rho,u,v,w,p,T)
        
        call calc_E<<<blocksE,threadsE>>>(id_muscl,rho,u,v,w,p,E)
        call calc_F<<<blocksF,threadsF>>>(id_muscl,rho,u,v,w,p,F)
        call calc_G<<<blocksG,threadsG>>>(id_muscl,rho,u,v,w,p,G)
        stat = cudaDeviceSynchronize()
        
        if (id_visc == 1) then
          call calc_Ev<<<blocksE,threadsE>>>(u,v,w,T,E)
          call calc_Fv<<<blocksF,threadsF>>>(u,v,w,T,F)
          call calc_Gv<<<blocksG,threadsG>>>(u,v,w,T,G)
        endif
        stat = cudaDeviceSynchronize()

        call calc_step3(E,F,G,Q3,Q_d)
        call set_bc(id_accuracy,Q_d)
      enddo
      Q = Q_d
      call print_vtk(t2,Q)
    enddo
  end subroutine RungeKutta
end module calc_time_dev_hybrid

