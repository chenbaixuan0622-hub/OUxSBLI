module calc_rescale
  use cudafor
  use mpi
  use mod_globals, only : id_gpumpi, nre1, nre2, rerank, nt, dt, gamma , R, Pr, u0, rho0, p0, M0, blt, start_rescale, rf, Taw
  use mod_constant, only : Cp, over_Cp, gamma_1, over_gamma_1, mu0_T0_S, over_T0
  use cpu_gpu_mpi
contains
  subroutine calc_mean(step, flag_re, nx, ny, nz, Jacobian, QJ, Qm)
    integer, intent(in)            :: step, flag_re, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny), QJ(5,nx,ny,nz)
    real(8), intent(inout), device :: Qm(ny*2)
    real(8) Q1, Q2, rhoinv, Jacobian_tmp, volinv, step1, step2
    integer i, k, istat
    volinv = 1.d0 / dble((nre2 - nre1 + 1) * (nz - 6))
    if (flag_re == 0) then
      !$cuf kernel do <<<*,*>>>
      do j = 1, ny
        Q1 = 0.d0; Q2 = 0.d0
        do k = 4, nz-3
          do i = nre1, nre2
            Jacobian_tmp = Jacobian(i,j)
            rhoinv = 1.d0 / QJ(1,i,j,k)
            Q1 = Q1 + QJ(2,i,j,k) * rhoinv
            Q2 = Q2 + QJ(3,i,j,k) * rhoinv
        enddo;enddo
        Qm(2*(j-1)+1) = Q1 * volinv
        Qm(2*(j-1)+2) = Q2 * volinv
      enddo
    else
      step1 = dble(step-1); step2 = 1.d0 / dble(step)
      !$cuf kernel do <<<*,*>>>
      do j = 1, ny
        Q1 = 0.d0; Q2 = 0.d0
        do k = 4, nz-3
          do i = nre1, nre2
            Jacobian_tmp = Jacobian(i,j)
            rhoinv = 1.d0 / QJ(1,i,j,k)
            Q1 = Q1 + QJ(2,i,j,k) * rhoinv
            Q2 = Q2 + QJ(3,i,j,k) * rhoinv
        enddo;enddo
        Qm(2*(j-1)+1) = (step1 * Qm(2*(j-1)+1) + Q1 * volinv) * step2
        Qm(2*(j-1)+2) = (step1 * Qm(2*(j-1)+2) + Q2 * volinv) * step2
      enddo
    endif
    istat = cudaDeviceSynchronize() ! This is necessary for reduction (+)
  end subroutine calc_mean


  subroutine copy(nx, ny, nz, QJ, Qre)
    integer, intent(in)          :: nx, ny, nz
    real(8), intent(in), device  :: QJ(5,nx,ny,nz)
    real(8), intent(out), device :: Qre(ny*(nz-6)*5)
    integer j, k, l, j_offset, k_offset
    !$cuf kernel do <<<*,*>>>
    do k = 1, nz-6
      k_offset = ny * 5 * (k-1)
      do j = 1, ny
        j_offset = 5 * (j-1)
        do l = 1, 5
          Qre(k_offset+j_offset+l) = QJ(l,nre2,j,k+3)
    enddo;enddo;enddo
  end subroutine copy


  subroutine step_rescale(num, myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJ, Qm, Qre)
    integer, intent(in)            :: num, myrank, step, nx, ny, nz
    integer, intent(inout)         :: flag_re, ireq, ireq2(2)
    real(8), intent(in), device    :: Jacobian(nx,ny), QJ(5,nx,ny,nz)
    real(8), intent(inout), device :: Qm(ny*2), Qre(ny*(nz-6)*5)
    integer ierr, j
    if (myrank == rerank) then
      call copy(nx, ny, nz, QJ, Qre)
      call CPUGPU_MPI_SEND(id_gpumpi, Qre, 5*ny*(nz-6), rerank+1, 0, MPI_COMM_WORLD, ireq2(1), ierr)
      if (num == 1) then
        call calc_mean(step, flag_re, nx, ny, nz, Jacobian, QJ, Qm)
        call CPUGPU_MPI_SEND(id_gpumpi, Qm, 2*ny, rerank+1, 1, MPI_COMM_WORLD, ireq2(2), ierr)
      endif
      !print *, "myrank=", myrank, "send Qre"
    endif
    if (myrank == 0) then
      call CPUGPU_MPI_RECV(id_gpumpi, Qre, 5*ny*(nz-6), rerank+1, 0, MPI_COMM_WORLD, ireq, ierr)
    endif
  end subroutine step_rescale


  subroutine wait_rescale(myrank, ireq, ireq2, istat, istat2)
    integer, intent(in)    :: myrank
    integer, intent(inout) :: ireq, ireq2(2), istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    integer ierr
    if (myrank == rerank) then
      call MPI_WAITALL(2, ireq2, istat2, ierr)
      !print *, "myrank=", myrank, "WAITALL send Qm, Qre"
    endif
    if (myrank == 0) then
      call MPI_WAIT(ireq, istat, ierr)
      !print *, "myrank=", myrank, "WAIT recv Qre"
    endif
  end subroutine wait_rescale


  subroutine rescale_recv_send(num, flag_re, nx, ny, nz, step, y, Jacobian, Qm_cpu)
    integer, intent(in)    :: num
    integer, intent(inout) :: flag_re
    integer, intent(in)    :: nx, ny, nz, step
    real(8), intent(in)    :: y(ny), Jacobian(nx,ny)
    real(8), intent(inout) :: Qm_cpu(ny*2)
    real(8)         :: Qre_cpu(ny*(nz-6)*5), bltre
    real(8), device ::     Qre(ny*(nz-6)*5), Qm(ny*2)
    integer stat, errorcode, ierr, ireq, ireqs(2), istat(MPI_STATUS_SIZE), istats(MPI_STATUS_SIZE,2), j
    real(8) t
    character(len=40) filename
    write(filename, "(a)") "data/rescaling.d"

    if (num == 1) then
      call CPUGPU_MPI_RECV(id_gpumpi, Qre, 5*ny*(nz-6), rerank, 0, MPI_COMM_WORLD, ireqs(1), ierr)
      call CPUGPU_MPI_RECV(id_gpumpi, Qm,  2*ny,        rerank, 1, MPI_COMM_WORLD, ireqs(2), ierr)
      call MPI_WAITALL(2, ireqs, istats, ierr)
      stat = cudaDeviceSynchronize()
      stat = cudaMemcpy(Qm_cpu, Qm, 2*ny, cudaMemcpyDeviceToHost)
      if (stat /= cudaSuccess) then
        print *, "Qm  cudaMemcpy failed:", trim(cudaGetErrorString(stat))
      endif
    else
      call CPUGPU_MPI_RECV(id_gpumpi, Qre, 5*ny*(nz-6), rerank, 0, MPI_COMM_WORLD, ireq, ierr)
      call MPI_WAIT(ireq, istat, ierr)
      stat = cudaDeviceSynchronize()
    endif

    stat = cudaMemcpy(Qre_cpu, Qre, 5*ny*(nz-6), cudaMemcpyDeviceToHost)
    if (stat /= cudaSuccess) then
      print *, "Qre cudaMemcpy failed:", trim(cudaGetErrorString(stat))
    endif
    stat = cudaDeviceSynchronize()
    call set_rescale(flag_re, step, nx, ny, nz-6, y, Jacobian, Qm_cpu, bltre, Qre_cpu)
    stat = cudaMemcpy(Qre, Qre_cpu, 5*ny*(nz-6), cudaMemcpyHostToDevice)
    call CPUGPU_MPI_SEND(id_gpumpi, Qre, 5*ny*(nz-6), 0, 0, MPI_COMM_WORLD, ireq, ierr)
    call MPI_WAIT(ireq, istat, ierr)
    if (flag_re == 1) then
      call MPI_BCAST(flag_re, 1, MPI_INTEGER, rerank+1, MPI_COMM_WORLD, ierr)
    endif
    if (num == 1) then
      t = nt * step * dt
      if (flag_re >= 1 .and. step >= start_rescale) then
        open(10, file=filename, position="append")
        write(10, "(2e12.4, a)") t, bltre/blt, "rescale"
        close(10)
      else
        open(10, file=filename, position="append")
        write(10, "(2e12.4, a)") t, bltre/blt, "cyclic"
        close(10)
      endif
    endif
  end subroutine rescale_recv_send


  subroutine write_Qm(ny, y, Qm)
    integer, intent(in) :: ny
    real(8), intent(in) :: y(ny), Qm(ny*2)
    real(8) u, v
    character(len=40) filename
    integer j, jj
    write(filename, "(a)") "recal/Qm.d"
    open(10, file=filename, status="replace", action="write")
    write(10, "(a)") "y    u     v"
    do j = 1, ny
      jj = 2 * (j-1)
      u = Qm(jj+1); v = Qm(jj+2)
      write(10, "(3e12.4)") y(j), u, v
    enddo
    close(10)
    write(filename, "(a)") "recal/Qm.dat"
    open(10, file=filename, status="replace", action="write", form="unformatted", access="stream")
    write(10) Qm
    close(10)
  end subroutine write_Qm


  subroutine set_rescale(flag_re, step, nx, ny, nz, y, Jacobian, Qm, bltre, Qre)
    integer, intent(inout) :: flag_re
    integer, intent(in)    :: step, nx, ny, nz
    real(8), intent(in)    :: y(ny), Jacobian(nx,ny), Qm(ny*2)
    real(8), intent(out)   :: bltre
    real(8), intent(inout) :: Qre(ny*nz*5) ! Q / J
    integer i, j, jup, jdown, jj, k, kh, l, j_offset, k_offset, ierr, errorcode
    integer, dimension(ny) :: jj_y, jj_e
    real(8) t, dudy, bltup, bltdown, taure, utre, utin, utre_nu, utin_nu, beta, mu, nu, ady, ade, one_ady, one_ade 
    ! mean properties at rescaling plane
    real(8), dimension(ny)    :: Um, Vm, Tm
    ! fluctuating properties at rescaling plane
    real(8), dimension(ny,nz) :: ufre, vfre, wfre, Tfre, pfre
    real(8), dimension(ny)    :: ypre, ypin, etre, etin
    ! fluctuating properties at both inner and outer region
    real(8), dimension(ny,nz) :: ufin, vfin, wfin, Tfin, pfin
    real(8), dimension(ny,nz) :: ufout, vfout, wfout, Tfout, pfout
    ! mean properties at both inner and outer region
    real(8), dimension(ny)    :: Umin, Vmin, Tmin, Umout, Vmout, Tmout
    ! weighting function
    real(8), dimension(ny)    :: weight
    ! properties at rescaling plane
    real(8) ure, vre, wre, rhore, Tre, pre
    ! rescaled properties at inlet
    real(8) uin, vin, win, rhoin, Tin, pin
    real(8) :: over_rhow = R * Taw / p0
    ! cache
    real(8) :: u_tmp, v_tmp, p_tmp, T_tmp, weight_tmp, one_weight, Jacobian_tmp, over_rho
    real(8) :: over_tanh4, over_blt = 1.d0 / blt
    ! for exception
    real(8) :: blt_min = 0.5d0 * blt, blt_max = 2.d0 * blt
    do k = 1, nz
      k_offset = ny * 5 * (k-1)
      do j = 1, ny
        j_offset = 5 * (j-1)
        do l = 1, 5
          i = k_offset + j_offset + l
          Qre(i) = Qre(i) * Jacobian(nre2,j)
    enddo;enddo;enddo

    do j = 1, ny
      u_tmp    = Qm(2*(j-1)+1)
      v_tmp    = Qm(2*(j-1)+2)
      T_tmp    = Taw - rf * u_tmp**2 * 0.5d0 * over_Cp
      Um(j)    = u_tmp
      Umin(j)  = u_tmp
      Umout(j) = u_tmp
      Vm(j)    = v_tmp
      Vmin(j)  = v_tmp
      Vmout(j) = v_tmp
      Tm(j)    = T_tmp
      Tmin(j)  = T_tmp
      Tmout(j) = T_tmp
    enddo

    ! blt up
    jup   = -1
    bltup = 0.d0
    do j = 2, ny
      if (Um(j) >= 0.99d0 * u0) then
        dudy  = (0.99d0 * u0 - Um(j-1)) / (-Um(j-1) + Um(j) + 1.d-15)
        bltup = y(j-1) + (-y(j-1) + y(j)) * dudy
        jup   = j
        exit
      endif
    enddo

    ! blt down
    jdown   = -1
    bltdown = 0.d0
    do j = ny-1, 1, -1
      if (Um(j) <= 0.99d0 * u0) then
        dudy    = (Um(j+1) - 0.99d0 * u0) / (-Um(j) + Um(j+1) + 1.d-15)
        bltdown = y(j+1) - (-y(j) + y(j+1)) * dudy
        jdown   = j
        exit
      endif
    enddo
    bltre = 0.5d0 * (bltup + bltdown)

    if (jup < 0 .or. jdown < 0) then
      print *, "No 99% doundary layer thickness: bltup=", bltup, ", bltdown=", bltdown
      call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
    endif
    if (flag_re >= 1 .and. step >= start_rescale) then
      if (abs(jup - jdown) > 1) then
        print *, "jup=", jup, ", jdown=", jdown, ", bltup=", bltup, ", bltdown=", bltdown
      endif
      if (bltre < blt_min) bltre = blt_min
      if (bltre > blt_max) bltre = blt_max
    endif

    if (bltre > blt) then
      flag_re = flag_re + 1
    endif

    if (flag_re >= 1 .and. step >= start_rescale) then
      ufin(:,:)  = 0.d0
      vfin(:,:)  = 0.d0
      wfin(:,:)  = 0.d0
      Tfin(:,:)  = 0.d0
      pfin(:,:)  = 0.d0
      ufout(:,:) = 0.d0
      vfout(:,:) = 0.d0
      wfout(:,:) = 0.d0
      Tfout(:,:) = 0.d0
      pfout(:,:) = 0.d0
      do k = 1, nz
        k_offset = ny * 5 * (k-1)
        do j = 1, ny
          j_offset = 5 * (j-1)
          rhore = Qre(k_offset+j_offset+1)
          over_rho = 1.d0 / rhore
          ure   = Qre(k_offset+j_offset+2) * over_rho
          vre   = Qre(k_offset+j_offset+3) * over_rho
          wre   = Qre(k_offset+j_offset+4) * over_rho
          pre   = gamma_1 * (Qre(k_offset+j_offset+5) - 0.5d0 * rhore * (ure**2 + vre**2 + wre**2)) 
          Tre   = pre / (rhore * R)
          ufre(j,k) = ure - Um(j)
          vfre(j,k) = vre - Vm(j)
          wfre(j,k) = wre
          Tfre(j,k) = Tre - Tm(j)
          pfre(j,k) = pre - p0
      enddo;enddo

      ! friction velocity
      mu    = mu0_T0_S / (Taw + 111.d0) * (Taw * over_T0)**1.5
      nu    = mu * over_rhow
      taure = mu * abs(-Um(1) + Um(2)) / (-y(1) + y(2))
      utre  = sqrt(taure * over_rhow)
      beta  = (bltre * over_blt)**0.1
      utin  = beta * utre

      utin_nu = utin / nu
      utre_nu = utre / nu
      over_tanh4 = 1.d0 / tanh(4.d0)
      do j = 1, ny
        ypin(j) = y(j) * utin_nu
        ypre(j) = y(j) * utre_nu
        etin(j) = y(j) * over_blt
        etre(j) = y(j) / bltre
        weight(j) = min(1.d0, 0.5d0 * (1.d0 + tanh(4.d0 * (etin(j) - 0.2d0) / (0.6d0 * etin(j) + 0.2d0)) * over_tanh4))
      enddo

      jj_y(:) = -1
      do j = 1, ny
        do jj = 2, ny
          if (ypre(jj) > ypin(j)) then
            ady     = (-ypre(jj-1) + ypin(j)) / (-ypre(jj-1) + ypre(jj))
            one_ady = 1.d0 - ady
            Umin(j) = beta * (one_ady * Um(jj-1) + ady * Um(jj))!Um(jj-1) + ady * (-Um(jj-1) + Um(jj))
            Vmin(j) =         one_ady * Vm(jj-1) + ady * Vm(jj) !Vm(jj-1) + ady * (-Vm(jj-1) + Vm(jj))
            Tmin(j) =         one_ady * Tm(jj-1) + ady * Tm(jj) !Tm(jj-1) + ady * (-Tm(jj-1) + Tm(jj))
            jj_y(j) = jj
            exit
          endif
      enddo;enddo

      jj_e(:) = -1
      do j = 1, ny
        do jj = 2, ny
          if (etre(jj) > etin(j)) then
            ade     = (-etre(jj-1) + etin(j)) / (-etre(jj-1) + etre(jj))
            one_ade = 1.d0 - ade
            Umout(j) = beta * (one_ade * Um(jj-1) + ade * Um(jj)) + (1.d0 - beta) * u0
            Vmout(j) =         one_ade * Vm(jj-1) + ade * Vm(jj) !Vm(jj-1) + ade * (-Vm(jj-1) + Vm(jj))
            Tmout(j) =         one_ade * Tm(jj-1) + ade * Tm(jj) !Tm(jj-1) + ade * (-Tm(jj-1) + Tm(jj))
            jj_e(j)  = jj
            exit
          endif
      enddo;enddo
      
      do k = 1, nz
        do j = 1, ny
          jj = jj_y(j)
          if (jj > 0) then
            ady     = (-ypre(jj-1) + ypin(j)) / (-ypre(jj-1) + ypre(jj))
            one_ady = 1.d0 - ady
            ufin(j,k) = beta * (one_ady * ufre(jj-1,k) + ady * ufre(jj,k))!ufre(jj-1,k) + ady * (-ufre(jj-1,k) + ufre(jj,k))
            vfin(j,k) = beta * (one_ady * vfre(jj-1,k) + ady * vfre(jj,k))!vfre(jj-1,k) + ady * (-vfre(jj-1,k) + vfre(jj,k))
            wfin(j,k) = beta * (one_ady * wfre(jj-1,k) + ady * wfre(jj,k))!wfre(jj-1,k) + ady * (-wfre(jj-1,k) + wfre(jj,k))
            Tfin(j,k) =         one_ady * Tfre(jj-1,k) + ady * Tfre(jj,k) !Tfre(jj-1,k) + ady * (-Tfre(jj-1,k) + Tfre(jj,k))
            pfin(j,k) =         one_ady * pfre(jj-1,k) + ady * pfre(jj,k) !pfre(jj-1,k) + ady * (-pfre(jj-1,k) + pfre(jj,k))
          endif
          jj = jj_e(j)
          if (jj > 0) then
            ade     = (-etre(jj-1) + etin(j)) / (-etre(jj-1) + etre(jj))
            one_ade = 1.d0 - ade
            ufout(j,k) = beta * (one_ade * ufre(jj-1,k) + ade * ufre(jj,k))!ufre(jj-1,k) + ade * (-ufre(jj-1,k) + ufre(jj,k))
            vfout(j,k) = beta * (one_ade * vfre(jj-1,k) + ade * vfre(jj,k))!vfre(jj-1,k) + ade * (-vfre(jj-1,k) + vfre(jj,k))
            wfout(j,k) = beta * (one_ade * wfre(jj-1,k) + ade * wfre(jj,k))!wfre(jj-1,k) + ade * (-wfre(jj-1,k) + wfre(jj,k))
            Tfout(j,k) =         one_ade * Tfre(jj-1,k) + ade * Tfre(jj,k) !Tfre(jj-1,k) + ade * (-Tfre(jj-1,k) + Tfre(jj,k))
            pfout(j,k) =         one_ade * pfre(jj-1,k) + ade * pfre(jj,k) !pfre(jj-1,k) + ade * (-pfre(jj-1,k) + pfre(jj,k))
          endif
      enddo;enddo
  
      ! re-introducing
      do k = 1, nz
        kh = mod(k+nz/2,nz) + 1
        k_offset = ny * 5 * (k-1)
        do j = 1, ny
          j_offset     = 5 * (j-1)
          weight_tmp   = weight(j)
          one_weight   = 1.d0 - weight_tmp
          Jacobian_tmp = 1.d0 / Jacobian(nre2,j)
          uin = (Umin(j) + ufin(j,kh)) * one_weight + (Umout(j) + ufout(j,kh)) * weight_tmp
          vin = (Vmin(j) + vfin(j,kh)) * one_weight + (Vmout(j) + vfout(j,kh)) * weight_tmp
          win =            wfin(j,kh)  * one_weight +             wfout(j,kh)  * weight_tmp
          Tin = (Tmin(j) + Tfin(j,kh)) * one_weight + (Tmout(j) + Tfout(j,kh)) * weight_tmp
          pin = (p0      + pfin(j,kh)) * one_weight + (p0       + pfout(j,kh)) * weight_tmp
          rhoin = pin / (R * Tin)
          Qre(k_offset+j_offset+1) = rhoin * Jacobian_tmp
          Qre(k_offset+j_offset+2) = rhoin * uin * Jacobian_tmp
          Qre(k_offset+j_offset+3) = rhoin * vin * Jacobian_tmp
          Qre(k_offset+j_offset+4) = rhoin * win * Jacobian_tmp
          Qre(k_offset+j_offset+5) = (pin * over_gamma_1 + 0.5d0 * rhoin * (uin**2 + vin**2 + win**2)) * Jacobian_tmp
      enddo;enddo
    else
      ! cyclic boundary condition !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do k = 1, nz
        k_offset = ny * 5 * (k-1)
        do j = 1, ny
          j_offset = 5 * (j-1)
          Jacobian_tmp = 1.d0 / Jacobian(nre2,j)
          do l = 1, 5
            Qre(k_offset+j_offset+l) = Qre(k_offset+j_offset+l) * Jacobian_tmp
      enddo;enddo;enddo
    endif
  end subroutine set_rescale
end module calc_rescale

