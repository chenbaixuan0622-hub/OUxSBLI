module calc_rescale
  use mpi
  use mod_globals, only : nre1, nre2, rerank, nt, dt, gamma , R, Pr, u0, rho0, p0, M0, blt, start_rescale
contains
  subroutine calc_mean(nx, ny, nz, Jacobian, QJ, Qm_cpu)
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(in), device :: Jacobian(ny), QJ(nx,ny,nz,5)
    real(8), intent(out)        :: Qm_cpu(ny,5)
    real(8) Q1, Q2, Q3, Q4, Q5
    real(8), device :: Qm(ny,5)
    integer i, k
    !$cuf kernel do <<<*,*>>>
    do j = 1, ny
      Q1 = 0.d0
      Q2 = 0.d0
      Q3 = 0.d0
      Q4 = 0.d0
      Q5 = 0.d0
      do k = 4, nz-3
        do i = nre1, nre2
          Q1 = Q1 + QJ(i,j,k,1) * Jacobian(j)
          Q2 = Q2 + QJ(i,j,k,2) / QJ(i,j,k,1)
          Q3 = Q3 + QJ(i,j,k,3) / QJ(i,j,k,1)
          Q4 = Q4 + QJ(i,j,k,4) / QJ(i,j,k,1)
          Q5 = Q5 + (gamma - 1.d0) * Jacobian(j) * (QJ(i,j,k,5) &
                    - 0.5d0 * (QJ(i,j,k,2)**2 + QJ(i,j,k,3)**2 + QJ(i,j,k,4)**2) / QJ(i,j,k,1))
      enddo;enddo
      Qm(j,1) = Q1 / dble((nre2 - nre1 + 1) * (nz - 6))
      Qm(j,2) = Q2 / dble((nre2 - nre1 + 1) * (nz - 6))
      Qm(j,3) = Q3 / dble((nre2 - nre1 + 1) * (nz - 6))
      Qm(j,4) = Q4 / dble((nre2 - nre1 + 1) * (nz - 6))
      Qm(j,5) = Q5 / dble((nre2 - nre1 + 1) * (nz - 6))
    enddo
    Qm_cpu = Qm
  end subroutine calc_mean

  subroutine copy(nx, ny, nz, QJ, Qre)
    integer, intent(in)          :: nx, ny, nz
    real(8), intent(in), device  :: QJ(nx,ny,nz,5)
    real(8), intent(out), device :: Qre(ny,nz-6,5)
    integer j, k, l
    !$cuf kernel do <<<*,*>>>
    do l = 1, 5
      do k = 4, nz-3
        do j = 1, ny
          Qre(j,k-3,l) = QJ(nre2,j,k,l)
    enddo;enddo;enddo
  end subroutine copy

  subroutine rescale_recv_send(nx, ny, nz, step, y, Jacobian)
    integer, intent(in) :: nx, ny, nz, step
    real(8), intent(in) :: y(ny), Jacobian(ny)
    real(8) Qre(ny,nz-6,5), Qm(ny,5)
    integer ierr, ireqs(2), istats(MPI_STATUS_SIZE,2)

    call MPI_IRECV(Qre, 5*ny*(nz-6), MPI_REAL8, rerank, 0, MPI_COMM_WORLD, ireqs(1), ierr)
    call MPI_IRECV(Qm,  5*ny,        MPI_REAL8, rerank, 1, MPI_COMM_WORLD, ireqs(2), ierr)
    call MPI_WAITALL(2, ireqs, istats, ierr)
    call set_rescale(step, nx, ny, nz-6, y, Jacobian, Qm, Qre)
    call MPI_SEND(Qre, 5*ny*(nz-6), MPI_REAL8, 0, 0, MPI_COMM_WORLD, ierr)
  end subroutine rescale_recv_send

  subroutine set_rescale(step, nx, ny, nz, y, Jacobian, Qm, Qre)
    integer, intent(in)    :: step, nx, ny, nz! nz-6
    real(8), intent(in)    :: y(ny), Jacobian(ny), Qm(ny,5)
    real(8), intent(inout) :: Qre(ny,nz,5) ! Q / J
    integer i, j, jj, k, kh, l
    real(8) :: mu0 = 1.716d-5, T0 = 273.2d0, S = 111.d0, Cp = gamma * R / (gamma - 1.d0)
    real(8) t, dudy, blt1, blt2, bltre, taure, utre, utin, beta, mu, nu, ady, ade 
    ! mean properties at rescaling plane
    real(8), dimension(ny)    :: Um, Vm, Wm, rhom, Tm, pm
    ! fluctuating properties at rescaling plane
    real(8), dimension(ny,nz) :: ufre, vfre, wfre, rhofre, Tfre, pfre
    real(8), dimension(ny)    :: ypre, ypin, etre, etin
    ! fluctuating properties at both inner and outer region
    real(8), dimension(ny,nz) :: ufin,  vfin,  wfin,  rhofin,  Tfin,   pfin,   ufout,  vfout,  wfout, rhofout, Tfout, pfout
    real(8)                   :: ufins, vfins, wfins, rhofins, ufouts, vfouts, wfouts, rhofouts
    ! mean properties at both inner and outer region
    real(8), dimension(ny)    :: Umin, Vmin, rhomin, pmin, Tmin, Umout, Vmout, rhomout, pmout, Tmout
    ! weighting function
    real(8), dimension(ny)    :: weight
    ! properties at rescaling plane
    real(8) ure, vre, wre, rhore, Tre, pre
    ! rescaled properties at inlet
    real(8) uin, vin, win, rhoin, Tin, pin
    character(len=40) filename
    write(filename, "(a)") "data/rescaling.d"

    t = nt * step * dt

    do l = 1, 5
      do k = 1, nz
        do j = 1, ny
          Qre(j,k,l) = Qre(j,k,l) * Jacobian(j)
    enddo;enddo;enddo

    do j = 1, ny
      rhom(j) = Qm(j,1)
        Um(j) = Qm(j,2)
        Vm(j) = Qm(j,3)
        Wm(j) = Qm(j,4)
        pm(j) = Qm(j,5)
        Tm(j) = pm(j) / (R * rhom(j))
    enddo

    ! check boundary layer thickness at rescaling plane
    bltre = 0.d0
    do j = 2, ny
      if (Um(j) >= 0.99d0 * u0) then
        dudy  = (Um(j) - 0.99d0 * u0) / (-Um(j-1) + Um(j))
        bltre =    y(j) - (-   y(j-1) +    y(j)) * dudy
        blt1  =  y(j-1)
        blt2  =    y(j)
        if (bltre < blt1 .or. blt2 < bltre) then
          bltre = blt1
        endif
        exit
      endif
    enddo

    if (bltre > blt .and. step >= start_rescale) then
      ! rescaling !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! calc fluctuating part   u'(x,y,z,t) = u(x,y,z,t) - U(x,y)
      ! U(x,y) average velocity in the spanwise direction and time
      open(10, file=filename, position="append")
      write(10, "(4e12.4, a)") t*1d3, blt1, bltre, blt2, "rescale"
      close(10)
      do k = 1, nz
        do j = 1, ny
          rhore = Qre(j,k,1)
          ure   = Qre(j,k,2) / rhore
          vre   = Qre(j,k,3) / rhore
          wre   = Qre(j,k,4) / rhore
          pre   = (gamma - 1.d0) * (Qre(j,k,5) - 0.5d0 * rhore * (ure**2 + vre**2 + wre**2)) 
          Tre   = pre / (rhore * R)
          ufre(j,k)   = ure   -   Um(j)
          vfre(j,k)   = vre   -   Vm(j)
          wfre(j,k)   = wre   -   Wm(j)
          rhofre(j,k) = rhore - rhom(j)
          Tfre(j,k)   = Tre   -   Tm(j)
          pfre(j,k)   = pre   -   pm(j)
          ! set initial values to avoid NaN
          ! inner region
          ! fluctuating
          ufin(j,k)   = 0.d0
          vfin(j,k)   = 0.d0
          wfin(j,k)   = 0.d0
          rhofin(j,k) = 0.d0
          Tfin(j,k)   = 0.d0
          pfin(j,k)   = 0.d0
          ! outer region
          ! fluctuating
          ufout(j,k)   = 0.d0
          vfout(j,k)   = 0.d0
          wfout(j,k)   = 0.d0
          rhofout(j,k) = 0.d0
          Tfout(j,k)   = 0.d0
          pfout(j,k)   = 0.d0
      enddo;enddo

      ! friction velocity
      mu    = mu0 * ((T0 + S) / (Tm(1) + S)) * (Tm(1) / T0)**1.5
      nu    = mu / rhom(1)
      taure = mu * abs(-Um(1) + Um(2)) / (-y(1) + y(2))
      utre  = sqrt(taure / rhom(1))
      beta  = (bltre / blt)**0.1
      utin  = beta * utre

      do j = 1, ny
        ypin(j) = y(j) * utin / nu
        ypre(j) = y(j) * utre / nu
        etin(j) = y(j) / blt
        etre(j) = y(j) / bltre
        ! mean
        Umin(j)    = Um(j)
        Vmin(j)    = Vm(j)
        rhomin(j)  = rhom(j)
        Tmin(j)    = Tm(j)
        pmin(j)    = pm(j)
        Umout(j)   = Um(j)
        Vmout(j)   = Vm(j)
        rhomout(j) = rhom(j)
        Tmout(j)   = Tm(j)
        pmout(j)   = pm(j)
      enddo

      do j = 1, ny
        do jj = 2, ny
          if (ypre(jj) > ypin(j)) then
            ady = (-ypre(jj-1) + ypin(j)) / (-ypre(jj-1) + ypre(jj))
            ! mean
            Umin(j)   = beta * (Um(jj-1) + ady * (-  Um(jj-1) +   Um(jj)))
            Vmin(j)   =         Vm(jj-1) + ady * (-  Vm(jj-1) +   Vm(jj))
            rhomin(j) =       rhom(jj-1) + ady * (-rhom(jj-1) + rhom(jj)) 
            Tmin(j)   =         Tm(jj-1) + ady * (-  Tm(jj-1) +   Tm(jj))
            pmin(j)   =         pm(jj-1) + ady * (-  pm(jj-1) +   pm(jj))
            do k = 1, nz
              ! fluctuating
              ufin(j,k)   = beta * (ufre(jj-1,k) + ady * (-  ufre(jj-1,k) +   ufre(jj,k)))
              vfin(j,k)   = beta * (vfre(jj-1,k) + ady * (-  vfre(jj-1,k) +   vfre(jj,k)))
              wfin(j,k)   = beta * (wfre(jj-1,k) + ady * (-  wfre(jj-1,k) +   wfre(jj,k)))
              rhofin(j,k) =       rhofre(jj-1,k) + ady * (-rhofre(jj-1,k) + rhofre(jj,k))
              Tfin(j,k)   =         Tfre(jj-1,k) + ady * (-  Tfre(jj-1,k) +   Tfre(jj,k))
              pfin(j,k)   =         pfre(jj-1,k) + ady * (-  pfre(jj-1,k) +   pfre(jj,k))
            enddo
            exit
          endif
      enddo;enddo

      do j = 1, ny
        do jj = 2, ny
          if (etre(jj) > etin(j)) then
            ade = (-etre(jj-1) + etin(j)) / (-etre(jj-1) + etre(jj))
            ! mean
            Umout(j)   = beta * (Um(jj-1) + ade * (-  Um(jj-1) +   Um(jj))) + (1.d0 - beta) * u0
            Vmout(j)   =         Vm(jj-1) + ade * (-  Vm(jj-1) +   Vm(jj))
            rhomout(j) =       rhom(jj-1) + ade * (-rhom(jj-1) + rhom(jj))
            Tmout(j)   =         Tm(jj-1) + ade * (-  Tm(jj-1) +   Tm(jj))
            pmout(j)   =         pm(jj-1) + ade * (-  pm(jj-1) +   pm(jj))
            do k = 1, nz
              ! fluctuating
              ufout(j,k)   = beta * (ufre(jj-1,k) + ade * (-  ufre(jj-1,k) +   ufre(jj,k)))
              vfout(j,k)   = beta * (vfre(jj-1,k) + ade * (-  vfre(jj-1,k) +   vfre(jj,k)))
              wfout(j,k)   = beta * (wfre(jj-1,k) + ade * (-  wfre(jj-1,k) +   wfre(jj,k)))
              rhofout(j,k) =       rhofre(jj-1,k) + ade * (-rhofre(jj-1,k) + rhofre(jj,k))
              Tfout(j,k)   =         Tfre(jj-1,k) + ade * (-  Tfre(jj-1,k) +   Tfre(jj,k))
              pfout(j,k)   =         pfre(jj-1,k) + ade * (-  pfre(jj-1,k) +   pfre(jj,k))
            enddo
            exit
          endif
      enddo;enddo
      
      ! weighting function
      do j = 1, ny
        weight(j) = min(1.d0, 0.5d0 * (1.d0 + tanh(4.d0 * (etin(j) - 0.2d0) / ((1.d0 - 0.4d0) * etin(j) + 0.2d0)) / tanh(4.d0)))
      enddo
  
      ! re-introducing
      do k = 1, nz
        kh = mod(k+nz/2,nz) + 1
        do j = 1, ny
          uin = (Umin(j) + ufin(j,kh)) * (1.d0 - weight(j)) + (Umout(j) + ufout(j,kh)) * weight(j)
          vin = (Vmin(j) + vfin(j,kh)) * (1.d0 - weight(j)) + (Vmout(j) + vfout(j,kh)) * weight(j)
          win =            wfin(j,kh)  * (1.d0 - weight(j)) +             wfout(j,kh)  * weight(j)
          Tin = (Tmin(j) + Tfin(j,kh)) * (1.d0 - weight(j)) + (Tmout(j) + Tfout(j,kh)) * weight(j)
          pin = (pmin(j) + pfin(j,kh)) * (1.d0 - weight(j)) + (pmout(j) + pfout(j,kh)) * weight(j)
          rhoin      = pin / (R * Tin)
          Qre(j,k,1) = rhoin / Jacobian(j)
          Qre(j,k,2) = rhoin * uin / Jacobian(j)
          Qre(j,k,3) = rhoin * vin / Jacobian(j)
          Qre(j,k,4) = rhoin * win / Jacobian(j)
          Qre(j,k,5) = (pin / (gamma - 1.d0) + 0.5d0 * rhoin * (uin**2 + vin**2 + win**2)) / Jacobian(j)
      enddo;enddo
    else
      open(10, file=filename, position="append")
      write(10, "(4e12.4, a)") t*1d3, blt1, bltre, blt2, "cyclic"
      close(10)
      ! cyclic boundary condition !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do k = 1, nz
        do j = 1, ny
          Qre(j,k,1) = Qre(j,k,1) / Jacobian(j)
          Qre(j,k,2) = Qre(j,k,2) / Jacobian(j)
          Qre(j,k,3) = Qre(j,k,3) / Jacobian(j)
          Qre(j,k,4) = Qre(j,k,4) / Jacobian(j)
          Qre(j,k,5) = Qre(j,k,5) / Jacobian(j)
      enddo;enddo
    endif
  end subroutine set_rescale
end module calc_rescale

