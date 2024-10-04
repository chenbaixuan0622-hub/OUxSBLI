module calc_rescale
  use mod_globals, only : nt, dt, gamma , R, u0, strat_rescale
contains
  subroutine set_rescale(step,nx,ny,nz,nre,blt,y,Jacobian,Qre)
    integer, intent(in)    :: step, nx, ny, nz, nre
    real(8), intent(in)    :: blt
    real(8), intent(in)    :: y(ny)
    real(8), intent(in)    :: Jacobian(nx,ny,nz)
    real(8), intent(inout) :: Qre(10,ny,nz,5) ! Q / J
    integer i, j, jj, k, l
    real(8) :: mu0 = 1.716d-5, T0 = 273.2d0, S = 111.d0
    real(8) t, bltre1, bltre2, bltre, taure, utre, utin, beta, mu, nu, ady, ade 
    ! mean properties at rescaling plane
    real(8), dimension(ny)    :: Um, Vm, Wm, pm, Tm
    ! fluctuating properties at rescaling plane
    real(8), dimension(ny,nz) :: ufre, vfre, wfre, pfre, Tfre
    real(8), dimension(ny)    :: ypre, ypin, etre, etin
    ! fluctuating properties at both inner and outer region
    real(8), dimension(ny,nz) :: ufin, vfin, wfin, pfin, Tfin, ufout, vfout, wfout, pfout, Tfout
    ! mean properties at both inner and outer region
    real(8), dimension(ny,nz) :: Umin, Vmin, Wmin, pmin, Tmin, Umout, Vmout, Wmout, pmout, Tmout
    ! weighting function
    real(8), dimension(ny)    :: weight
    ! properties at rescaling plane
    real(8) ure, vre, wre, pre, Tre, rhore
    ! rescaled properties at inlet
    real(8) uin, vin, win, pin, Tin, rhoin
    character(len=40) filename
    write(filename, "(a)") "data/rescaling.d"

    t = nt * step * dt

    do l = 1, 5
      do k = 1, nz
        do j = 1, ny
          do i = 1, 10
            Qre(i,j,k,l) = Qre(i,j,k,l) * Jacobian(nre+i,j,k)
    enddo;enddo;enddo;enddo

    do j = 1, ny
      rho   = sum(Qre(:,j,:,1)) / dble(size(Qre(:,j,:,1)))
      Um(j) = sum(Qre(:,j,:,2) / Qre(:,j,:,1)) / dble(size(Qre(:,j,:,1)))
      Vm(j) = sum(Qre(:,j,:,3) / Qre(:,j,:,1)) / dble(size(Qre(:,j,:,1)))
      Wm(j) = sum(Qre(:,j,:,4) / Qre(:,j,:,1)) / dble(size(Qre(:,j,:,1)))
      pm(j) = sum((gamma - 1.d0) * (Qre(:,j,:,5) &
              - 0.5d0 * (Qre(:,j,:,2)**2 + Qre(:,j,:,3)**2 + Qre(:,j,:,4)**2) / Qre(:,j,:,1))) / dble(size(Qre(:,j,:,1)))
      Tm(j) = pm(j) / (R * rho)
    enddo

    ! check boundary layer thickness at rescaling plane
    bltre = 0.d0
    do j = 2, ny
      if (Um(j) >= 0.99d0 * u0) then
        bltre = y(j) - (-y(j-1) + y(j)) * (Um(j) - 0.99d0 * u0) / (-Um(j-1) + Um(j) + 1.d-20)
        exit
      endif
    enddo

    if (bltre > blt .and. step >= start_rescale) then
      ! rescaling !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! calc fluctuating part   u'(x,y,z,t) = u(x,y,z,t) - U(x,y)
      ! U(x,y) average velocity in the spanwise direction and time
      open(10, file=filename, position="append")
      write(10, "(2e12.4, a)") t, bltre, "rescale"
      close(10)
      do k = 1, nz
        do j = 1, ny
          rhore = Qre(1,j,k,1)
          ure   = Qre(1,j,k,2) / rhore
          vre   = Qre(1,j,k,3) / rhore
          wre   = Qre(1,j,k,4) / rhore
          pre   = (gamma - 1.d0) * (Qre(1,j,k,5) - 0.5d0 * rhore * (ure**2 + vre**2 + wre**2)) 
          Tre   = pre / (rhore * R)
          ufre(j,k) = ure - Um(j)
          vfre(j,k) = vre - Vm(j)
          wfre(j,k) = wre - Wm(j)
          pfre(j,k) = pre - pm(j)
          Tfre(j,k) = Tre - Tm(j)
      enddo;enddo

      ! friction velocity
      rhore = pm(1) / (R * Tm(1))
      mu    = mu0 * ((T0 + S) / (Tm(1) + S)) * (Tm(1) / T0)**1.5
      nu    = mu / rhore
      taure = mu * abs(-Um(1) + Um(2)) / (-y(1) + y(2))
      utre  = sqrt(taure / rhore)
      utin  = utre * (bltre / blt)**0.1
      beta  = utin / utre

      do j = 1, ny
        ypin(j) = y(j) * utin / nu
        ypre(j) = y(j) * utre / nu
        etin(j) = y(j) / blt
        etre(j) = y(j) / bltre
      enddo

      ! set temporal values to avoid NaN
      do k = 1, nz
        do j = 1, ny
          ! inner region
          ! mean
          Umin(j,k) = Um(j)
          Vmin(j,k) = Vm(j)
          Wmin(j,k) = Wm(j)
          pmin(j,k) = pm(j)
          Tmin(j,k) = Tm(j)
          ! fluctuating
          ufin(j,k) = 0.d0
          vfin(j,k) = 0.d0
          wfin(j,k) = 0.d0
          pfin(j,k) = 0.d0
          Tfin(j,k) = 0.d0
          ! outer region
          ! mean
          Umout(j,k) = Um(j)
          Vmout(j,k) = Vm(j)
          Wmout(j,k) = Wm(j)
          pmout(j,k) = pm(j)
          Tmout(j,k) = Tm(j)
          ! fluctuating
          ufout(j,k) = 0.d0
          vfout(j,k) = 0.d0
          wfout(j,k) = 0.d0
          pfout(j,k) = 0.d0
          Tfout(j,k) = 0.d0
      enddo;enddo

      do j = 1, ny
        do jj = 2, ny
          if (ypre(jj) > ypin(j)) then
            ady = (-ypre(jj-1) + ypin(j)) / (-ypre(jj-1) + ypre(jj))
            do k = 1, nz
              ! mean
              Umin(j,k) = beta * (Um(jj-1) + ady * (-Um(jj-1) + Um(jj)))
              Vmin(j,k) =         Vm(jj-1) + ady * (-Vm(jj-1) + Vm(jj))
              Wmin(j,k) = 0.d0
              pmin(j,k) =         pm(jj-1) + ady * (-pm(jj-1) + pm(jj)) 
              Tmin(j,k) =         Tm(jj-1) + ady * (-Tm(jj-1) + Tm(jj))
              ! fluctuating
              ufin(j,k) = beta * (ufre(jj-1,k) + ady * (-ufre(jj-1,k) + ufre(jj,k)))
              vfin(j,k) = beta * (vfre(jj-1,k) + ady * (-vfre(jj-1,k) + vfre(jj,k)))
              wfin(j,k) = beta * (wfre(jj-1,k) + ady * (-wfre(jj-1,k) + wfre(jj,k)))
              pfin(j,k) =         pfre(jj-1,k) + ady * (-pfre(jj-1,k) + pfre(jj,k))
              Tfin(j,k) =         Tfre(jj-1,k) + ady * (-Tfre(jj-1,k) + Tfre(jj,k))
            enddo
            exit
          endif
      enddo;enddo

      do j = 1, ny
        do jj = 2, ny
          if (etre(jj) > etin(j)) then
            ade = (-etre(jj-1) + etin(j)) / (-etre(jj-1) + etre(jj))
            do k = 1, nz
              ! mean
              Umout(j,k) = beta * (Um(jj-1) + ade * (-Um(jj-1) + Um(jj))) + (1.d0 - beta) * u0
              Vmout(j,k) =         Vm(jj-1) + ade
              Wmout(j,k) = 0.d0
              pmout(j,k) =         pm(jj-1) + ade * (-pm(jj-1) + pm(jj))
              Tmout(j,k) =         Tm(jj-1) + ade * (-Tm(jj-1) + Tm(jj))
              ! fluctuating
              ufout(j,k) = beta * (ufre(jj-1,k) + ade * (-ufre(jj-1,k) + ufre(jj,k)))
              vfout(j,k) = beta * (vfre(jj-1,k) + ade * (-vfre(jj-1,k) + vfre(jj,k)))
              wfout(j,k) = beta * (wfre(jj-1,k) + ade * (-wfre(jj-1,k) + wfre(jj,k)))
              pfout(j,k) =         pfre(jj-1,k) + ade * (-pfre(jj-1,k) + pfre(jj,k))
              Tfout(j,k) =         Tfre(jj-1,k) + ade * (-Tfre(jj-1,k) + Tfre(jj,k))
            enddo
            exit
          endif
      enddo;enddo
      
      ! weighting function
      do j = 1, ny
        weight(j) = min(1.d0, 0.5d0 * (1.d0 + tanh(4.d0 * (etin(j) - 0.2d0) / ((1.d0 - 0.4d0) * etin(j) + 0.2d0)) / tanh(4.d0)))
      enddo

      ! re-introducing
      do k = 2, nz-1
        do j = 1, ny
          uin   = (Umin(j,k) + ufin(j,k)) * (1.d0 - weight(j)) + (Umout(j,k) + ufout(j,k)) * weight(j)
          vin   = (Vmin(j,k) + vfin(j,k)) * (1.d0 - weight(j)) + (Vmout(j,k) + vfout(j,k)) * weight(j)
          win   = (Wmin(j,k) + wfin(j,k)) * (1.d0 - weight(j)) + (Wmout(j,k) + wfout(j,k)) * weight(j)
          pin   = (pmin(j,k) + pfin(j,k)) * (1.d0 - weight(j)) + (pmout(j,k) + pfout(j,k)) * weight(j)
          Tin   = (Tmin(j,k) + Tfin(j,k)) * (1.d0 - weight(j)) + (Tmout(j,k) + Tfout(j,k)) * weight(j)
          rhoin = pin / (R * Tin)
          Qre(1,j,k,1) = rhoin / Jacobian(nre,j,k)
          Qre(1,j,k,2) = rhoin * uin / Jacobian(nre,j,k)
          Qre(1,j,k,3) = rhoin * vin / Jacobian(nre,j,k)
          Qre(1,j,k,4) = rhoin * win / Jacobian(nre,j,k)
          Qre(1,j,k,5) = (pin / (gamma - 1.d0) + 0.5d0 * rhoin * (uin**2 + vin**2 + win**2)) / Jacobian(nre,j,k)
      enddo;enddo
    else
      open(10, file=filename, position="append")
      write(10, "(2e12.4, a)") t, bltre, "cyclic"
      close(10)
      ! cyclic boundary condition !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do l = 1, 5
        do k = 1, nz
          do j = 1, ny
            ! inlet
            Qre(1,j,k,l) = Qre(1,j,k,l) / Jacobian(nre,j,k)
      enddo;enddo;enddo
    endif
  end subroutine set_rescale
end module calc_rescale

