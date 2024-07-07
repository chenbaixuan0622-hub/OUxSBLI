module calc_rescale
  use cudafor
  use mpi
  use mod_globals, only : gamma , R, u0
contains
  subroutine calc_mean(step,nx,ny,nz,Q,Um,Vm,Wm,pm,Tm)
    integer, intent(in)                       :: step, nx, ny, nz
    real(8), intent(in), dimension(2,ny,nz,5) :: Q
    real(8), intent(inout), dimension(2,ny)   :: Um, Vm, Wm, pm, Tm
    integer i, j, k
    real(8) p, rho
    real(8), dimension(2,ny) :: usum, vsum, wsum, psum, Tsum
    usum(:,:) = 0.d0; vsum(:,:) = 0.d0; wsum(:,:) = 0.d0; psum(:,:) = 0.d0; Tsum(:,:) = 0.d0
    ! span-wise direction
    do k = 1, nz
      do j = 1, ny
        do i = 1, 2
          rho       = Q(i,j,k,1)
          usum(i,j) = usum(i,j) + Q(i,j,k,2) / rho 
          vsum(i,j) = vsum(i,j) + Q(i,j,k,3) / rho
          wsum(i,j) = wsum(i,j) + Q(i,j,k,4) / rho
          p         = (gamma - 1.d0) * (Q(i,j,k,5) &
                      - 0.5d0 * (Q(i,j,k,2)**2 + Q(i,j,k,3)**2 + Q(i,j,k,4)**2) / rho) 
          psum(i,j) = psum(i,j) + p
          Tsum(i,j) = Tsum(i,j) + p / (R * rho)
    enddo;enddo;enddo

    ! time direction
    do j = 1, ny
      do i = 1, 2
        Um(i,j) = ((dble(step) - 1.d0) * Um(i,j) + usum(i,j) / dble(nz)) / dble(step)
        Vm(i,j) = ((dble(step) - 1.d0) * Vm(i,j) + vsum(i,j) / dble(nz)) / dble(step)
        Wm(i,j) = ((dble(step) - 1.d0) * Wm(i,j) + wsum(i,j) / dble(nz)) / dble(step)
        pm(i,j) = ((dble(step) - 1.d0) * pm(i,j) + psum(i,j) / dble(nz)) / dble(step)
        Tm(i,j) = ((dble(step) - 1.d0) * Tm(i,j) + Tsum(i,j) / dble(nz)) / dble(step)
    enddo;enddo
  end subroutine calc_mean

  subroutine set_rescale(step,myrank,nx,ny,nz,nre,blt,y,Jacobian,Um,Vm,Wm,pm,Tm,Qre)
    integer, intent(in)                     :: step, myrank, nx, ny, nz, nre
    real(8), intent(in)                     :: blt
    real(8), intent(in)                     :: y(ny)
    real(8), intent(in)                     :: Jacobian(nx,ny,nz)
    real(8), intent(inout), dimension(2,ny) :: Um, Vm, Wm, pm, Tm
    real(8), intent(inout)                  :: Qre(2,ny,nz,5) ! Q / J
    integer i, j, jj, k, l, ierr, status(MPI_STATUS_SIZE)
    real(8) :: mu0 = 1.716d-5, T0 = 273.2d0, S = 111.d0
    real(8) bltre1, bltre2, bltre, taure, utre, utin, beta, mu, nu, ady, ade 
    ! mean properties at rescaling plane
    ! fluctuating properties at rescaling plane
    real(8), dimension(2,ny,nz) :: ufre, vfre, wfre, pfre, Tfre
    real(8), dimension(ny)      :: ypre, ypin, etre, etin
    ! fluctuating properties at both inner and outer region
    real(8), dimension(2,ny,nz) :: ufin, vfin, wfin, pfin, Tfin, ufout, vfout, wfout, pfout, Tfout
    ! mean properties at both inner and outer region
    real(8), dimension(2,ny,nz) :: Umin, Vmin, Wmin, pmin, Tmin, Umout, Vmout, Wmout, pmout, Tmout
    ! weighting function
    real(8), dimension(ny)      :: weight
    ! properties at rescaling plane
    real(8) ure, vre, wre, pre, Tre, rhore
    ! rescaled properties at inlet
    real(8) uin, vin, win, pin, Tin, rhoin

    do l = 1, 5
      do k = 1, nz
        do j = 1, ny
          do i = 1, 2
            Qre(i,j,k,l) = Qre(i,j,k,l) * Jacobian(nre+i,j,k)
    enddo;enddo;enddo;enddo

    call calc_mean(step,nx,ny,nz,Qre,Um,Vm,Wm,pm,Tm)

    ! check boundary layer thickness at rescaling plane
    if (myrank == 3) then
      do j = 2, ny
        if (Um(1,j) >= 0.99d0 * u0 .and. Um(2,j) >= 0.99d0 * u0) then
          bltre1 = y(j) - (-y(j-1) + y(j)) * (Um(1,j) - 0.99d0 * u0) / (-Um(1,j-1) + Um(1,j) + 1.d-20)
          bltre2 = y(j) - (-y(j-1) + y(j)) * (Um(2,j) - 0.99d0 * u0) / (-Um(2,j-1) + Um(2,j) + 1.d-20)
          ! ensure bltre is not NaN
          if (bltre1 == bltre1 .and. bltre2 == bltre2) then
            bltre = 0.5d0 * (bltre1 + bltre2)
            exit
          endif
        endif
      enddo
    endif

    if (bltre >= blt .and. myrank == 3 .and. step >= 1000) then
      print *, "rescale", " blt=", real(bltre)
      ! rescaling !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! calc fluctuating part   u'(x,y,z,t) = u(x,y,z,t) - U(x,y)
      ! U(x,y) average velocity in the spanwise direction and time
      do k = 1, nz
        do j = 1, ny
          do i = 1, 2 ! 2 rescaleing planes are required for 4th-order accuracy flux
            rhore = Qre(i,j,k,1)
            ure   = Qre(i,j,k,2) / rhore
            vre   = Qre(i,j,k,3) / rhore
            wre   = Qre(i,j,k,4) / rhore
            pre   = (gamma - 1.d0) * (Qre(i,j,k,5) - 0.5d0 * rhore * (ure**2 + vre**2 + wre**2)) 
            Tre   = pre / (rhore * R)
            ufre(i,j,k) = ure - Um(i,j)
            vfre(i,j,k) = vre - Vm(i,j)
            wfre(i,j,k) = wre - Wm(i,j)
            pfre(i,j,k) = pre - pm(i,j)
            Tfre(i,j,k) = Tre - Tm(i,j)
      enddo;enddo;enddo

      ! friction velocity
      rhore = pm(1,1) / (R * Tm(1,1))
      mu    = mu0 * ((T0 + S) / (Tm(1,1) + S)) * (Tm(1,1) / T0)**1.5
      nu    = mu / rhore
      taure = mu * abs(Um(1,2)) / (-y(1) + y(2))
      utre  = sqrt(taure / rhore)
      utin  = utre * (bltre / blt)**0.1
      beta  = utin / utre

      do j = 1, ny
        ypin(j) = y(j) * utin / nu
        ypre(j) = y(j) * utre / nu
        etin(j) = y(j) / blt
        etre(j) = y(j) / bltre
      enddo

      do j = 1, ny
        do jj = 2, ny
          if (ypre(jj) > ypin(j)) then
            ady = (-ypre(jj-1) + ypin(j)) / (-ypre(jj-1) + ypre(jj))
            do k = 1, nz
              do i = 1, 2
                ! mean
                Umin(i,j,k) = beta * (Um(i,jj-1) + ady * (-Um(i,jj-1) + Um(i,jj)))
                Vmin(i,j,k) =         Vm(i,jj-1) + ady * (-Vm(i,jj-1) + Vm(i,jj))
                Wmin(i,j,k) = 0.d0
                pmin(i,j,k) =         pm(i,jj-1) + ady * (-pm(i,jj-1) + pm(i,jj)) 
                Tmin(i,j,k) =         Tm(i,jj-1) + ady * (-Tm(i,jj-1) + Tm(i,jj))
                ! fluctuating
                ufin(i,j,k) = beta * (ufre(i,jj-1,k) + ady * (-ufre(i,jj-1,k) + ufre(i,jj,k)))
                vfin(i,j,k) = beta * (vfre(i,jj-1,k) + ady * (-vfre(i,jj-1,k) + vfre(i,jj,k)))
                wfin(i,j,k) = beta * (wfre(i,jj-1,k) + ady * (-wfre(i,jj-1,k) + wfre(i,jj,k)))
                pfin(i,j,k) =         pfre(i,jj-1,k) + ady * (-pfre(i,jj-1,k) + pfre(i,jj,k))
                Tfin(i,j,k) =         Tfre(i,jj-1,k) + ady * (-Tfre(i,jj-1,k) + Tfre(i,jj,k))
            enddo;enddo
            exit
          endif
      enddo;enddo

      do j = 1, ny
        do jj = 2, ny
          if (etre(jj) > etin(j)) then
            ade = (-etre(jj-1) + etin(j)) / (-etre(jj-1) + etre(jj))
            do k = 1, nz
              do i = 1, 2
                ! mean
                Umout(i,j,k) = beta * (Um(i,jj-1) + ade * (-Um(i,jj-1) + Um(i,jj))) + (1.d0 - beta) * u0
                Vmout(i,j,k) =         Vm(i,jj-1) + ade
                Wmout(i,j,k) = 0.d0
                pmout(i,j,k) =         pm(i,jj-1) + ade * (-pm(i,jj-1) + pm(i,jj))
                Tmout(i,j,k) =         Tm(i,jj-1) + ade * (-Tm(i,jj-1) + Tm(i,jj))
                ! fluctuating
                ufout(i,j,k) = beta * (ufre(i,jj-1,k) + ade * (-ufre(i,jj-1,k) + ufre(i,jj,k)))
                vfout(i,j,k) = beta * (vfre(i,jj-1,k) + ade * (-vfre(i,jj-1,k) + vfre(i,jj,k)))
                wfout(i,j,k) = beta * (wfre(i,jj-1,k) + ade * (-wfre(i,jj-1,k) + wfre(i,jj,k)))
                pfout(i,j,k) =         pfre(i,jj-1,k) + ade * (-pfre(i,jj-1,k) + pfre(i,jj,k))
                Tfout(i,j,k) =         Tfre(i,jj-1,k) + ade * (-Tfre(i,jj-1,k) + Tfre(i,jj,k))
            enddo;enddo
            exit
          endif
      enddo;enddo
      
      ! weighting function
      do j = 1, ny
        weight(j) = min(1.d0, 0.5d0 * (1.d0 + tanh(4.d0 * (etin(j) - 0.2d0) / ((1.d0 - 0.4d0) * etin(j) + 0.2d0)) / tanh(4.d0)))
      enddo

      ! re-introducing
      do k = 2, nz-1
        do j = 2, ny-1
          do i = 1, 2
            uin   = (Umin(i,j,k) + ufin(i,j,k)) * (1.d0 - weight(j)) + (Umout(i,j,k) + ufout(i,j,k)) * weight(j)
            vin   = (Vmin(i,j,k) + vfin(i,j,k)) * (1.d0 - weight(j)) + (Vmout(i,j,k) + vfout(i,j,k)) * weight(j)
            win   = (Wmin(i,j,k) + wfin(i,j,k)) * (1.d0 - weight(j)) + (Wmout(i,j,k) + wfout(i,j,k)) * weight(j)
            pin   = (pmin(i,j,k) + pfin(i,j,k)) * (1.d0 - weight(j)) + (pmout(i,j,k) + pfout(i,j,k)) * weight(j)
            Tin   = (Tmin(i,j,k) + Tfin(i,j,k)) * (1.d0 - weight(j)) + (Tmout(i,j,k) + Tfout(i,j,k)) * weight(j)
            rhoin = pin / (R * Tin)
            Qre(i,j,k,1) = rhoin / Jacobian(i,j,k)
            Qre(i,j,k,2) = rhoin * uin / Jacobian(i,j,k)
            Qre(i,j,k,3) = rhoin * vin / Jacobian(i,j,k)
            Qre(i,j,k,4) = rhoin * win / Jacobian(i,j,k)
            Qre(i,j,k,5) = (pin / (gamma - 1.d0) + 0.5d0 * rhoin * (uin**2 + vin**2 + win**2)) / Jacobian(i,j,k)
      enddo;enddo;enddo
    else
      ! cyclic boundary condition !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (myrank == 3) then
        print *, "cyclic ", " blt=", real(bltre)
      endif
      do l = 1, 5
        do k = 1, nz
          do j = 1, ny
            ! inlet
            Qre(1,j,k,l) = Qre(1,j,k,l) / Jacobian(nre+1,j,k)
            Qre(2,j,k,l) = Qre(2,j,k,l) / Jacobian(nre+2,j,k)
      enddo;enddo;enddo
    endif
  end subroutine set_rescale
end module calc_rescale

