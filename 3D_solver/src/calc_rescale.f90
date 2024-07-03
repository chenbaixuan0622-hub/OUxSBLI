module calc_rescale
  use openacc
  use mpi
  use mod_globals, only : gamma , R, u0
contains
  subroutine calc_mean(step,nx,ny,nz,nre,QJ,Jacobians,ure,vre,wre,rhore,Tre)
    integer, intent(in), value                          :: step, nx, ny, nz, nre ! nre: rescale plane
    real(8), intent(in), dimension(nx,ny,nz,5), device  :: QJ
    real(8), intent(in), dimension(nx,ny), device       :: Jacobians
    real(8), intent(inout), dimension(ny), device       :: ure, vre, wre, rhore, Tre
    integer j, k
    real(8) pJ
    real(8), dimension(ny), device :: us, vs, ws, rhos, Ts
    ! get mean value for span
    us(:) = 0.d0; vs(:) = 0.d0; ws(:) = 0.d0; rhos(:) = 0.d0; Ts(:) = 0.d0
    !$acc parallel loop reduction(+:us, vs, ws, rhos, Ts)
    do k = 2, nz-1
      do j = 1, ny
        us(j)   = us(j) +   QJ(nre,j,k,2) / QJ(nre,j,k,1) 
        vs(j)   = vs(j) +   QJ(nre,j,k,3) / QJ(nre,j,k,1)
        ws(j)   = ws(j) +   QJ(nre,j,k,4) / QJ(nre,j,k,1)
        rhos(j) = rhos(j) + QJ(nre,j,k,1) * Jacobian(nre,j)
        pJ = (gamma - 1.d0) * (QJ(nre,j,k,5) - 0.5d0 * (QJ(nre,j,k,2)**2 + QJ(nre,j,k,3)**2 + QJ(nre,j,k,4)**2) / QJ(nre,j,k,1))
        Ts(j)   = Ts(j) + pJ / (R * QJ(nre,j,k,1))
    enddo;enddo

    ! get mean value for time
    !$cuf kernel do(1) <<<*,*>>>
    do j = 1, ny
      ure(j) =   ((dble(step) - 1.d0) * ure(j)   + us(j)   / dble(nz-2)) / dble(step)
      vre(j) =   ((dble(step) - 1.d0) * vre(j)   + vs(j)   / dble(nz-2)) / dble(step)
      wre(j) =   ((dble(step) - 1.d0) * wre(j)   + ws(j)   / dble(nz-2)) / dble(step)
      rhore(j) = ((dble(step) - 1.d0) * rhore(j) + rhos(j) / dble(nz-2)) / dble(step)
      Tre(j) =   ((dble(step) - 1.d0) * Tre(j)   + Ts(j)   / dble(nz-2)) / dble(step)
    enddo
  end subroutine calc_mean

  subroutine set_rescale(nx,ny,nz,y,Um,Vm,Wm,pm,Tm,QJ)
    integer, intent(in), value                   :: nx, ny, nz
    real(8), intent(in), device                  :: y(ny)
    real(8), intent(in), dimension(2,ny), device :: Um, Vm, Wm, pm, Tm ! mean properties at rescaling plane
    real(8), intent(inout), device               :: QJ(nx,ny,nz,5)
    integer i, j, jj, k, l
    real(8) blt99
    real(8) :: blt = 2.d-3
    ! fluctuating properties at rescaling plane
    real(8), dimension(nx,ny,nz), device :: ufre, vfre, wfre, pfre, Tfre

    ! calc boundary layer thickness at rescaling plane
    do j = 2, ny
      if (u(j) >= 0.99d0 * u0) then
        blt99 = y(j) - (-y(j-1) + y(j)) * (Um(j) - 0.99d0 * u0) / (-Um(j-1) + Um(j))
        exit
      endif
    enddo

    if (blt99 <= blt) then
      ! cyclic boundary condition !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do l = 1, 5
        do k = 2, nz-1
          do j = 1, ny
            ! inlet
            QJ(1,j,k,l) = QJ(nx,j,k,l)
            QJ(2,j,k,l) = QJ(nx,j,k,l)
      enddo;enddo;enddo;enddo
    else
      ! rescaling !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! calc fluctuating part   u'(x,y,z,t) = u(x,y,z,t) - U(x,y)
      ! U(x,y) average velocity in the spanwise direction and time
      do k = 2, nz-1
        do j = 2, ny-1
          do i = 1, 2 ! 2 rescaleing planes are required for 4th-order accuracy flux
            ufre(i,j,k) = u(i+offsetre,j,k) - Um(i,j)
            vfre(i,j,k) = v(i+offsetre,j,k) - Vm(i,j)
            wfre(i,j,k) = w(i+offsetre,j,k) - Wm(i,j)
            pfre(i,j,k) = p(i+offsetre,j,k) - pm(i,j)
            Tfre(i,j,k) = T(i+offsetre,j,k) - Tm(i,j)
      enddo;enddo;enddo

      ! friction velocity
      rhore = pr(1) / (R * Tr(1))
      ! linear interpolation
      taure = mu(Tr(1)) * abs(ur(2)) / (-y(1) + y(2))
      utre  = sqrt(taur / rhor)
      utin  = utre * (bltre / blt)**0.1
      ! rescaling factor
      beta = utin / utre

      ! inner region
      ! y plus
      do j = 1, ny
        ypin(j) = 
        ypre(j) = 
      enddo
      do j = 2, ny
        do jj = 2, ny
          if (ypre(jj) > ypin(j)) then
            ady = (-ypre(jj-1) + ypin(j)) / (-ypre(jj-1) + ypre(jj))
            do k = 2, nz-1
              do i = 1, 2
                ! mean
                Umin(i,j,k) = beta * (Um(i,jj-1,k) + ady * (-Um(i,jj-1,k) + Um(i,jj,k)))
                Vmin(i,j,k) =         Vm(i,jj-1,k) + ady * (-Vm(i,jj-1,k) + Vm(i,jj,k))
                Wmin(i,j,k) = 0.d0
                pmin(i,j,k) =         pm(i,jj-1,k) + ady * (-pm(i,jj-1,k) + pm(i,jj,k)) 
                Tmin(i,j,k) =         Tm(i,jj-1,k) + ady * (-Tm(i,jj-1,k) + Tm(i,jj,k))
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

      ! outer region
      ! eta
      do j = 1, ny
        etin(j) =
        etre(j) = 
      enddo
      do j = 2, ny
        do jj = 2, ny
          if (ypre(jj) > ypin(j)) then
            ade = (-etre(jj-1) + etin(j)) / (-etre(jj-1) + etre(jj))
            do k = 2, nz-1
              do i = 1, 2
                ! mean
                Umout(i,j,k) = beta * (Um(i,jj-1,k) + ade * (-Um(i,jj-1,k) + Um(i,jj,k))) + (1.d0 - beta) * u0
                Vmout(i,j,k) =         Vm(i,jj-1,k) + ade
                Wmout(i,j,k) = 0.d0
                pmout(i,j,k) =         pm(i,jj-1,k) + ade * (-pm(i,jj-1,k) + pm(i,jj,k))
                Tmout(i,j,k) =         Tm(i,jj-1,k) + ade * (-Tm(i,jj-1,k) + Tm(i,jj,k))
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
      do j = 2, ny
        weight(j) = min(1.d0, &
                   0.5d0 * (1.d0 + tanh(4.d0 * (etin(j) - 0.2d0) / ((1.d0 - 0.4d0) * etin(j) + 0.2d0)) / tanh(4.d0)))
      enddo

      ! re-introducing
      do k = 2, nz-1
        do j = 1, ny
          do i = 1, 2
            uin(i,j,k) = (Umin(i,j,k) + ufin(i,j,k)) * (1.d0 - weight(j)) + (Umout(i,j,k) + ufout(i,j,k)) * weight(j)
            vin(i,j,k) = (Vmin(i,j,k) + vfin(i,j,k)) * (1.d0 - weight(j)) + (Vmout(i,j,k) + vfout(i,j,k)) * weight(j)
            win(i,j,k) = (Wmin(i,j,k) + wfin(i,j,k)) * (1.d0 - weight(j)) + (Wmout(i,j,k) + wfout(i,j,k)) * weight(j)
            pin(i,j,k) = (pmin(i,j,k) + pfin(i,j,k)) * (1.d0 - weight(j)) + (pmout(i,j,k) + pfout(i,j,k)) * weight(j)
            Tin(i,j,k) = (Tmin(i,j,k) + Tfin(i,j,k)) * (1.d0 - weight(j)) + (Tmout(i,j,k) + Tfout(i,j,k)) * weight(j)
            rhoin(i,j,k) = pin(i,j,k) / (R * Tin(i,j,k))
      enddo;enddo;enddo
    endif
  end subroutine set_rescale
end module calc_rescale

