module calc_visc
  use mod_globals, only : accuracy, offset, id_visc, id_turbulence, id_dim, gamma, R, Re, Pr, nx, ny, nz, dxi, dyi, dzi
  use calc_Sutherland
  implicit none
contains
  attributes(device) function u_x(d,mu,a1,a2) result(ans)
    real(8), intent(in), value :: d, mu, a1, a2  
    real(8) :: ans
    ans = d * mu * (-a1 + a2)
  end function u_x

  attributes(device) function u_y(d,mu1,mu2,u1,u2,u3,u4,u5,u6) result(ans)
    real(8), intent(in), value :: d, mu1, mu2, u1, u2, u3, u4, u5, u6
    real(8) :: ans
    !y
    !!!!!!!!!!!!!
    ! u6     u5 !
    !    mu2    !
    ! u1     u4 !
    !    mu1    !
    ! u2     u3 !
    !!!!!!!!!!!!!x
    ans = d * 0.25d0 * (mu1 * (-u2 + u1 - u3 + u4) &
            & + mu2 * (-u1 + u6 - u4 + u5)) 
  end function u_y

  attributes(global) subroutine calc_Ev(xix, Jacobian, u, v, w, T, mut, E)
    real(8), intent(in), dimension(nx,ny,nz), device :: xix, Jacobian, u, v, w, T, mut
    real(8), intent(inout), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    integer i, j, k
    real(8) :: mux = 0.d0
    real(8) :: muy1 = 0.d0
    real(8) :: muy2 = 0.d0
    real(8) :: muz1 = 0.d0
    real(8) :: muz2 = 0.d0
    real(8) :: kappa = 0.d0
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    real(8) ux, uy, uz, vx, vy, vz, wx, wy, wz, txx, txy, txz, kappal, kappat, xixJ
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset - 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset

    ! calc viscosity
    if (id_visc == 1) then
      call calc_mu(T(i,j,k),T(i+1,j,k),mux)
      call calc_mu(T(i,j,k),T(i,j-1,k),T(i+1,j,k),T(i+1,j-1,k),muy1)
      call calc_mu(T(i,j,k),T(i,j+1,k),T(i+1,j,k),T(i+1,j+1,k),muy2)
      call calc_mu(T(i,j,k),T(i,j,k-1),T(i+1,j,k),T(i+1,j,k-1),muz1)
      call calc_mu(T(i,j,k),T(i,j,k+1),T(i+1,j,k),T(i+1,j,k+1),muz2)
      call calc_kappa(T(i,j,k),T(i+1,j,k),kappal)
    endif
    kappat = 0.5d0 * (mut(i,j,k) + mut(i+1,j,k)) * Cp / Pr
    kappa = kappal + kappat

    ! x direction
    mux = mux + 0.5d0 * (mut(i,j,k) + mut(i+1,j,k))

    ! y direction
    muy1 = muy1 + 0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i+1,j-1,k) + mut(i+1,j,k)) 
    muy2 = muy2 + 0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i+1,j,k) + mut(i+1,j+1,k))

    ! z direction
    muz1 = muz1 + 0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i+1,j,k-1) + mut(i+1,j,k))
    muz2 = muz2 + 0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i+1,j,k) + mut(i+1,j,k+1))

    ux = u_x(dxi,mux,u(i,j,k),u(i+1,j,k))
    vx = u_x(dxi,mux,v(i,j,k),v(i+1,j,k))
    wx = u_x(dxi,mux,w(i,j,k),w(i+1,j,k))
    uy = u_y(dyi,muy1,muy2,u(i,j,k),u(i,j-1,k),u(i+1,j-1,k),u(i+1,j,k),u(i+1,j+1,k),u(i,j+1,k))
    vy = u_y(dyi,muy1,muy2,v(i,j,k),v(i,j-1,k),v(i+1,j-1,k),v(i+1,j,k),v(i+1,j+1,k),v(i,j+1,k))
    uz = u_y(dzi,muz1,muz2,u(i,j,k),u(i,j,k-1),u(i+1,j,k-1),u(i+1,j,k),u(i+1,j,k+1),u(i,j,k+1))
    wz = u_y(dzi,muz1,muz2,w(i,j,k),w(i,j,k-1),w(i+1,j,k-1),w(i+1,j,k),w(i+1,j,k+1),w(i,j,k+1))

    xixJ = (xix(i,j,k) + xix(i+1,j,k)) / (Jacobian(i,j,k) + Jacobian(i+1,j,k))
    txx = (2.d0 * (2.d0 * ux - vy - wz) / 3.d0) * xixJ
    txy = (uy + vx) * xixJ
    txz = (wx + uz) * xixJ

    if (id_dim == 1) then
      E(i-offset+1,j-offset,k-offset,2) = E(i-offset+1,j-offset,k-offset,2) - txx
      E(i-offset+1,j-offset,k-offset,3) = E(i-offset+1,j-offset,k-offset,3) - txy
      E(i-offset+1,j-offset,k-offset,4) = E(i-offset+1,j-offset,k-offset,4) - txz
      E(i-offset+1,j-offset,k-offset,5) = E(i-offset+1,j-offset,k-offset,5) & 
      & - txx * 0.5d0 * (u(i,j,k) + u(i+1,j,k)) &
      & - txy * 0.5d0 * (v(i,j,k) + v(i+1,j,k)) &
      & - txz * 0.5d0 * (w(i,j,k) + w(i+1,j,k)) - kappa * (-T(i,j,k) + T(i+1,j,k)) * xixJ
    else
      E(i-offset+1,j-offset,k-offset,2) = E(i-offset+1,j-offset,k-offset,2) - txx / Re
      E(i-offset+1,j-offset,k-offset,3) = E(i-offset+1,j-offset,k-offset,3) - txy / Re
      E(i-offset+1,j-offset,k-offset,4) = E(i-offset+1,j-offset,k-offset,4) - txz / Re
      E(i-offset+1,j-offset,k-offset,5) = E(i-offset+1,j-offset,k-offset,5) &
      & -(txx * 0.5d0 * (u(i,j,k) + u(i+1,j,k)) &
      & + txy * 0.5d0 * (v(i,j,k) + v(i+1,j,k)) &
      & + txz * 0.5d0 * (w(i,j,k) + w(i+1,j,k)) + (-T(i,j,k) + T(i+1,j,k)) * xixJ / Pr) / Re
    endif
  end subroutine calc_Ev
  
  attributes(global) subroutine calc_Fv(etay, Jacobian, u, v, w, T, mut, F)
    real(8), intent(in), dimension(nx,ny,nz), device :: etay, Jacobian, u, v, w, T, mut
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    integer i, j, k
    real(8) :: muy = 0.d0
    real(8) :: muz1 = 0.d0
    real(8) :: muz2 = 0.d0
    real(8) :: mux1 = 0.d0
    real(8) :: mux2 = 0.d0
    real(8) :: kappa = 0.d0
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    real(8) ux, uy, uz, vx, vy, vz, wx, wy, wz, tyx, tyy, tyz, kappal, kappat, etayJ
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset - 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset

    ! calc viscosity
    if (id_visc == 1) then
      call calc_mu(T(i,j,k),T(i,j+1,k),muy)
      call calc_mu(T(i,j,k),T(i,j,k-1),T(i,j+1,k),T(i,j+1,k-1),muz1)
      call calc_mu(T(i,j,k),T(i,j,k+1),T(i,j+1,k),T(i,j+1,k+1),muz2)
      call calc_mu(T(i,j,k),T(i-1,j,k),T(i,j+1,k),T(i-1,j+1,k),mux1)
      call calc_mu(T(i,j,k),T(i+1,j,k),T(i,j+1,k),T(i+1,j+1,k),mux2)
      call calc_kappa(T(i,j,k),T(i,j+1,k),kappal)
    endif
    kappat = 0.5d0 * (mut(i,j,k) + mut(i,j+1,k)) * Cp / Pr
    kappa = kappal + kappat

    ! y direction
    muy = muy + 0.5d0 * (mut(i,j,k) + mut(i,j+1,k))

    ! z direction
    muz1 = muz1 + 0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i,j+1,k-1) + mut(i,j+1,k))
    muz2 = muz2 + 0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i,j+1,k) + mut(i,j+1,k+1))

    ! x direction
    mux1 = mux1 + 0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j+1,k) + mut(i,j+1,k))
    mux2 = mux2 + 0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j+1,k) + mut(i+1,j+1,k))

    vy = u_x(dyi,muy,v(i,j,k),v(i,j+1,k))
    wy = u_x(dyi,muy,w(i,j,k),w(i,j+1,k))
    uy = u_x(dyi,muy,u(i,j,k),u(i,j+1,k))
    vz = u_y(dzi,muz1,muz2,v(i,j,k),v(i,j,k-1),v(i,j+1,k-1),v(i,j+1,k),v(i,j+1,k+1),v(i,j,k+1)) 
    wz = u_y(dzi,muz1,muz2,w(i,j,k),w(i,j,k-1),w(i,j+1,k-1),w(i,j+1,k),w(i,j+1,k+1),w(i,j,k+1)) 
    ux = u_y(dxi,mux1,mux2,u(i,j,k),u(i-1,j,k),u(i-1,j+1,k),u(i,j+1,k),u(i+1,j+1,k),u(i+1,j,k)) 
    vx = u_y(dxi,mux1,mux2,v(i,j,k),v(i-1,j,k),v(i-1,j+1,k),v(i,j+1,k),v(i+1,j+1,k),v(i+1,j,k)) 

    etayJ = (etay(i,j,k) + etay(i,j+1,k)) / (Jacobian(i,j,k) + Jacobian(i,j+1,k))
    tyx = (uy + vx) * etayJ
    tyy = (2.d0 * (2.d0 * vy - wz - ux) / 3.d0) * etayJ
    tyz = (vz + wy) * etayJ
    if (id_dim == 1) then
      F(i-offset,j-offset+1,k-offset,2) = F(i-offset,j-offset+1,k-offset,2) - tyx
      F(i-offset,j-offset+1,k-offset,3) = F(i-offset,j-offset+1,k-offset,3) - tyy
      F(i-offset,j-offset+1,k-offset,4) = F(i-offset,j-offset+1,k-offset,4) - tyz
      F(i-offset,j-offset+1,k-offset,5) = F(i-offset,j-offset+1,k-offset,5) &
      & - tyx * 0.5d0 * (u(i,j,k) + u(i,j+1,k)) &
      & - tyy * 0.5d0 * (v(i,j,k) + v(i,j+1,k)) &
      & - tyz * 0.5d0 * (w(i,j,k) + w(i,j+1,k)) - kappa * (-T(i,j,k) + T(i,j+1,k)) * etayJ
    else
      F(i-offset,j-offset+1,k-offset,2) = F(i-offset,j-offset+1,k-offset,2) - tyx / Re
      F(i-offset,j-offset+1,k-offset,3) = F(i-offset,j-offset+1,k-offset,3) - tyy / Re
      F(i-offset,j-offset+1,k-offset,4) = F(i-offset,j-offset+1,k-offset,4) - tyz / Re
      F(i-offset,j-offset+1,k-offset,5) = F(i-offset,j-offset+1,k-offset,5) &
      & -(tyx * 0.5d0 * (u(i,j,k) + u(i,j+1,k)) &
      & + tyy * 0.5d0 * (v(i,j,k) + v(i,j+1,k)) &
      & + tyz * 0.5d0 * (w(i,j,k) + w(i,j+1,k)) + (-T(i,j,k) + T(i,j+1,k)) * etayJ / Pr) / Re
    endif
  end subroutine calc_Fv
  
  attributes(global) subroutine calc_Gv(Jacobian, u, v, w, T, mut, G)
    real(8), intent(in), dimension(nx,ny,nz), device :: Jacobian, u, v, w, T, mut
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    integer i, j, k
    real(8) :: muz = 0.d0
    real(8) :: mux1 = 0.d0
    real(8) :: mux2 = 0.d0
    real(8) :: muy1 = 0.d0
    real(8) :: muy2 = 0.d0
    real(8) :: kappa = 0.d0
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    real(8) ux, uy, uz, vx, vy, vz, wx, wy, wz, tzx, tzy, tzz, kappal, kappat, zetaJ
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset - 1
    
    ! calc viscosity
    if (id_visc == 1) then
      call calc_mu(T(i,j,k),T(i,j,k+1),muz)
      call calc_mu(T(i,j,k),T(i-1,j,k),T(i,j,k+1),T(i-1,j,k+1),mux1)
      call calc_mu(T(i,j,k),T(i+1,j,k),T(i,j,k+1),T(i+1,j,k+1),mux2)
      call calc_mu(T(i,j,k),T(i,j-1,k),T(i,j,k+1),T(i,j-1,k+1),muy1)
      call calc_mu(T(i,j,k),T(i,j+1,k),T(i,j,k+1),T(i,j+1,k+1),muy2)
      call calc_kappa(T(i,j,k),T(i,j,k+1),kappal)
    endif
    kappat = 0.5d0 * (mut(i,j,k) + mut(i,j,k+1)) * Cp / Pr
    kappa = kappal + kappat

    ! z direction
    muz = muz + 0.5d0 * (mut(i,j,k) + mut(i,j,k+1))
    
    ! x direction
    mux1 = mux1 + 0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j,k+1) + mut(i,j,k+1)) 
    mux2 = mux2 + 0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j,k+1) + mut(i+1,j,k+1))
    
    ! y direction
    muy1 = muy1 + 0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i,j-1,k+1) + mut(i,j,k+1))
    muy2 = muy2 + 0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i,j,k+1) + mut(i,j+1,k+1))

    wz = u_x(dzi,muz,w(i,j,k),w(i,j,k+1)) 
    uz = u_x(dzi,muz,u(i,j,k),u(i,j,k+1)) 
    vz = u_x(dzi,muz,v(i,j,k),v(i,j,k+1)) 
    wx = u_y(dxi,mux1,mux2,w(i,j,k),w(i-1,j,k),w(i-1,j,k+1),w(i,j,k+1),w(i+1,j,k+1),w(i+1,j,k))
    ux = u_y(dxi,mux1,mux2,u(i,j,k),u(i-1,j,k),u(i-1,j,k+1),u(i,j,k+1),u(i+1,j,k+1),u(i+1,j,k))
    vy = u_y(dyi,muy1,muy2,v(i,j,k),v(i,j-1,k),v(i,j-1,k+1),v(i,j,k+1),v(i,j+1,k+1),v(i,j+1,k)) 
    wy = u_y(dyi,muy1,muy2,w(i,j,k),w(i,j-1,k),w(i,j-1,k+1),w(i,j,k+1),w(i,j+1,k+1),w(i,j+1,k))  

    zetaJ = 2.d0 / (Jacobian(i,j,k) + Jacobian(i,j,k+1))
    tzx = (wx + uz) * zetaJ
    tzy = (vz + wy) * zetaJ
    tzz = (2.d0 * (2.d0 * wz - ux - vy) / 3.d0) * zetaJ 
    if (id_dim == 1) then
      G(i-offset,j-offset,k-offset+1,2) = G(i-offset,j-offset,k-offset+1,2) - tzx
      G(i-offset,j-offset,k-offset+1,3) = G(i-offset,j-offset,k-offset+1,3) - tzy
      G(i-offset,j-offset,k-offset+1,4) = G(i-offset,j-offset,k-offset+1,4) - tzz
      G(i-offset,j-offset,k-offset+1,5) = G(i-offset,j-offset,k-offset+1,5) &
      & - tzx * 0.5d0 * (u(i,j,k) + u(i,j,k+1)) &
      & - tzy * 0.5d0 * (v(i,j,k) + v(i,j,k+1)) &
      & - tzz * 0.5d0 * (w(i,j,k) + w(i,j,k+1)) - kappa * (-T(i,j,k) + T(i,j,k+1)) * zetaJ
    else
      G(i-offset,j-offset,k-offset+1,2) = G(i-offset,j-offset,k-offset+1,2) - tzx / Re
      G(i-offset,j-offset,k-offset+1,3) = G(i-offset,j-offset,k-offset+1,3) - tzy / Re
      G(i-offset,j-offset,k-offset+1,4) = G(i-offset,j-offset,k-offset+1,4) - tzz / Re
      G(i-offset,j-offset,k-offset+1,5) = G(i-offset,j-offset,k-offset+1,5) &
      & -(tzx * 0.5d0 * (u(i,j,k) + u(i,j,k+1)) &
      & + tzy * 0.5d0 * (v(i,j,k) + v(i,j,k+1)) &
      & + tzz * 0.5d0 * (w(i,j,k) + w(i,j,k+1)) + (-T(i,j,k) + T(i,j,k+1)) * zetaJ / Pr) / Re
    endif
  end subroutine calc_Gv
end module calc_visc

