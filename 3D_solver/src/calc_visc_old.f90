module calc_visc
  use mod_globals, only : accuracy, offset, id_visc, id_turbulence, gamma, R, Pr, Prt
  use calc_Sutherland
  implicit none
  interface u_x
    module procedure u_x2, u_x4
  end interface
  
  interface u_y
    module procedure u_y2, u_y4
  end interface
contains
  attributes(device) function u_x2(d,mu,a1,a2) result(ans)
    real(8), intent(in), value :: d, mu, a1, a2  
    real(8) :: ans
    ans = d * mu * (-a1 + a2)
  end function u_x2

  attributes(device) function u_x4(d,mu,a1,a2,a3,a4) result(ans)
    real(8), intent(in), value :: d, mu, a1, a2, a3, a4
    real(8) :: ans
    ans = d * mu * (a1 - 27.d0 * a2 + 27.d0 * a3 - a4) / 24.d0
  end function u_x4

  attributes(device) function u_y2(d1,d2,mu1,mu2,u1,u2,u3,u4,u5,u6) result(ans)
    real(8), intent(in), value :: d1, d2, mu1, mu2, u1, u2, u3, u4, u5, u6
    real(8) :: ans
    !y
    !!!!!!!!!!!!!
    ! u6     u5 !
    !    mu2    !
    ! u1     u4 !
    !    mu1    !
    ! u2     u3 !
    !!!!!!!!!!!!!x
    ans = 0.25d0 * (d1 * mu1 * (-u2 + u1 - u3 + u4) &
                & + d2 * mu2 * (-u1 + u6 - u4 + u5)) 
  end function u_y2

  attributes(device) function u_y4(d1,d2,mu1,mu2,u1,u2,u3,u4,u5,u6,u7,u8,u9,u10) result(ans)
    real(8), intent(in), value :: d1, d2, mu1, mu2, u1, u2, u3, u4, u5, u6, u7, u8, u9, u10
    real(8) :: ans
    !y
    !!!!!!!!!!!!!
    ! u9     u8 !
    !           !
    ! u10    u7 !
    !    mu2    !
    ! u1     u6 !
    !    mu1    !
    ! u2     u5 !
    !           !
    ! u3     u4 !
    !!!!!!!!!!!!!x
    ans = (d1 * mu1 * ((u3 - 27.d0 * u2 + 27.d0 * u1 - u10) + &
                       (u4 - 27.d0 * u5 + 27.d0 * u6 - u7)) + &
           d2 * mu2 * ((u2 - 27.d0 * u1 + 27.d0 * u10- u9)  + &
                       (u5 - 27.d0 * u6 + 27.d0 * u7 - u8))) / 96.d0
  end function u_y4

  attributes(device) function u_m(a1,a2,a3,a4) result(ans)
    real(8), intent(in), value :: a1, a2, a3, a4
    real(8) :: ans
    ans = 0.06250 * (-a1 + 9.d0 * a2 + 9.d0 * a3 - a4)
  end function u_m

  attributes(global) subroutine calc_Ev(nx, ny, nz, dx, dy, dz, rho, u, v, w, T, p, mut, qc2, E)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device      :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1), device      :: dy ! 1 / dy
    real(8), intent(in), dimension(nz-1), device      :: dz ! 1 / dz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, T, p, mut, qc2
    real(8), intent(inout), device                    :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    integer i, j, k
    real(8) :: mux = 0.d0
    real(8) :: muy1 = 0.d0
    real(8) :: muy2 = 0.d0
    real(8) :: muz1 = 0.d0
    real(8) :: muz2 = 0.d0
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    real(8) um, vm, wm, ux, uy, uz, vx, vy, wx, wz, txx, txy, txz, Tx, H1, H2, kappa
    real(8) :: txxsgs = 0.d0, txysgs = 0.d0, txzsgs = 0.d0, Hsgs = 0.d0
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset - 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset

    ! calc viscosity
    if (1 <= id_visc) then
      call calc_mu(T(i,j,k),T(i+1,j,k),mux)
      call calc_mu(T(i,j,k),T(i,j-1,k),T(i+1,j,k),T(i+1,j-1,k),muy1)
      call calc_mu(T(i,j,k),T(i,j+1,k),T(i+1,j,k),T(i+1,j+1,k),muy2)
      call calc_mu(T(i,j,k),T(i,j,k-1),T(i+1,j,k),T(i+1,j,k-1),muz1)
      call calc_mu(T(i,j,k),T(i,j,k+1),T(i+1,j,k),T(i+1,j,k+1),muz2)
    endif
    kappa = Cp * mux / Pr

    if (2 <= i .and. i <= nx-2 .and. id_visc == 2) then
      ux = u_x(dx(i),mux,u(i-1,j,k),u(i,j,k),u(i+1,j,k),u(i+2,j,k))
      vx = u_x(dx(i),mux,v(i-1,j,k),v(i,j,k),v(i+1,j,k),v(i+2,j,k))
      wx = u_x(dx(i),mux,w(i-1,j,k),w(i,j,k),w(i+1,j,k),w(i+2,j,k))
      Tx = u_x(dx(i),kappa,T(i-1,j,k),T(i,j,k),T(i+1,j,k),T(i+2,j,k))
      um = u_m(u(i-1,j,k),u(i,j,k),u(i+1,j,k),u(i+2,j,k))
      vm = u_m(v(i-1,j,k),v(i,j,k),v(i+1,j,k),v(i+2,j,k))
      wm = u_m(w(i-1,j,k),w(i,j,k),w(i+1,j,k),w(i+2,j,k))
    else
      ux = u_x(dx(i),mux,u(i,j,k),u(i+1,j,k))
      vx = u_x(dx(i),mux,v(i,j,k),v(i+1,j,k))
      wx = u_x(dx(i),mux,w(i,j,k),w(i+1,j,k))
      Tx = u_x(dx(i),kappa,T(i,j,k),T(i+1,j,k))
      um = 0.5d0 * (u(i,j,k) + u(i+1,j,k))
      vm = 0.5d0 * (v(i,j,k) + v(i+1,j,k))
      wm = 0.5d0 * (w(i,j,k) + w(i+1,j,k))
    endif
    if (3 <= j .and. j <= ny-2 .and. id_visc == 2) then
      uy = u_y(dy(j-1),dy(j),muy1,muy2,u(i,j,k),u(i,j-1,k),u(i,j-2,k),u(i+1,j-2,k),u(i+1,j-1,k),u(i+1,j,k),&
                                       u(i+1,j+1,k),u(i+1,j+2,k),u(i,j+2,k),u(i,j+1,k))
      vy = u_y(dy(j-1),dy(j),muy1,muy2,v(i,j,k),v(i,j-1,k),v(i,j-2,k),v(i+1,j-2,k),v(i+1,j-1,k),v(i+1,j,k),&
                                       v(i+1,j+1,k),v(i+1,j+2,k),v(i,j+2,k),v(i,j+1,k))
    else
      uy = u_y(dy(j-1),dy(j),muy1,muy2,u(i,j,k),u(i,j-1,k),u(i+1,j-1,k),u(i+1,j,k),u(i+1,j+1,k),u(i,j+1,k))
      vy = u_y(dy(j-1),dy(j),muy1,muy2,v(i,j,k),v(i,j-1,k),v(i+1,j-1,k),v(i+1,j,k),v(i+1,j+1,k),v(i,j+1,k))
    endif
    if (3 <= k .and. k <= nz-2 .and. id_visc == 2) then
      uz = u_y(dz(k-1),dz(k),muz1,muz2,u(i,j,k),u(i,j,k-1),u(i,j,k-2),u(i+1,j,k-2),u(i+1,j,k-1),u(i+1,j,k),&
                                       u(i+1,j,k+1),u(i+1,j,k+2),u(i,j,k+2),u(i,j,k+1))
      wz = u_y(dz(k-1),dz(k),muz1,muz2,w(i,j,k),w(i,j,k-1),w(i,j,k-2),w(i+1,j,k-2),w(i+1,j,k-1),w(i+1,j,k),&
                                       w(i+1,j,k+1),w(i+1,j,k+2),w(i,j,k+2),w(i,j,k+1))
    else
      uz = u_y(dz(k-1),dz(k),muz1,muz2,u(i,j,k),u(i,j,k-1),u(i+1,j,k-1),u(i+1,j,k),u(i+1,j,k+1),u(i,j,k+1))
      wz = u_y(dz(k-1),dz(k),muz1,muz2,w(i,j,k),w(i,j,k-1),w(i+1,j,k-1),w(i+1,j,k),w(i+1,j,k+1),w(i,j,k+1))
    endif

    txx = 2.d0 * (2.d0 * ux - vy - wz) / 3.d0
    txy = uy + vx
    txz = wx + uz

    if (id_turbulence == 1) then
      ! x direction
      mux = 0.5d0 * (mut(i,j,k) + mut(i+1,j,k))
      ! y direction
      muy1 = 0.25d0 * (mut(i,j-1,k) + mut(i,j,k)   + mut(i+1,j-1,k) + mut(i+1,j,k)) 
      muy2 = 0.25d0 * (mut(i,j,k)   + mut(i,j+1,k) + mut(i+1,j,k)   + mut(i+1,j+1,k))
      ! z direction
      muz1 = 0.25d0 * (mut(i,j,k-1) + mut(i,j,k)   + mut(i+1,j,k-1) + mut(i+1,j,k))
      muz2 = 0.25d0 * (mut(i,j,k)   + mut(i,j,k+1) + mut(i+1,j,k)   + mut(i+1,j,k+1))
    
      ux = u_x(dx(i),mux,u(i,j,k),u(i+1,j,k))
      vx = u_x(dx(i),mux,v(i,j,k),v(i+1,j,k))
      wx = u_x(dx(i),mux,w(i,j,k),w(i+1,j,k))
      uy = u_y(dy(j-1),dy(j),muy1,muy2,u(i,j,k),u(i,j-1,k),u(i+1,j-1,k),u(i+1,j,k),u(i+1,j+1,k),u(i,j+1,k))
      vy = u_y(dy(j-1),dy(j),muy1,muy2,v(i,j,k),v(i,j-1,k),v(i+1,j-1,k),v(i+1,j,k),v(i+1,j+1,k),v(i,j+1,k))
      uz = u_y(dz(k-1),dz(k),muz1,muz2,u(i,j,k),u(i,j,k-1),u(i+1,j,k-1),u(i+1,j,k),u(i+1,j,k+1),u(i,j,k+1))
      wz = u_y(dz(k-1),dz(k),muz1,muz2,w(i,j,k),w(i,j,k-1),w(i+1,j,k-1),w(i+1,j,k),w(i+1,j,k+1),w(i,j,k+1))
      ! sgs visc
      txxsgs = 2.d0 * (2.d0 * ux - vy - wz - (rho(i,j,k) + rho(i+1,j,k)) * qc2(i,j,k)) / 3.d0
      txysgs = uy + vx
      txzsgs = wx + uz
      ! sgs kappa
      kappa = kappa + 0.5d0 * (mut(i,j,k) + mut(i+1,j,k)) * Cp / Prt
      ! sgs enthalpy
      H1 = (gamma * p(i,j,k)   / (rho(i,j,k)   * (gamma - 1.d0))) + 0.5d0 * (u(i,j,k)**2   + v(i,j,k)**2   + w(i,j,k)**2)   + qc2(i,j,k) 
      H2 = (gamma * p(i+1,j,k) / (rho(i+1,j,k) * (gamma - 1.d0))) + 0.5d0 * (u(i+1,j,k)**2 + v(i+1,j,k)**2 + w(i+1,j,k)**2) + qc2(i+1,j,k)
      Hsgs = -mux * (-H1 + H2) * dx(i) / Prt
    endif
    E(i-offset+1,j-offset,k-offset,2) = E(i-offset+1,j-offset,k-offset,2) - (txx+txxsgs)
    E(i-offset+1,j-offset,k-offset,3) = E(i-offset+1,j-offset,k-offset,3) - (txy+txysgs)
    E(i-offset+1,j-offset,k-offset,4) = E(i-offset+1,j-offset,k-offset,4) - (txz+txzsgs)
    E(i-offset+1,j-offset,k-offset,5) = E(i-offset+1,j-offset,k-offset,5) &
    & -(txx * um + txy * vm + txz * wm + Tx + Hsgs)
  end subroutine calc_Ev
  
  attributes(global) subroutine calc_Fv(nx, ny, nz, dy, dx, dz, rho, u, v, w, T, p, mut, qc2, F)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(ny-1), device      :: dy ! 1 / dy
    real(8), intent(in), dimension(nx-1), device      :: dx ! 1 / dx
    real(8), intent(in), dimension(nz-1), device      :: dz ! 1 / dz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, T, p, mut, qc2
    real(8), intent(inout), device                    :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    integer i, j, k
    real(8) :: muy = 0.d0
    real(8) :: muz1 = 0.d0
    real(8) :: muz2 = 0.d0
    real(8) :: mux1 = 0.d0
    real(8) :: mux2 = 0.d0
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    real(8) um, vm, wm, ux, uy, vx, vy, vz, wy, wz, tyx, tyy, tyz, Ty, H1, H2, kappa
    real(8) :: tyxsgs = 0.d0, tyysgs = 0.d0, tyzsgs = 0.d0, Hsgs = 0.d0
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset - 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset

    ! calc viscosity
    if (1 <= id_visc) then
      call calc_mu(T(i,j,k),T(i,j+1,k),muy)
      call calc_mu(T(i,j,k),T(i,j,k-1),T(i,j+1,k),T(i,j+1,k-1),muz1)
      call calc_mu(T(i,j,k),T(i,j,k+1),T(i,j+1,k),T(i,j+1,k+1),muz2)
      call calc_mu(T(i,j,k),T(i-1,j,k),T(i,j+1,k),T(i-1,j+1,k),mux1)
      call calc_mu(T(i,j,k),T(i+1,j,k),T(i,j+1,k),T(i+1,j+1,k),mux2)
    endif
    kappa = Cp * muy / Pr

    if (2 <= j .and. j <= ny-2 .and. id_visc == 2) then
      vy = u_x(dy(j),muy,v(i,j-1,k),v(i,j,k),v(i,j+1,k),v(i,j+2,k))
      wy = u_x(dy(j),muy,w(i,j-1,k),w(i,j,k),w(i,j+1,k),w(i,j+2,k))
      uy = u_x(dy(j),muy,u(i,j-1,k),u(i,j,k),u(i,j+1,k),u(i,j+2,k))
      Ty = u_x(dy(j),kappa,T(i,j-1,k),T(i,j,k),T(i,j+1,k),T(i,j+2,k))
      um = u_m(u(i,j-1,k),u(i,j,k),u(i,j+1,k),u(i,j+2,k))
      vm = u_m(v(i,j-1,k),v(i,j,k),v(i,j+1,k),v(i,j+2,k))
      wm = u_m(w(i,j-1,k),w(i,j,k),w(i,j+1,k),w(i,j+2,k))
    else
      vy = u_x(dy(j),muy,v(i,j,k),v(i,j+1,k))
      wy = u_x(dy(j),muy,w(i,j,k),w(i,j+1,k))
      uy = u_x(dy(j),muy,u(i,j,k),u(i,j+1,k))
      Ty = u_x(dy(j),kappa,T(i,j,k),T(i,j+1,k))
      um = 0.5d0 * (u(i,j,k) + u(i,j+1,k))
      vm = 0.5d0 * (v(i,j,k) + v(i,j+1,k))
      wm = 0.5d0 * (w(i,j,k) + w(i,j+1,k))
    endif
    if (3 <= k .and. k <= nz-2 .and. id_visc == 2) then
      vz = u_y(dz(k-1),dz(k),muz1,muz2,v(i,j,k),v(i,j,k-1),v(i,j,k-2),v(i,j+1,k-2),v(i,j+1,k-1),v(i,j+1,k),&
                                       v(i,j+1,k+1),v(i,j+1,k+2),v(i,j,k+2),v(i,j,k+1))
      wz = u_y(dz(k-1),dz(k),muz1,muz2,w(i,j,k),w(i,j,k-1),w(i,j,k-2),w(i,j+1,k-2),w(i,j+1,k-1),w(i,j+1,k),&
                                       w(i,j+1,k+1),w(i,j+1,k+2),w(i,j,k+2),w(i,j,k+1))
    else
      vz = u_y(dz(k-1),dz(k),muz1,muz2,v(i,j,k),v(i,j,k-1),v(i,j+1,k-1),v(i,j+1,k),v(i,j+1,k+1),v(i,j,k+1)) 
      wz = u_y(dz(k-1),dz(k),muz1,muz2,w(i,j,k),w(i,j,k-1),w(i,j+1,k-1),w(i,j+1,k),w(i,j+1,k+1),w(i,j,k+1)) 
    endif
    if (3 <= i .and. i <= nx-2 .and. id_visc == 2) then
      ux = u_y(dx(i-1),dx(i),mux1,mux2,u(i,j,k),u(i-1,j,k),u(i-2,j,k),u(i-2,j+1,k),u(i-1,j+1,k),u(i,j+1,k),&
                                       u(i+1,j+1,k),u(i+2,j+1,k),u(i+2,j,k),u(i+1,j,k))
      vx = u_y(dx(i-1),dx(i),mux1,mux2,v(i,j,k),v(i-1,j,k),v(i-2,j,k),v(i-2,j+1,k),v(i-1,j+1,k),v(i,j+1,k),&
                                       v(i+1,j+1,k),v(i+2,j+1,k),v(i+2,j,k),u(i+1,j,k))
    else
      ux = u_y(dx(i-1),dx(i),mux1,mux2,u(i,j,k),u(i-1,j,k),u(i-1,j+1,k),u(i,j+1,k),u(i+1,j+1,k),u(i+1,j,k)) 
      vx = u_y(dx(i-1),dx(i),mux1,mux2,v(i,j,k),v(i-1,j,k),v(i-1,j+1,k),v(i,j+1,k),v(i+1,j+1,k),v(i+1,j,k)) 
    endif

    tyx = uy + vx
    tyy = 2.d0 * (2.d0 * vy - wz - ux) / 3.d0
    tyz = vz + wy

    if (id_turbulence == 1) then
      ! y direction
      muy = 0.5d0 * (mut(i,j,k) + mut(i,j+1,k))
      ! z direction
      muz1 = 0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i,j+1,k-1) + mut(i,j+1,k))
      muz2 = 0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i,j+1,k) + mut(i,j+1,k+1))
      ! x direction
      mux1 = 0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j+1,k) + mut(i,j+1,k))
      mux2 = 0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j+1,k) + mut(i+1,j+1,k))

      vy = u_x(dy(j),muy,v(i,j,k),v(i,j+1,k))
      wy = u_x(dy(j),muy,w(i,j,k),w(i,j+1,k))
      uy = u_x(dy(j),muy,u(i,j,k),u(i,j+1,k))
      vz = u_y(dz(k-1),dz(k),muz1,muz2,v(i,j,k),v(i,j,k-1),v(i,j+1,k-1),v(i,j+1,k),v(i,j+1,k+1),v(i,j,k+1)) 
      wz = u_y(dz(k-1),dz(k),muz1,muz2,w(i,j,k),w(i,j,k-1),w(i,j+1,k-1),w(i,j+1,k),w(i,j+1,k+1),w(i,j,k+1)) 
      ux = u_y(dx(i-1),dx(i),mux1,mux2,u(i,j,k),u(i-1,j,k),u(i-1,j+1,k),u(i,j+1,k),u(i+1,j+1,k),u(i+1,j,k)) 
      vx = u_y(dx(i-1),dx(i),mux1,mux2,v(i,j,k),v(i-1,j,k),v(i-1,j+1,k),v(i,j+1,k),v(i+1,j+1,k),v(i+1,j,k)) 
      ! sgs visc
      tyxsgs = uy + vx
      tyysgs = 2.d0 * (2.d0 * vy - wz - ux - (rho(i,j,k) + rho(i,j+1,k)) * qc2(i,j,k)) / 3.d0
      tyzsgs = vz + wy
      ! sgs kappa
      kappa = kappa + 0.5d0 * (mut(i,j,k) + mut(i,j+1,k)) * Cp / Prt
      ! sgs enthalpy
      H1 = (gamma * p(i,j,k)   / (rho(i,j,k)   * (gamma - 1.d0))) + 0.5d0 * (u(i,j,k)**2   + v(i,j,k)**2   + w(i,j,k)**2)   + qc2(i,j,k) 
      H2 = (gamma * p(i,j+1,k) / (rho(i,j+1,k) * (gamma - 1.d0))) + 0.5d0 * (u(i,j+1,k)**2 + v(i,j+1,k)**2 + w(i,j+1,k)**2) + qc2(i,j+1,k)
      Hsgs = -muy * (-H1 + H2) * dy(j) / Prt
    endif
    F(i-offset,j-offset+1,k-offset,2) = F(i-offset,j-offset+1,k-offset,2) - (tyx+tyxsgs)
    F(i-offset,j-offset+1,k-offset,3) = F(i-offset,j-offset+1,k-offset,3) - (tyy+tyysgs)
    F(i-offset,j-offset+1,k-offset,4) = F(i-offset,j-offset+1,k-offset,4) - (tyz+tyzsgs)
    F(i-offset,j-offset+1,k-offset,5) = F(i-offset,j-offset+1,k-offset,5) &
    & -(tyx * um + tyy * vm + tyz * wm + Ty + Hsgs)
  end subroutine calc_Fv
  
  attributes(global) subroutine calc_Gv(nx, ny, nz, dx, dy, dz, rho, u, v, w, T, p, mut, qc2, G)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device      :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1), device      :: dy ! 1 / dy
    real(8), intent(in), dimension(nz-1), device      :: dz ! 1 / dz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, T, p, mut, qc2
    real(8), intent(inout), device                    :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    integer i, j, k
    real(8) :: muz = 0.d0
    real(8) :: mux1 = 0.d0
    real(8) :: mux2 = 0.d0
    real(8) :: muy1 = 0.d0
    real(8) :: muy2 = 0.d0
    real(8) :: Cp = gamma * R / (gamma - 1.d0)
    real(8) um, vm, wm, ux, uz, vy, vz, wx, wy, wz, tzx, tzy, tzz, Tz, H1, H2, kappa
    real(8) :: tzxsgs = 0.d0, tzysgs = 0.d0, tzzsgs = 0.d0, Hsgs = 0.d0
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset - 1
    
    ! calc viscosity
    if (1 <= id_visc) then
      call calc_mu(T(i,j,k),T(i,j,k+1),muz)
      call calc_mu(T(i,j,k),T(i-1,j,k),T(i,j,k+1),T(i-1,j,k+1),mux1)
      call calc_mu(T(i,j,k),T(i+1,j,k),T(i,j,k+1),T(i+1,j,k+1),mux2)
      call calc_mu(T(i,j,k),T(i,j-1,k),T(i,j,k+1),T(i,j-1,k+1),muy1)
      call calc_mu(T(i,j,k),T(i,j+1,k),T(i,j,k+1),T(i,j+1,k+1),muy2)
    endif
    kappa = Cp * muz / Pr

    if (2 <= k .and. k <= nz-2 .and. id_visc == 2) then
      wz = u_x(dz(k),muz,w(i,j,k-1),w(i,j,k),w(i,j,k+1),w(i,j,k+2)) 
      uz = u_x(dz(k),muz,u(i,j,k-1),u(i,j,k),u(i,j,k+1),u(i,j,k+2))
      vz = u_x(dz(k),muz,v(i,j,k-1),v(i,j,k),v(i,j,k+1),v(i,j,k+2))
      Tz = u_x(dz(k),kappa,T(i,j,k-1),T(i,j,k),T(i,j,k+1),T(i,j,k+2))
      um = u_m(u(i,j,k-1),u(i,j,k),u(i,j,k+1),u(i,j,k+2))
      vm = u_m(v(i,j,k-1),v(i,j,k),v(i,j,k+1),v(i,j,k+2))
      wm = u_m(w(i,j,k-1),w(i,j,k),w(i,j,k+1),w(i,j,k+2))
    else
      wz = u_x(dz(k),muz,w(i,j,k),w(i,j,k+1)) 
      uz = u_x(dz(k),muz,u(i,j,k),u(i,j,k+1)) 
      vz = u_x(dz(k),muz,v(i,j,k),v(i,j,k+1)) 
      Tz = u_x(dz(k),kappa,T(i,j,k),T(i,j,k+1))
      um = 0.5d0 * (u(i,j,k)+ u(i,j,k+1))
      vm = 0.5d0 * (v(i,j,k)+ v(i,j,k+1))
      wm = 0.5d0 * (w(i,j,k)+ w(i,j,k+1))
    endif
    if (3 <= i .and. i <= nx-2 .and. id_visc == 2) then
      wx = u_y(dx(i-1),dx(i),mux1,mux2,w(i,j,k),w(i-1,j,k),w(i-2,j,k),w(i-2,j,k+1),w(i-1,j,k+1),w(i,j,k+1),&
                                       w(i+1,j,k+1),w(i+2,j,k+1),w(i+2,j,k),w(i+1,j,k))
      ux = u_y(dx(i-1),dx(i),mux1,mux2,u(i,j,k),u(i-1,j,k),u(i-2,j,k),u(i-2,j,k+1),u(i-1,j,k+1),u(i,j,k+1),&
                                       u(i+1,j,k+1),u(i+2,j,k+1),u(i+2,j,k),u(i+1,j,k))
    else
      wx = u_y(dx(i-1),dx(i),mux1,mux2,w(i,j,k),w(i-1,j,k),w(i-1,j,k+1),w(i,j,k+1),w(i+1,j,k+1),w(i+1,j,k))
      ux = u_y(dx(i-1),dx(i),mux1,mux2,u(i,j,k),u(i-1,j,k),u(i-1,j,k+1),u(i,j,k+1),u(i+1,j,k+1),u(i+1,j,k))
    endif
    if (3 <= j .and. j <= ny-2 .and. id_visc == 2) then
      vy = u_y(dy(j-1),dy(j),muy1,muy2,v(i,j,k),v(i,j-1,k),v(i,j-2,k),v(i,j-2,k+1),v(i,j-1,k+1),v(i,j,k+1),&
                                       v(i,j+1,k+1),v(i,j+2,k+1),v(i,j+2,k),v(i,j+1,k))
      wy = u_y(dy(j-1),dy(j),muy1,muy2,w(i,j,k),w(i,j-1,k),w(i,j-2,k),w(i,j-2,k+1),w(i,j-1,k+1),w(i,j,k+1),&
                                       w(i,j+1,k+1),w(i,j+2,k+1),w(i,j+2,k),w(i,j+1,k))
    else
      vy = u_y(dy(j-1),dy(j),muy1,muy2,v(i,j,k),v(i,j-1,k),v(i,j-1,k+1),v(i,j,k+1),v(i,j+1,k+1),v(i,j+1,k)) 
      wy = u_y(dy(j-1),dy(j),muy1,muy2,w(i,j,k),w(i,j-1,k),w(i,j-1,k+1),w(i,j,k+1),w(i,j+1,k+1),w(i,j+1,k))  
    endif

    tzx = wx + uz
    tzy = vz + wy
    tzz = 2.d0 * (2.d0 * wz - ux - vy) / 3.d0 

    if (id_turbulence == 1) then
      ! z direction
      muz = 0.5d0 * (mut(i,j,k) + mut(i,j,k+1))
      ! x direction
      mux1 = 0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j,k+1) + mut(i,j,k+1)) 
      mux2 = 0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j,k+1) + mut(i+1,j,k+1))
      ! y direction
      muy1 = 0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i,j-1,k+1) + mut(i,j,k+1))
      muy2 = 0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i,j,k+1) + mut(i,j+1,k+1))

      wz = u_x(dz(k),muz,w(i,j,k),w(i,j,k+1)) 
      uz = u_x(dz(k),muz,u(i,j,k),u(i,j,k+1)) 
      vz = u_x(dz(k),muz,v(i,j,k),v(i,j,k+1)) 
      wx = u_y(dx(i-1),dx(i),mux1,mux2,w(i,j,k),w(i-1,j,k),w(i-1,j,k+1),w(i,j,k+1),w(i+1,j,k+1),w(i+1,j,k))
      ux = u_y(dx(i-1),dx(i),mux1,mux2,u(i,j,k),u(i-1,j,k),u(i-1,j,k+1),u(i,j,k+1),u(i+1,j,k+1),u(i+1,j,k))
      vy = u_y(dy(j-1),dy(j),muy1,muy2,v(i,j,k),v(i,j-1,k),v(i,j-1,k+1),v(i,j,k+1),v(i,j+1,k+1),v(i,j+1,k)) 
      wy = u_y(dy(j-1),dy(j),muy1,muy2,w(i,j,k),w(i,j-1,k),w(i,j-1,k+1),w(i,j,k+1),w(i,j+1,k+1),w(i,j+1,k))  
      ! sgs visc
      tzxsgs = wx + uz
      tzysgs = vz + wy
      tzzsgs = 2.d0 * (2.d0 * wz - ux - vy - (rho(i,j,k) + rho(i,j,k+1)) * qc2(i,j,k)) / 3.d0
      ! sgs kappa
      kappa = kappa + 0.5d0 * (mut(i,j,k) + mut(i,j,k+1)) * Cp / Prt
      ! sgs enthalpy
      H1 = (gamma * p(i,j,k)   / (rho(i,j,k)   * (gamma - 1.d0))) + 0.5d0 * (u(i,j,k)**2   + v(i,j,k)**2   + w(i,j,k)**2)   + qc2(i,j,k)
      H2 = (gamma * p(i,j,k+1) / (rho(i,j,k+1) * (gamma - 1.d0))) + 0.5d0 * (u(i,j,k+1)**2 + v(i,j,k+1)**2 + w(i,j,k+1)**2) + qc2(i,j,k+1)
      Hsgs = -muz * (-H1 + H2) * dz(k) / Prt
    endif
    G(i-offset,j-offset,k-offset+1,2) = G(i-offset,j-offset,k-offset+1,2) - (tzx+tzxsgs)
    G(i-offset,j-offset,k-offset+1,3) = G(i-offset,j-offset,k-offset+1,3) - (tzy+tzysgs)
    G(i-offset,j-offset,k-offset+1,4) = G(i-offset,j-offset,k-offset+1,4) - (tzz+tzzsgs)
    G(i-offset,j-offset,k-offset+1,5) = G(i-offset,j-offset,k-offset+1,5) &
    & -(tzx * um + tzy * vm + tzz * wm + Tz + Hsgs)
  end subroutine calc_Gv
end module calc_visc

