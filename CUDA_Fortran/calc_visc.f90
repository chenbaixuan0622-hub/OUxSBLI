module calc_visc
  use mod_globals, only : id_turbulence
  use calc_Sutherland
  implicit none
contains
  attributes(device) function Smagorinsky(Cs,delta,ux,uy,uz,vx,vy,vz,wx,wy,wz) result(ans)
    real(8), intent(in), value :: Cs, delta, ux, uy, uz, vx, vy, vz, wx, wy, wz
    real(8) :: Dxx, Dyy, Dzz, Dxy, Dyz, Dzx, ans
    Dxx = ux
    Dyy = vy
    Dzz = wz
    Dxy = 0.5d0 * (uy + vx)
    Dyz = 0.5d0 * (vz + wy)
    Dzx = 0.5d0 * (wx + uz)
    ans = ((Cs * delta)**2) * sqrt(2.d0 * (Dxx**2 + Dyy**2 + Dzz**2) + 4.d0 * (Dxy**2 + Dyz**2 + Dzx**2))
  end function Smagorinsky

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

  attributes(global) subroutine calc_Ev(id, nx, ny, nz, dX, dY, dZ, delta, Cs, rho, u, v, w, T, Ev)
    integer(kind=2), intent(in), value :: id
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: dX, dY, dZ, delta, Cs
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w, T
    real(8), intent(out), dimension(nx-1,ny-2,nz-2,5), device :: Ev
    integer i, j, k
    real(8) mux, muy1, muy2, muz1, muz2, kappa
    real(8) ux, uy, uz, vx, vy, vz, wx, wy, wz, txx, txy, txz
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1

    ! x direction
    call calc_mu(T(i,j,k),T(i+1,j,k),mux)
    ux = u_x(dX,mux,u(i,j,k),u(i+1,j,k))
    vx = u_x(dX,mux,v(i,j,k),v(i+1,j,k))
    wx = u_x(dX,mux,w(i,j,k),w(i+1,j,k))

    ! y direction
    call calc_mu(T(i,j,k),T(i,j-1,k),T(i+1,j,k),T(i+1,j-1,k),muy1)
    call calc_mu(T(i,j,k),T(i,j+1,k),T(i+1,j,k),T(i+1,j+1,k),muy2)
    uy = u_y(dY,muy1,muy2,u(i,j,k),u(i,j-1,k),u(i+1,j-1,k),u(i+1,j,k),u(i+1,j+1,k),u(i,j+1,k))
    vy = u_y(dY,muy1,muy2,v(i,j,k),v(i,j-1,k),v(i+1,j-1,k),v(i+1,j,k),v(i+1,j+1,k),v(i,j+1,k))

    ! z direction
    call calc_mu(T(i,j,k),T(i,j,k-1),T(i+1,j,k),T(i+1,j,k-1),muz1)
    call calc_mu(T(i,j,k),T(i,j,k+1),T(i+1,j,k),T(i+1,j,k+1),muz2)
    uz = u_y(dZ,muz1,muz2,u(i,j,k),u(i,j,k-1),u(i+1,j,k-1),u(i+1,j,k),u(i+1,j,k+1),u(i,j,k+1))
    wz = u_y(dZ,muz1,muz2,w(i,j,k),w(i,j,k-1),w(i+1,j,k-1),w(i+1,j,k),w(i+1,j,k+1),w(i,j,k+1))

    txx = 2.d0 * (2.d0 * ux - vy - wz) / 3.d0
    txy = uy + vx
    txz = wx + uz
    call calc_kappa(T(i,j,k),T(i+1,j,k),kappa)
    Ev(i,j-1,k-1,1) = 0.d0
    Ev(i,j-1,k-1,2) = txx
    Ev(i,j-1,k-1,3) = txy
    Ev(i,j-1,k-1,4) = txz
    Ev(i,j-1,k-1,5) = txx * 0.5d0 *(u(i,j,k) + u(i+1,j,k)) + txy * 0.5d0 *(v(i,j,k) + v(i+1,j,k)) &
    & + txz * 0.5d0 *(w(i,j,k) + w(i+1,j,k)) + kappa * (-T(i,j,k) + T(i+1,j,k))
  end subroutine calc_Ev
  
  attributes(global) subroutine calc_Fv(id, nx, ny, nz, dX, dY, dZ, delta, Cs, rho, u, v, w, T, Fv)
    integer(kind=2), intent(in), value :: id
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: dX, dY, dZ, delta, Cs
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w, T
    real(8), intent(out), device :: Fv(nx-2,ny-1,nz-2,5)
    integer i, j, k
    real(8) muy, muz1, muz2, mux1, mux2, kappa
    real(8) ux, uy, uz, vx, vy, vz, wx, wy, wz, tyx, tyy, tyz
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1

    ! y direction
    call calc_mu(T(i,j,k),T(i,j+1,k),muy)
    vy = u_x(dY,muy,v(i,j,k),v(i,j+1,k))
    wy = u_x(dY,muy,w(i,j,k),w(i,j+1,k))
    uy = u_x(dY,muy,u(i,j,k),u(i,j+1,k))

    ! z direction
    call calc_mu(T(i,j,k),T(i,j,k-1),T(i,j+1,k),T(i,j+1,k-1),muz1)
    call calc_mu(T(i,j,k),T(i,j,k+1),T(i,j+1,k),T(i,j+1,k+1),muz2)
    vz = u_y(dZ,muz1,muz2,v(i,j,k),v(i,j,k-1),v(i,j+1,k-1),v(i,j+1,k),v(i,j+1,k+1),v(i,j,k+1)) 
    wz = u_y(dZ,muz1,muz2,w(i,j,k),w(i,j,k-1),w(i,j+1,k-1),w(i,j+1,k),w(i,j+1,k+1),w(i,j,k+1)) 

    ! x direction
    call calc_mu(T(i,j,k),T(i-1,j,k),T(i,j+1,k),T(i-1,j+1,k),mux1)
    call calc_mu(T(i,j,k),T(i+1,j,k),T(i,j+1,k),T(i+1,j+1,k),mux2)
    ux = u_y(dX,mux1,mux2,u(i,j,k),u(i-1,j,k),u(i-1,j+1,k),u(i,j+1,k),u(i+1,j+1,k),u(i+1,j,k)) 
    vx = u_y(dX,mux1,mux2,v(i,j,k),v(i-1,j,k),v(i-1,j+1,k),v(i,j+1,k),v(i+1,j+1,k),v(i+1,j,k)) 

    tyx = uy + vx
    tyy = 2.d0 * (2.d0 * vy - wz - ux) / 3.d0
    tyz = vz + wy
    call calc_kappa(T(i,j,k),T(i,j+1,k),kappa)
    Fv(i-1,j,k-1,1) = 0.d0
    Fv(i-1,j,k-1,2) = tyx
    Fv(i-1,j,k-1,3) = tyy
    Fv(i-1,j,k-1,4) = tyz
    Fv(i-1,j,k-1,5) = tyx * 0.5d0 * (u(i,j,k) + u(i,j+1,k)) + tyy * 0.5d0 * (v(i,j,k) + v(i,j+1,k)) &
    & + tyz * 0.5d0 * (w(i,j,k) + w(i,j+1,k)) + kappa * (-T(i,j,k) + T(i,j+1,k))
  end subroutine calc_Fv
  
  attributes(global) subroutine calc_Gv(id, nx, ny, nz, dX, dY, dZ, delta, Cs, rho, u, v, w, T, Gv)
    integer(kind=2), intent(in), value :: id
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(in), value :: dX, dY, dZ, delta, Cs
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w, T
    real(8), intent(out), device :: Gv(nx-2,ny-2,nz-1,5)
    integer i, j, k
    real(8) muz, mux1, mux2, muy1, muy2, kappa
    real(8) ux, uy, uz, vx, vy, vz, wx, wy, wz, tzx, tzy, tzz
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    
    ! z direction
    call calc_mu(T(i,j,k),T(i,j,k+1),muz)
    wz = u_x(dZ,muz,w(i,j,k),w(i,j,k+1)) 
    uz = u_x(dZ,muz,u(i,j,k),u(i,j,k+1)) 
    vz = u_x(dZ,muz,v(i,j,k),v(i,j,k+1)) 
    
    ! x direction
    call calc_mu(T(i,j,k),T(i-1,j,k),T(i,j,k+1),T(i-1,j,k+1),mux1)
    call calc_mu(T(i,j,k),T(i+1,j,k),T(i,j,k+1),T(i+1,j,k+1),mux1)
    wx = u_y(dX,mux1,mux2,w(i,j,k),w(i-1,j,k),w(i-1,j,k+1),w(i,j,k+1),w(i+1,j,k+1),w(i+1,j,k))
    ux = u_y(dX,mux1,mux2,u(i,j,k),u(i-1,j,k),u(i-1,j,k+1),u(i,j,k+1),u(i+1,j,k+1),u(i+1,j,k))
    
    ! y direction
    call calc_mu(T(i,j,k),T(i,j-1,k),T(i,j,k+1),T(i,j-1,k+1),muy1)
    call calc_mu(T(i,j,k),T(i,j+1,k),T(i,j,k+1),T(i,j+1,k+1),muy1)
    vy = u_y(dY,muy1,muy2,v(i,j,k),v(i,j-1,k),v(i,j-1,k+1),v(i,j,k+1),v(i,j+1,k+1),v(i,j+1,k)) 
    wy = u_y(dY,muy1,muy2,w(i,j,k),w(i,j-1,k),w(i,j-1,k+1),w(i,j,k+1),w(i,j+1,k+1),w(i,j+1,k))  

    tzx = wx + uz
    tzy = vz + wy
    tzz = 2.d0 * (2.d0 * wz - ux - vy) / 3.d0
    call calc_kappa(T(i,j,k),T(i,j,k+1),kappa)
    Gv(i-1,j-1,k,1) = 0.d0
    Gv(i-1,j-1,k,2) = tzx
    Gv(i-1,j-1,k,3) = tzy
    Gv(i-1,j-1,k,4) = tzz
    Gv(i-1,j-1,k,5) = tzx * 0.5d0 * (u(i,j,k) + u(i,j,k+1)) + tzy * 0.5d0 * (v(i,j,k) + v(i,j,k+1)) &
    & + tzz * 0.5d0 * (w(i,j,k) + w(i,j,k+1)) + kappa * (-T(i,j,k) + T(i,j,k+1))
  end subroutine calc_Gv
end module calc_visc

