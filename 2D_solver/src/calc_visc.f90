module calc_visc
  use mod_globals, only : accuracy, offset, id_turbulence, nx, ny
  use calc_Sutherland
  implicit none
contains
  attributes(device) function u_x(d,mu,a1,a2) result(ans)
    real(8), intent(in), value :: d, mu, a1, a2  
    real(8) :: ans
    ans = d * mu * (-a1 + a2)
  end function u_x

  attributes(device) function u_y(d1,d2,mu1,mu2,u1,u2,u3,u4,u5,u6) result(ans)
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
  end function u_y

  attributes(global) subroutine calc_Ev(dx, xix, dy, Jacobian, u, v, T, E)
    real(8), intent(in), dimension(nx), device    :: dx, xix
    real(8), intent(in), dimension(ny), device    :: dy
    real(8), intent(in), dimension(nx,ny), device :: Jacobian, u, v, T
    real(8), intent(inout), device                :: E(nx-accuracy+1,ny-accuracy,4)
    integer i, j
    real(8) mux, muy1, muy2, kappa
    real(8) ux, uy, vx, vy, txx, txy, xixJ, dx1, dy1, dy2
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset

    ! x direction
    call calc_mu(T(i,j),T(i+1,j),mux)
    dx1 = 0.5d0 * (dx(i) + dx(i+1))
    ux = u_x(dx1,mux,u(i,j),u(i+1,j))
    vx = u_x(dx1,mux,v(i,j),v(i+1,j))

    ! y direction
    call calc_mu(T(i,j),T(i,j-1),T(i+1,j),T(i+1,j-1),muy1)
    call calc_mu(T(i,j),T(i,j+1),T(i+1,j),T(i+1,j+1),muy2)
    dy1 = 0.5d0 * (dy(j-1) + dy(j))
    dy2 = 0.5d0 * (dy(j) + dy(j+1))
    uy = u_y(dy1,dy2,muy1,muy2,u(i,j),u(i,j-1),u(i+1,j-1),u(i+1,j),u(i+1,j+1),u(i,j+1))
    vy = u_y(dy1,dy2,muy1,muy2,v(i,j),v(i,j-1),v(i+1,j-1),v(i+1,j),v(i+1,j+1),v(i,j+1))

    xixJ = (xix(i) + xix(i+1)) / (Jacobian(i,j) + Jacobian(i+1,j))
    txx = (2.d0 * (2.d0 * ux - vy) / 3.d0) * xixJ
    txy = (uy + vx) * xixJ
    call calc_kappa(T(i,j),T(i+1,j),kappa)
    E(i-offset,j-offset,2) = E(i-offset,j-offset,2) - txx
    E(i-offset,j-offset,3) = E(i-offset,j-offset,3) - txy
    E(i-offset,j-offset,4) = E(i-offset,j-offset,4) - txx * 0.5d0 *(u(i,j) + u(i+1,j)) - txy * 0.5d0 *(v(i,j) + v(i+1,j)) &
    & - kappa * (-T(i,j) + T(i+1,j)) * dx1 * xixJ
  end subroutine calc_Ev
  
  attributes(global) subroutine calc_Fv(dy, etay, dx, Jacobian, u, v, T, F)
    real(8), intent(in), dimension(ny), device    :: dy, etay
    real(8), intent(in), dimension(nx), device    :: dx
    real(8), intent(in), dimension(nx,ny), device :: Jacobian, u, v, T
    real(8), intent(inout), device                :: F(nx-accuracy,ny-accuracy+1,4)
    integer i, j
    real(8) muy, mux1, mux2, kappa
    real(8) ux, uy, vx, vy, tyx, tyy, etayJ, dy1, dx1, dx2
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset

    ! y direction
    call calc_mu(T(i,j),T(i,j+1),muy)
    dy1 = 0.5d0 * (dy(j) + dy(j+1))
    vy = u_x(dy1,muy,v(i,j),v(i,j+1))
    uy = u_x(dy1,muy,u(i,j),u(i,j+1))

    ! x direction
    call calc_mu(T(i,j),T(i-1,j),T(i,j+1),T(i-1,j+1),mux1)
    call calc_mu(T(i,j),T(i+1,j),T(i,j+1),T(i+1,j+1),mux2)
    dx1 = 0.5d0 * (dx(i-1) + dx(i))
    dx2 = 0.5d0 * (dx(i) + dx(i+1))
    ux = u_y(dx1,dx2,mux1,mux2,u(i,j),u(i-1,j),u(i-1,j+1),u(i,j+1),u(i+1,j+1),u(i+1,j)) 
    vx = u_y(dx1,dx2,mux1,mux2,v(i,j),v(i-1,j),v(i-1,j+1),v(i,j+1),v(i+1,j+1),v(i+1,j)) 

    etayJ = (etay(j) + etay(j+1)) / (Jacobian(i,j) + Jacobian(i,j+1))
    tyx = (uy + vx) * etayJ
    tyy = (2.d0 * (2.d0 * vy - ux) / 3.d0) * etayJ
    call calc_kappa(T(i,j),T(i,j+1),kappa)
    F(i-offset,j-offset,2) = F(i-offset,j-offset,2) - tyx
    F(i-offset,j-offset,3) = F(i-offset,j-offset,3) - tyy
    F(i-offset,j-offset,4) = F(i-offset,j-offset,4) - tyx * 0.5d0 * (u(i,j) + u(i,j+1)) - tyy * 0.5d0 * (v(i,j) + v(i,j+1)) &
    & - kappa * (-T(i,j) + T(i,j+1)) * dy1 * etayJ
  end subroutine calc_Fv
end module calc_visc

