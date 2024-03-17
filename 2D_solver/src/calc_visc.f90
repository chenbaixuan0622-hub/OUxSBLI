module calc_visc
  use mod_globals, only : accuracy, id_turbulence, nx, ny, dxi, dyi
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

  attributes(global) subroutine calc_Ev(u, v, T, E)
    real(8), intent(in), dimension(nx,ny), device :: u, v, T
    real(8), intent(inout), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    integer i, j, offset
    real(8) mux, muy1, muy2, kappa
    real(8) ux, uy, vx, vy, txx, txy
    offset = accuracy / 2
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset - 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset

    ! x direction
    call calc_mu(T(i,j),T(i+1,j),mux)
    ux = u_x(dxi,mux,u(i,j),u(i+1,j))
    vx = u_x(dxi,mux,v(i,j),v(i+1,j))

    ! y direction
    call calc_mu(T(i,j),T(i,j-1),T(i+1,j),T(i+1,j-1),muy1)
    call calc_mu(T(i,j),T(i,j+1),T(i+1,j),T(i+1,j+1),muy2)
    uy = u_y(dyi,muy1,muy2,u(i,j),u(i,j-1),u(i+1,j-1),u(i+1,j),u(i+1,j+1),u(i,j+1))
    vy = u_y(dyi,muy1,muy2,v(i,j),v(i,j-1),v(i+1,j-1),v(i+1,j),v(i+1,j+1),v(i,j+1))

    txx = 2.d0 * (2.d0 * ux - vy) / 3.d0
    txy = uy + vx
    call calc_kappa(T(i,j),T(i+1,j),kappa)
    E(i-offset+1,j-offset,2) = E(i-offset+1,j-offset,2) - txx
    E(i-offset+1,j-offset,3) = E(i-offset+1,j-offset,3) - txy
    E(i-offset+1,j-offset,4) = E(i-offset+1,j-offset,4) - txx * 0.5d0 *(u(i,j) + u(i+1,j)) + txy * 0.5d0 *(v(i,j) + v(i+1,j)) &
    & + kappa * (-T(i,j) + T(i+1,j))
  end subroutine calc_Ev
  
  attributes(global) subroutine calc_Fv(u, v, T, F)
    real(8), intent(in), dimension(nx,ny), device :: u, v, T
    real(8), intent(inout), device :: F(nx-accuracy,ny-accuracy+1,4)
    integer i, j, offset
    real(8) muy, mux1, mux2, kappa
    real(8) ux, uy, vx, vy, tyx, tyy
    offset = accuracy / 2
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y

    ! y direction
    call calc_mu(T(i,j),T(i,j+1),muy)
    vy = u_x(dyi,muy,v(i,j),v(i,j+1))
    uy = u_x(dyi,muy,u(i,j),u(i,j+1))

    ! x direction
    call calc_mu(T(i,j),T(i-1,j),T(i,j+1),T(i-1,j+1),mux1)
    call calc_mu(T(i,j),T(i+1,j),T(i,j+1),T(i+1,j+1),mux2)
    ux = u_y(dxi,mux1,mux2,u(i,j),u(i-1,j),u(i-1,j+1),u(i,j+1),u(i+1,j+1),u(i+1,j)) 
    vx = u_y(dxi,mux1,mux2,v(i,j),v(i-1,j),v(i-1,j+1),v(i,j+1),v(i+1,j+1),v(i+1,j)) 

    tyx = uy + vx
    tyy = 2.d0 * (2.d0 * vy - ux) / 3.d0
    call calc_kappa(T(i,j),T(i,j+1),kappa)
    F(i-offset,j-offset+1,2) = F(i-offset,j-offset+1,2) - tyx
    F(i-offset,j-offset+1,3) = F(i-offset,j-offset+1,3) - tyy
    F(i-offset,j-offset+1,4) = F(i-offset,j-offset+1,4) - tyx * 0.5d0 * (u(i,j) + u(i,j+1)) + tyy * 0.5d0 * (v(i,j) + v(i,j+1)) &
    & + kappa * (-T(i,j) + T(i,j+1))
  end subroutine calc_Fv
end module calc_visc

