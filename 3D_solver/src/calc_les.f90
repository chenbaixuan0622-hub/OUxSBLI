module calc_les
  use mod_globals, only : id_turbulence, nx, ny, nz, dx, dy, dz, dxi, dyi, dzi
  implicit none
  interface nu_sgs
    module procedure Smagorinsky, selective_mixed_scale
  end interface
contains
  attributes(device) subroutine calc_differential(u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
    real(8), intent(in), dimension(3,3,3), device :: u, v, w
    real(8), intent(out) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    dudx = 0.5d0 * (-u(1,2,2) + u(3,2,2)) * dxi
    dudy = 0.5d0 * (-u(2,1,2) + u(2,3,2)) * dyi
    dudz = 0.5d0 * (-u(2,2,1) + u(2,2,3)) * dzi
    dvdx = 0.5d0 * (-v(1,2,2) + v(3,2,2)) * dxi
    dvdy = 0.5d0 * (-v(2,1,2) + v(2,3,2)) * dyi
    dvdz = 0.5d0 * (-v(2,2,1) + v(2,2,3)) * dzi
    dwdx = 0.5d0 * (-w(1,2,2) + w(3,2,2)) * dxi
    dwdy = 0.5d0 * (-w(2,1,2) + w(2,3,2)) * dyi
    dwdz = 0.5d0 * (-w(2,2,1) + w(2,2,3)) * dzi
  end subroutine calc_differential

  attributes(device) function vorticity(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz) result(omega2)
    real(8), intent(in), value :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    real(8) omega2
    omega2 = (dwdy-dvdz)**2 + (dudz-dwdx)**2 + (dvdx-dudy)**2
  end function vorticity

  attributes(device) function strain_velocity_tensor(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz) result(S2)
    real(8), intent(in), value :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    real(8) S2, S11, S22, S33, S12, S23, S31
    S11 = dudx
    S22 = dvdy
    S33 = dwdz
    S12 = 0.5d0 * (dvdx + dudy)
    S23 = 0.5d0 * (dwdy + dvdz)
    S31 = 0.5d0 * (dudz + dwdx)
    S2 = 2.d0 * (S11**2 + S22**2 + S33**2) + 4.d0 * (S12**2 + S23**2 + S31**2)
  end function strain_velocity_tensor

  attributes(device) function test_filter(x) result(xhat)
    real(8), intent(in), dimension(3,3,3), device :: x
    real(8) xhat
    xhat = x(2,2,2) + 0.25d0 * (x(1,2,2) - 2.d0 * x(2,2,2) + x(3,2,2) &
    & + x(2,1,2) - 2.d0 * x(2,2,2) + x(2,3,2) &
    & + x(2,2,1) - 2.d0 * x(2,2,2) + x(2,2,3))
  end function test_filter

  attributes(device) function stride_filter(x) result(xhat)
    real(8), intent(in), dimension(5,5,5), device :: x
    real(8), intent(out), dimension(3,3,3), device :: xhat
    integer i, j, k
    real(8), dimension(3,3,3) :: xs
    do k = 1, 3
      do j = 1, 3
        do i = 1, 3
          xs(:,:,:) = x(i:i+2,j:j+2,k:k+2)
          xhat(i,j,k) = test_filter(xs)
        enddo
      enddo
    enddo
  end function stride_filter

  attributes(device) function selective_mixed_scale(u,v,w,uhat,vhat,what) result(nut)
    real(8), intent(in), dimension(3,3,3), device :: u, v, w, uhat, vhat, what
    real(8) dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    real(8) dudxh, dudyh, dudzh, dvdxh, dvdyh, dvdzh, dwdxh, dwdyh, dwdzh
    real(8) nut, ftheta, S2, qc2, omega2, omegam2, domega, tan2theta2, a1, a2, delta
    real(8) :: pi = acos(-1.d0)
    real(8) :: Cm = 0.06d0
    real(8) :: alpha = 0.5d0
    call calc_differential(u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
    call calc_differential(uhat,vhat,what,dudxh,dudyh,dudzh,dvdxh,dvdyh,dvdzh,dwdxh,dwdyh,dwdzh)
    omega2 = vorticity(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz) 
    omegam2 = vorticity(dudxh,dudyh,dudzh,dvdxh,dvdyh,dvdzh,dwdxh,dwdyh,dwdzh)
    domega = sqrt(omegam2) - sqrt(omega2)
    tan2theta2 = (2.d0 * sqrt(omegam2) * sqrt(omega2) - omegam2 - omega2 + domega**2) &
    & / (2.d0 * sqrt(omegam2) * sqrt(omega2) + omegam2 + omega2 - domega**2) 
    if (2.d0 * atan(sqrt(tan2theta2)) > 20.d0 * pi / 180.d0) then
      ftheta = 1.d0
    else
      ftheta = (tan2theta2 / (tan(10.d0 * pi / 180.d0)**2))**2
    endif
    S2 = strain_velocity_tensor(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
    qc2 = 0.5d0 * ((u(2,2,2) - uhat(2,2,2))**2 + (v(2,2,2) - vhat(2,2,2))**2 + (w(2,2,2) - what(2,2,2))**2)
    ! aspect ratio
    a1 = 1.d0
    a2 = 1.d0
    delta = cosh(sqrt(4.d0 * ((log(a1))**2 - log(a1) * log(a2) + (log(a2))**2) / 27.d0)) * (dx * dy * dz)**(1.d0/3.d0) 
    nut = ftheta * Cm * (S2**(0.5d0 * alpha)) * (qc2**(0.5d0 * (1.d0 - alpha))) * (delta**(1.d0 - alpha))
  end function selective_mixed_scale

  attributes(device) function Smagorinsky(u,v,w) result(nut)
    real(8), intent(in), dimension(3,3,3), device :: u, v, w
    real(8) dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    real(8) nut, S2, delta
    real(8) :: Cs = 0.1d0
    call calc_differential(u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
    S2 = strain_velocity_tensor(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
    delta = (dx * dy * dz)**(1.d0/3.d0)
    nut = ((Cs * delta)**2) * sqrt(S2)
  end function Smagorinsky

  attributes(global) subroutine calc_mut(rho,u,v,w,mut)
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w
    real(8), intent(out), dimension(nx,ny,nz), device :: mut
    integer i, j, k
    real(8), dimension(3,3,3), device :: us, vs, ws, uhat, vhat, what
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    !us = u(i-1:i+1,j-1:j+1,w-1:w+1)
    !vs = v(i-1:i+1,j-1:j+1,w-1:w+1)
    !ws = w(i-1:i+1,j-1:j+1,w-1:w+1)
    us() = u()
    us() = u()
    us() = u()
    if (id_turbulence == 1) then
      mut(i,j,k) = rho(i,j,k) * Smagorinsky(us,vs,ws) 
    else
      uhat(:,:,:) = stride_filter(u(i-2:i+2,j-2:j+2,k-2:k+2))
      vhat(:,:,:) = stride_filter(v(i-2:i+2,j-2:j+2,k-2:k+2))
      what(:,:,:) = stride_filter(w(i-2:i+2,j-2:j+2,k-2:k+2))
      mut(i,j,k) = rho(i,j,k) * selective_mixed_scale(us,vs,ws,uhat,vhat,what)
    endif
  end subroutine calc_mut
end module calc_les

