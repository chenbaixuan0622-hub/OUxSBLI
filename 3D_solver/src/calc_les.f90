module calc_les
  implicit none
contains
  attributes(device) subroutine calc_differential(u,v,w,dx,dy,dz,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
    real(8), intent(in), dimension(3,3,3), device :: u, v, w
    real(8), intent(in), value                    :: dx, dy, dz ! 1 /dx, 1 / dy, 1 / dz
    real(8), intent(out) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    dudx = 0.5d0 * (-u(1,2,2) + u(3,2,2)) * dx
    dudy = 0.5d0 * (-u(2,1,2) + u(2,3,2)) * dy
    dudz = 0.5d0 * (-u(2,2,1) + u(2,2,3)) * dz
    dvdx = 0.5d0 * (-v(1,2,2) + v(3,2,2)) * dx
    dvdy = 0.5d0 * (-v(2,1,2) + v(2,3,2)) * dy
    dvdz = 0.5d0 * (-v(2,2,1) + v(2,2,3)) * dz
    dwdx = 0.5d0 * (-w(1,2,2) + w(3,2,2)) * dx
    dwdy = 0.5d0 * (-w(2,1,2) + w(2,3,2)) * dy
    dwdz = 0.5d0 * (-w(2,2,1) + w(2,2,3)) * dz
  end subroutine calc_differential

  attributes(device) function vorticity(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz) result(omega2)
    real(8), intent(in), value :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    real(8) omega2
    omega2 = (dwdy-dvdz)**2 + (dudz-dwdx)**2 + (dvdx-dudy)**2
  end function vorticity

  attributes(device) function strain_velocity_tensor(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz) result(S2)
    real(8), intent(in), value :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    real(8) S2
    S2 = 2.d0 * (dudx**2 + dvdy**2 + dwdz**2) &
              + ((dvdx + dudy)**2 + (dwdy + dvdz)**2 + (dudz + dwdx)**2)
  end function strain_velocity_tensor

  attributes(device) function SMS(u,v,w,uh,vh,wh,dx,dy,dz) result(nut)
    real(8), intent(in), dimension(3,3,3), device :: u, v, w, uh, vh, wh
    real(8), intent(in), value                    :: dx, dy, dz ! 1 / dx, 1 / dy, 1 / dz
    real(8) dudx,  dudy,  dudz,  dvdx,  dvdy,  dvdz,  dwdx,  dwdy,  dwdz
    real(8) dudxh, dudyh, dudzh, dvdxh, dvdyh, dvdzh, dwdxh, dwdyh, dwdzh
    real(8) nut, ftheta, S2, qc2, omega2, omegam2, domega, tan2theta2, a1, a2, delta
    real(8) :: pi    = acos(-1.d0)
    real(8) :: Cm    = 0.06d0
    real(8) :: alpha = 0.5d0
    call calc_differential(u, v, w, dx, dy, dz, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
    call calc_differential(uh,vh,wh,dx, dy, dz, dudxh,dudyh,dudzh,dvdxh,dvdyh,dvdzh,dwdxh,dwdyh,dwdzh)
    omega2  = vorticity(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz) 
    omegam2 = vorticity(dudxh,dudyh,dudzh,dvdxh,dvdyh,dvdzh,dwdxh,dwdyh,dwdzh)
    domega  = sqrt(omegam2) - sqrt(omega2)
    tan2theta2 = (2.d0 * sqrt(omegam2) * sqrt(omega2) - omegam2 - omega2 + domega**2) &
               / (2.d0 * sqrt(omegam2) * sqrt(omega2) + omegam2 + omega2 - domega**2) 
    if (2.d0 * atan(sqrt(tan2theta2)) > 20.d0 * pi / 180.d0) then
      ftheta = 1.d0
    else
      ftheta = (tan2theta2 / (tan(10.d0 * pi / 180.d0)**2))**2
    endif
    S2 = strain_velocity_tensor(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
    qc2 = 0.5d0 * ((u(2,2,2) - uh(2,2,2))**2 + (v(2,2,2) - vh(2,2,2))**2 + (w(2,2,2) - wh(2,2,2))**2)
    ! aspect ratio
    a1 = min(dx, dy, dz) / dx
    a2 = min(dx, dy, dz) / dy
    delta = cosh(sqrt(4.d0 * ((log(a1))**2 - log(a1) * log(a2) + (log(a2))**2) / 27.d0)) * (dx * dy * dz)**(-1.d0/3.d0) 
    nut = ftheta * Cm * (S2**(0.5d0 * alpha)) * (qc2**(0.5d0 * (1.d0 - alpha))) * (delta**(1.d0 - alpha))
  end function SMS

  attributes(device) function test_filter(x) result(xh)
    real(8), intent(in), dimension(3,3,3), device :: x
    real(8) xh
    xh = x(2,2,2) + 0.25d0 * (x(1,2,2) - 2.d0 * x(2,2,2) + x(3,2,2) &
                            + x(2,1,2) - 2.d0 * x(2,2,2) + x(2,3,2) &
                            + x(2,2,1) - 2.d0 * x(2,2,2) + x(2,2,3))
  end function test_filter

  attributes(device) function stride_filter(x) result(xh)
    real(8), intent(in), dimension(5,5,5), device   :: x
    real(8), intent(out), dimension(3,3,3), device  :: xh
    integer i, j, k
    real(8), dimension(3,3,3) :: xs
    do k = 1, 3
      do j = 1, 3
        do i = 1, 3
          xs(:,:,:) = x(i:i+2,j:j+2,k:k+2)
          xh(i,j,k) = test_filter(xs)
    enddo;enddo;enddo
  end function stride_filter

  attributes(global) subroutine calc_mut(nx,ny,nz,dx,dy,dz,rho,u,v,w,mut,qc2)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), device                       :: dx(nx-1), dy(ny-1), dz(nz-1)
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w
    real(8), intent(out), dimension(nx,ny,nz), device :: mut, qc2
    integer i, j, k
    real(8), dimension(3,3,3), device :: u3, v3, w3, uh, vh, wh
    real(8), dimension(5,5,5), device :: u5, v5, w5
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    u3 = u(i-1:i+1,j-1:j+1,k-1:k+1)
    v3 = v(i-1:i+1,j-1:j+1,k-1:k+1)
    w3 = w(i-1:i+1,j-1:j+1,k-1:k+1)
    if (3<=i .and. i<=nx-2 .and. 3<=j .and. j<=ny-2 .and. 3<=k .and. k <=nz-2) then
      u5 = u(i-2:i+2,j-2:j+2,k-2:k+2)
      v5 = v(i-2:i+2,j-2:j+2,k-2:k+2)
      w5 = w(i-2:i+2,j-2:j+2,k-2:k+2)
      uh = stride_filter(u5)
      vh = stride_filter(v5)
      wh = stride_filter(w5)
    else
      uh = u3
      vh = v3
      wh = w3
    endif
    mut(i,j,k) = rho(i,j,k) * SMS(u3,v3,w3,uh,vh,wh,dx(i),dy(j),dz(k))
    qc2(i,j,k) = 0.5d0 * ((u(i,j,k) - uh(2,2,2))**2 + (v(i,j,k) - vh(2,2,2))**2 + (w(i,j,k) - wh(2,2,2))**2)
  end subroutine calc_mut
end module calc_les

