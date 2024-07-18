module calc_flux
  use mod_globals, only : id_scheme, id_sensor, id_muscl, accuracy, offset, gamma, threshold
  use calc_keep
  use calc_slau
  use calc_roe
  use calc_hybrid
  implicit none
  interface minmod
    module procedure minmod2, minmod3
  end interface

  interface MUSCL
    module procedure MUSCL3rdnonTVD, MUSCL3rdMinmod, MUSCL4th
  end interface

contains
  attributes(device) function minmod2(x,y) result(ans)
    real(8), intent(in), value :: x, y
    real(8) :: ans, sgn
    sgn = sign(1.d0, x)
    ans = sgn * max(min(abs(x), sgn * y), 0.d0)
  end function minmod2

  attributes(device) function minmod3(x,y,z) result(ans)
    real(8), intent(in), value :: x, y, z
    real(8) :: ans, sgn
    sgn = sign(1.d0, x)
    ans = sgn * max(min(abs(x), sgn * y, sgn * z), 0.d0)
  end function minmod3

  attributes(device) function d33(sigma,d1,d2,d3) result(ans)
    real(8), intent(in), value :: sigma, d1, d2, d3 
    real(8) :: ans, da, db, dc
    da = minmod(d1, sigma * d2, sigma * d3)
    db = minmod(d2, sigma * d1, sigma * d3)
    dc = minmod(d3, sigma * d1, sigma * d2)
    ans = da - 2.d0 * db + dc
  end function d33

  attributes(device) function MUSCL3rdnonTVD(id_tvd,eps,k,a2,a3,d) result(alr)
    integer(kind=2), intent(in), value :: id_tvd
    real(8), intent(in), value         :: eps, k, a2, a3
    real(8), intent(in), device        :: d(3)
    real(8) :: b, alr(2)
    b      = (3.d0 - k) / (1.d0 -k)
    alr(1) = a2 + 0.25d0 * eps * ((1.d0 - k) * d(1) + (1.d0 + K) * d(2))
    alr(2) = a3 - 0.25d0 * eps * ((1.d0 - k) * d(3) + (1.d0 + K) * d(2))
  end function MUSCL3rdnonTVD

  attributes(device) function MUSCL3rdMinmod(id_tvd,eps,k,a2,a3,d) result(alr)
    integer(kind=4), intent(in), value :: id_tvd
    real(8), intent(in), value         :: eps, k, a2, a3
    real(8), intent(in), device        :: d(3)
    real(8) :: b, dt1, dt2, dt3, dt4, alr(2)
    b   = (3.d0 - k) / (1.d0 -k)
    dt1 = minmod(d(1), b * d(2))
    dt2 = minmod(d(2), b * d(1))
    dt3 = minmod(d(3), b * d(2))
    dt4 = minmod(d(2), b * d(3))
    alr(1) = a2 + 0.25d0 * eps * ((1.d0 - k) * dt1 + (1.d0 + K) * dt2)
    alr(2) = a3 - 0.25d0 * eps * ((1.d0 - k) * dt3 + (1.d0 + K) * dt4)
  end function MUSCL3rdMinmod

  attributes(device) function MUSCL4th(id_tvd,eps,k,a2,a3,d) result(alr)
    integer(kind=8), intent(in), value :: id_tvd
    real(8), intent(in), value         :: eps, k, a2, a3
    real(8), intent(in), device        :: d(5)
    real(8) delta1, delta2, delta3, dpl, dml, dpr, dmr, alr(2)
    real(8) :: sigma = 2.d0, w = 4.d0
    delta1 = d(2) - eps * d33(sigma, d(1), d(2), d(3)) / 6.d0
    delta2 = d(3) - eps * d33(sigma, d(2), d(3), d(4)) / 6.d0
    delta3 = d(4) - eps * d33(sigma, d(3), d(4), d(5)) / 6.d0
    dpl    = minmod(delta1, w * delta2)
    dml    = minmod(delta2, w * delta1)
    alr(1) = a2 + (dml + 2.d0 * dpl) / 6.d0
    dpr    = minmod(delta3, w * delta2)
    dmr    = minmod(delta2, w * delta3)
    alr(2) = a3 - (dpr + 2.d0 * dmr) / 6.d0
  end function MUSCL4th

  attributes(device) function delta(points,eps,k,a) result(alr)
    use mod_globals, only : id_tvd
    integer, intent(in), value  :: points
    real(8), intent(in), value  :: eps, k
    real(8), intent(in), device :: a(points)
    real(8) :: alr(2), d(points-1)
    integer i
    do i = 1, points-1
      d(i) = -a(i) + a(i+1)
    enddo
    alr  = MUSCL(id_tvd,eps,k,a(points/2),a(points/2+1),d)
  end function delta

  attributes(device) subroutine calc_points(points,eps1,eps2,eps3,k,rho,p,V,rho2,p2,V2)
    integer, intent(in), value   :: points
    real(8), intent(in), value   :: eps1, eps2, eps3, k
    real(8), intent(in), device  :: rho(points), p(points), V(points,3)
    real(8), intent(out), device :: rho2(2), p2(2), V2(2,3)
    real(8) Vtemp(points)
    rho2(:) = delta(points,eps1,k,rho)
    p2(:)   = delta(points,eps2,k,p)
    Vtemp   = V(:,1)
    V2(:,1) = delta(points,eps3,k,Vtemp)
    Vtemp   = V(:,2)
    V2(:,2) = delta(points,eps3,k,Vtemp)
    Vtemp   = V(:,3)
    V2(:,3) = delta(points,eps3,k,Vtemp)
  end subroutine calc_points

  attributes(device) function Fmuscl(id,points,rho,p,V,Normal,fd) result(F)
    integer, intent(in), value                       :: id, points
    real(8), intent(in), dimension(points), device   :: rho, p
    real(8), intent(in), dimension(points,3), device :: V
    real(8), intent(in), dimension(5), device        :: Normal
    real(8), intent(in), value                       :: fd
    real(8) :: sensor, eps, k = 1.d0 / 3.d0
    real(8) rho2(2), p2(2), V2(2,3), rho4(4), p4(4), V4(4,3), F(5), Fkeep(5), Fslau(5)
    if (points /= 4) then
      rho4(:) = rho(2:5) 
      p4(:)   =   p(2:5)
      V4(:,:) =   V(2:5,:)
    endif
    if (id_sensor == 1) then
      sensor = fd
    elseif (id_sensor == 2) then
      sensor = (1.d0 - Albada(rho4,p4,V4))
    elseif (id_sensor == 3) then
      sensor = fd * (1.d0 - Albada(rho4,p4,V4))
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !KEEP MUSCL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (id_scheme == 2) then
      eps = sensor
      call calc_points(points,eps,eps,eps,k,rho,p,V,rho2,p2,V2)
      F = KEEP2(id,rho2,p2,V2,Normal)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !SLAU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 3) then
      call calc_points(points,1.d0,1.d0,1.d0,k,rho,p,V,rho2,p2,V2)
      F = SLAU(id,rho2,p2,V2,Normal)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !KEEPUP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 4) then
      eps = sensor
      call calc_points(points,1.d0,1.d0,1.d0,0.d0,rho,p,V,rho2,p2,V2)
      F = KEEPUP(id,rho4,p4,V4,rho2,p2,V2,Normal,eps)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Hybrid weighted!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 5) then
      call calc_points(points,1.d0,1.d0,1.d0,k,rho,p,V,rho2,p2,V2)
      Fkeep = KEEP4(id,rho4,p4,V4,Normal)
      Fslau = SLAU(id,rho2,p2,V2,Normal)
      F = (1.d0 - sensor) * Fkeep(:) + sensor * Fslau(:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Hybrid threshold!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 6) then
      if (sensor < threshold) then
        F = KEEP4(id,rho4,p4,V4,Normal)
      else
        call calc_points(points,1.d0,1.d0,1.d0,k,rho,p,V,rho2,p2,V2)
        F = SLAU(id,rho2,p2,V2,Normal)
      endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !KEEP Rho!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 7) then
      call calc_points(points,1.d0,1.d0,1.d0,k,rho,p,V,rho2,p2,V2)
      F = KEEPRho(id,rho4,p4,V4,Normal,rho2,p2,V2,sensor)
    endif
  end function Fmuscl

  attributes(global) subroutine calc_E(nx, ny, nz, rho, u, v, w, p, fd, E)
    use mod_globals, only : id_tvd
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p, fd
    real(8), intent(out), device                      :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    integer i, j, k
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4
    real(8), dimension(4,3) :: V4
    real(8), dimension(6)   :: rho6, p6
    real(8), dimension(6,3) :: V6
    real(8), device         :: Normal(5) = (/0.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(8) fdx
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    fdx = max(fd(i,j,k), fd(i+1,j,k))
    !if (3 <= i .and. i <= nx-3 .and. kind(id_tvd) == 8) then
    !  rho6(:) = rho(i-2:i+3,j,k)
    !  p6(:)   =   p(i-2:i+3,j,k)
    !  V6(:,1) =   u(i-2:i+3,j,k)
    !  V6(:,2) =   v(i-2:i+3,j,k)
    !  V6(:,3) =   w(i-2:i+3,j,k)
    !  E(i,j-offset,k-offset,:) = Fmuscl(1,6,rho6,p6,V6,Normal,fdx)
    if (2 <= i .and. i <= nx-2) then
      rho4(:) = rho(i-1:i+2,j,k)
      p4(:)   =   p(i-1:i+2,j,k)
      V4(:,1) =   u(i-1:i+2,j,k)
      V4(:,2) =   v(i-1:i+2,j,k)
      V4(:,3) =   w(i-1:i+2,j,k)
      if (id_scheme == 1) then
        E(i,j-offset,k-offset,:) = KEEP4(1,rho4,p4,V4,Normal)
      else
        E(i,j-offset,k-offset,:) = Fmuscl(1,4,rho4,p4,V4,Normal,fdx)
      endif
    ! use 3rd-order SLAU at wall
    elseif (i == 1) then
      rho4(:) = (/rho(i,j,k),   rho(i,j,k), rho(i+1,j,k), rho(i+2,j,k)/)
      p4(:)   = (/p(i,j,k),       p(i,j,k),   p(i+1,j,k),   p(i+2,j,k)/)
      V4(:,1) = (/u(i,j,k),       u(i,j,k),   u(i+1,j,k),   u(i+2,j,k)/)
      V4(:,2) = (/v(i,j,k),       v(i,j,k),   v(i+1,j,k),   v(i+2,j,k)/)
      V4(:,3) = (/w(i,j,k),       w(i,j,k),   w(i+1,j,k),   w(i+2,j,k)/)
      call calc_points(4,1.d0,1.d0,1.d0,1.d0/3.d0,rho4,p4,V4,rho2,p2,V2)
      E(i,j-offset,k-offset,:) = SLAU(1,rho2,p2,V2,Normal)
    else
      rho4(:) = (/rho(i-1,j,k), rho(i,j,k), rho(i+1,j,k), rho(i+1,j,k)/)
      p4(:)   = (/p(i-1,j,k),     p(i,j,k),   p(i+1,j,k),   p(i+1,j,k)/)
      V4(:,1) = (/u(i-1,j,k),     u(i,j,k),   u(i+1,j,k),   u(i+1,j,k)/)
      V4(:,2) = (/v(i-1,j,k),     v(i,j,k),   v(i+1,j,k),   v(i+1,j,k)/)
      V4(:,3) = (/w(i-1,j,k),     w(i,j,k),   w(i+1,j,k),   w(i+1,j,k)/)
      call calc_points(4,1.d0,1.d0,1.d0,1.d0/3.d0,rho4,p4,V4,rho2,p2,V2)
      E(i,j-offset,k-offset,:) = SLAU(1,rho2,p2,V2,Normal)
    endif
  end subroutine calc_E

  attributes(global) subroutine calc_F(nx, ny, nz, rho, u, v, w, p, fd, F)
    use mod_globals, only : id_tvd
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p, fd
    real(8), intent(out), device                      :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    integer i, j, k
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4
    real(8), dimension(4,3) :: V4
    real(8), dimension(6)   :: rho6, p6
    real(8), dimension(6,3) :: V6
    real(8), device         :: Normal(5) = (/0.d0, 0.d0, 1.d0, 0.d0, 0.d0/)
    real(8) fdy
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    fdy = max(fd(i,j,k), fd(i,j+1,k))
    !if (3 <= j .and. j <= ny-3 .and. kind(id_tvd) == 8) then
    !  rho6(:) = rho(i,j-2:j+3,k)
    !  p6(:)   =   p(i,j-2:j+3,k)
    !  V6(:,1) =   u(i,j-2:j+3,k)
    !  V6(:,2) =   v(i,j-2:j+3,k)
    !  V6(:,3) =   w(i,j-2:j+3,k)
    !  F(i-offset,j,k-offset,:) = Fmuscl(2,6,rho6,p6,V6,Normal,fdy)
    if (2 <= j .and. j <= ny-2) then
      rho4(:) = rho(i,j-1:j+2,k)
      p4(:)   =   p(i,j-1:j+2,k)
      V4(:,1) =   u(i,j-1:j+2,k)
      V4(:,2) =   v(i,j-1:j+2,k)
      V4(:,3) =   w(i,j-1:j+2,k)
      if (id_scheme == 1) then
        F(i-offset,j,k-offset,:) = KEEP4(2,rho4,p4,V4,Normal)
      else
        F(i-offset,j,k-offset,:) = Fmuscl(2,4,rho4,p4,V4,Normal,fdy)
      endif
    ! use 3rd-order SLAU at wall
    elseif (j == 1) then
      rho4(:) = (/rho(i,j,k),   rho(i,j,k), rho(i,j+1,k), rho(i,j+2,k)/)
      p4(:)   = (/p(i,j,k),       p(i,j,k),   p(i,j+1,k),   p(i,j+2,k)/)
      V4(:,1) = (/u(i,j,k),       u(i,j,k),   u(i,j+1,k),   u(i,j+2,k)/)
      V4(:,2) = (/v(i,j,k),       v(i,j,k),   v(i,j+1,k),   v(i,j+2,k)/)
      V4(:,3) = (/w(i,j,k),       w(i,j,k),   w(i,j+1,k),   w(i,j+2,k)/)
      call calc_points(4,1.d0,1.d0,1.d0,1.d0/3.d0,rho4,p4,V4,rho2,p2,V2)
      F(i-offset,j,k-offset,:) = SLAU(2,rho2,p2,V2,Normal)
    else
      rho4(:) = (/rho(i,j-1,k), rho(i,j,k), rho(i,j+1,k), rho(i,j+1,k)/)
      p4(:)   = (/p(i,j-1,k),     p(i,j,k),   p(i,j+1,k),   p(i,j+1,k)/)
      V4(:,1) = (/u(i,j-1,k),     u(i,j,k),   u(i,j+1,k),   u(i,j+1,k)/)
      V4(:,2) = (/v(i,j-1,k),     v(i,j,k),   v(i,j+1,k),   v(i,j+1,k)/)
      V4(:,3) = (/w(i,j-1,k),     w(i,j,k),   w(i,j+1,k),   w(i,j+1,k)/)
      call calc_points(4,1.d0,1.d0,1.d0,1.d0/3.d0,rho4,p4,V4,rho2,p2,V2)
      F(i-offset,j,k-offset,:) = SLAU(2,rho2,p2,V2,Normal)
    endif
  end subroutine calc_F

  attributes(global) subroutine calc_G(nx, ny, nz, rho, u, v, w, p, fd, G)
    use mod_globals, only : id_tvd
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p, fd
    real(8), intent(out), device                      :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    integer i, j, k
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4
    real(8), dimension(4,3) :: V4
    real(8), dimension(6)   :: rho6, p6
    real(8), dimension(6,3) :: V6
    real(8), device         :: Normal(5) = (/0.d0, 0.d0, 0.d0, 1.d0, 0.d0/)
    real(8) :: fdz
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    fdz = max(fd(i,j,k), fd(i,j,k+1))
    !if (3 <= k .and. k <= nz-3 .and. kind(id_tvd) == 8) then
    !  rho6(:) = rho(i,j,k-2:k+3)
    !  p6(:)   =   p(i,j,k-2:k+3)
    !  V6(:,1) =   u(i,j,k-2:k+3)
    !  V6(:,2) =   v(i,j,k-2:k+3)
    !  V6(:,3) =   w(i,j,k-2:k+3)
    !  G(i-offset,j-offset,k,:) = Fmuscl(3,6,rho6,p6,V6,Normal,fdz)
    if (2 <= k .and. k <= nz-2) then
      rho4(:) = rho(i,j,k-1:k+2)
      p4(:)   =   p(i,j,k-1:k+2)
      V4(:,1) =   u(i,j,k-1:k+2)
      V4(:,2) =   v(i,j,k-1:k+2)
      V4(:,3) =   w(i,j,k-1:k+2)
      if (id_scheme == 1) then
        G(i-offset,j-offset,k,:) = KEEP4(3,rho4,p4,V4,Normal)
      else
        G(i-offset,j-offset,k,:) = Fmuscl(3,4,rho4,p4,V4,Normal,fdz)
      endif
    ! use 3rd-order SLAU at wall
    elseif (k == 1) then
      rho4(:) = (/rho(i,j,k),   rho(i,j,k), rho(i,j,k+1), rho(i,j,k+2)/)
      p4(:)   = (/p(i,j,k),       p(i,j,k),   p(i,j,k+1),   p(i,j,k+2)/)
      V4(:,1) = (/u(i,j,k),       u(i,j,k),   u(i,j,k+1),   u(i,j,k+2)/)
      V4(:,2) = (/v(i,j,k),       v(i,j,k),   v(i,j,k+1),   v(i,j,k+2)/)
      V4(:,3) = (/w(i,j,k),       w(i,j,k),   w(i,j,k+1),   w(i,j,k+2)/)
      call calc_points(4,1.d0,1.d0,1.d0,1.d0/3.d0,rho4,p4,V4,rho2,p2,V2)
      G(i-offset,j-offset,k,:) = SLAU(3,rho2,p2,V2,Normal)
    else
      rho4(:) = (/rho(i,j,k-1), rho(i,j,k), rho(i,j,k+1), rho(i,j,k+1)/)
      p4(:)   = (/p(i,j,k-1),     p(i,j,k),   p(i,j,k+1),   p(i,j,k+1)/)
      V4(:,1) = (/u(i,j,k-1),     u(i,j,k),   u(i,j,k+1),   u(i,j,k+1)/)
      V4(:,2) = (/v(i,j,k-1),     v(i,j,k),   v(i,j,k+1),   v(i,j,k+1)/)
      V4(:,3) = (/w(i,j,k-1),     w(i,j,k),   w(i,j,k+1),   w(i,j,k+1)/)
      call calc_points(4,1.d0,1.d0,1.d0,1.d0/3.d0,rho4,p4,V4,rho2,p2,V2)
      G(i-offset,j-offset,k,:) = SLAU(3,rho2,p2,V2,Normal)
    endif
  end subroutine calc_G
end module calc_flux

