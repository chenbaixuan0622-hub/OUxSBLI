module calc_flux
  use mod_globals, only : id_scheme, id_sensor, id_muscl, accuracy, offset, gamma, threshold
  use calc_keep
  use calc_slau
  use calc_roe
  use calc_hybrid
  use calc_muscl
  implicit none
contains
  attributes(device) function Fmuscl4(id,rho,p,V,Normal,fd) result(F)
    use mod_globals, only : id_slau
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(4), device   :: rho, p
    real(8), intent(in), dimension(4,3), device :: V
    real(8), intent(in), dimension(5), device   :: Normal
    real(8), intent(in), value                  :: fd
    real(8) :: sensor, eps, k = 1.d0 / 3.d0, wiggle
    real(8) rho2(2), p2(2), V2(2,3), F(5)
    if (id_sensor == 1) then
      sensor = fd
    elseif (id_sensor == 2) then
      sensor = (1.d0 - Albada(rho,p,V))
    elseif (id_sensor == 3) then
      sensor = fd * (1.d0 - Albada(rho,p,V))
    endif
    ! for HR-SLAU
    wiggle = wiggle_detector(p)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !KEEP MUSCL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (id_scheme == 2) then
      eps = sensor
      call calc_4points(eps,eps,eps,k,rho,p,V,rho2,p2,V2)
      F = KEEP2(id,rho2,p2,V2,Normal)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !SLAU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 3) then
      call calc_4points(1.d0,1.d0,1.d0,k,rho,p,V,rho2,p2,V2)
      F = SLAU(id_slau,id,rho2,p2,V2,Normal,wiggle,sensor)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !KEEPUP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 4) then
      eps = sensor
      call calc_4points(1.d0,1.d0,1.d0,0.d0,rho,p,V,rho2,p2,V2)
      F = KEEPUP(id,rho,p,V,rho2,p2,V2,Normal,eps)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Hybrid weighted!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 5) then
      call calc_4points(1.d0,1.d0,1.d0,k,rho,p,V,rho2,p2,V2)
      F = (1.d0 - sensor) * KEEP4(id,rho,p,V,Normal) &
          + sensor * SLAU(id_slau,id,rho2,p2,V2,Normal,wiggle,sensor)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Hybrid threshold!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 6) then
      if (sensor < threshold) then
        F = KEEP4(id,rho,p,V,Normal)
      else
        call calc_4points(1.d0,1.d0,1.d0,k,rho,p,V,rho2,p2,V2)
        F = SLAU(id_slau,id,rho2,p2,V2,Normal,wiggle,sensor)
      endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Hybrid Sigmoid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 7) then
      call calc_4points(1.d0,1.d0,1.d0,k,rho,p,V,rho2,p2,V2)
      F = (1.d0 - sigmoid(sensor)) * KEEP4(id,rho,p,V,Normal) &
          + sigmoid(sensor) * SLAU(id_slau,id,rho2,p2,V2,Normal,wiggle,sensor)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !KEEP Rho!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 8) then
      call calc_4points(1.d0,1.d0,1.d0,k,rho,p,V,rho2,p2,V2)
      F = KEEPRho(id,rho,p,V,Normal,rho2,p2,V2,sensor)
    endif
  end function Fmuscl4

  attributes(device) function Fmuscl6(id,rho,p,V,Normal,fd) result(F)
    use mod_globals, only : id_slau
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(6), device   :: rho, p
    real(8), intent(in), dimension(6,3), device :: V
    real(8), intent(in), dimension(5), device   :: Normal
    real(8), intent(in), value                  :: fd
    real(8) :: sensor, eps, k = 1.d0 / 3.d0, wiggle
    real(8) rho2(2), p2(2), V2(2,3), rho4(4), p4(4), V4(4,3), F(5)
    rho4(:) = rho(2:5)
    p4(:)   =   p(2:5)
    V4(:,:) =   V(2:5,:)
    if (id_sensor == 1) then
      sensor = fd
    elseif (id_sensor == 2) then
      sensor = (1.d0 - Albada(rho4,p4,V4))
    elseif (id_sensor == 3) then
      sensor = fd * (1.d0 - Albada(rho4,p4,V4))
    endif
    ! for HR-SLAU
    wiggle = wiggle_detector(p4)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !KEEP MUSCL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (id_scheme == 2) then
      eps = sensor
      call calc_6points(eps,eps,eps,k,rho,p,V,rho2,p2,V2)
      F = KEEP2(id,rho2,p2,V2,Normal)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !SLAU!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 3) then
      call calc_6points(1.d0,1.d0,1.d0,k,rho,p,V,rho2,p2,V2)
      F = SLAU(id_slau,id,rho2,p2,V2,Normal,wiggle,sensor)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !KEEPUP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 4) then
      eps = sensor
      call calc_6points(1.d0,1.d0,1.d0,0.d0,rho,p,V,rho2,p2,V2)
      F = KEEPUP(id,rho4,p4,V4,rho2,p2,V2,Normal,eps)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Hybrid weighted!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 5) then
      call calc_6points(1.d0,1.d0,1.d0,k,rho,p,V,rho2,p2,V2)
      F = (1.d0 - sensor) * KEEP4(id,rho4,p4,V4,Normal) &
          + sensor * SLAU(id_slau,id,rho2,p2,V2,Normal,wiggle,sensor)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Hybrid threshold!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 6) then
      if (sensor < threshold) then
        F = KEEP4(id,rho4,p4,V4,Normal)
      else
        call calc_6points(1.d0,1.d0,1.d0,k,rho,p,V,rho2,p2,V2)
        F = SLAU(id_slau,id,rho2,p2,V2,Normal,wiggle,sensor)
      endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Hybrid Sigmoid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 7) then
      call calc_6points(1.d0,1.d0,1.d0,k,rho,p,V,rho2,p2,V2)
      F = (1.d0 - sigmoid(sensor)) * KEEP4(id,rho4,p4,V4,Normal) &
          + sigmoid(sensor) * SLAU(id_slau,id,rho2,p2,V2,Normal,wiggle,sensor)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !KEEP Rho!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif (id_scheme == 8) then
      call calc_6points(1.d0,1.d0,1.d0,k,rho,p,V,rho2,p2,V2)
      F = KEEPRho(id,rho4,p4,V4,Normal,rho2,p2,V2,sensor)
    endif
  end function Fmuscl6

  attributes(global) subroutine calc_E(nx, ny, nz, rho, u, v, w, p, fd, E)
    use mod_globals, only : id_tvd
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p, fd
    real(8), intent(out), device                      :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    integer i, j, k
    integer(kind=2) id_slau_wall
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
    if (3 <= i .and. i <= nx-3 .and. kind(id_tvd) == 8 .and. id_scheme <= 8) then
      rho6(:) = rho(i-2:i+3,j,k)
      p6(:)   =   p(i-2:i+3,j,k)
      V6(:,1) =   u(i-2:i+3,j,k)
      V6(:,2) =   v(i-2:i+3,j,k)
      V6(:,3) =   w(i-2:i+3,j,k)
      E(i,j-offset,k-offset,:) = Fmuscl6(1,rho6,p6,V6,Normal,fdx)
    elseif (2 <= i .and. i <= nx-2 .and. id_scheme <= 8) then
      rho4(:) = rho(i-1:i+2,j,k)
      p4(:)   =   p(i-1:i+2,j,k)
      V4(:,1) =   u(i-1:i+2,j,k)
      V4(:,2) =   v(i-1:i+2,j,k)
      V4(:,3) =   w(i-1:i+2,j,k)
      if (id_scheme == 1) then
        E(i,j-offset,k-offset,:) = KEEP4(1,rho4,p4,V4,Normal)
      else
        E(i,j-offset,k-offset,:) = Fmuscl4(1,rho4,p4,V4,Normal,fdx)
      endif
    elseif (id_scheme == 9) then
      rho2(:) = rho(i:i+1,j,k)
      p2(:)   =   p(i:i+1,j,k)
      V2(:,1) =   u(i:i+1,j,k)
      V2(:,2) =   v(i:i+1,j,k)
      V2(:,3) =   w(i:i+1,j,k)
      E(i,j-offset,k-offset,:) = KEEP2(1,rho2,p2,V2,Normal)
    ! use 3rd-order SLAU at wall
    elseif (i == 1) then
      rho4(:) = (/rho(i,j,k),   rho(i,j,k), rho(i+1,j,k), rho(i+2,j,k)/)
      p4(:)   = (/p(i,j,k),       p(i,j,k),   p(i+1,j,k),   p(i+2,j,k)/)
      V4(:,1) = (/u(i,j,k),       u(i,j,k),   u(i+1,j,k),   u(i+2,j,k)/)
      V4(:,2) = (/v(i,j,k),       v(i,j,k),   v(i+1,j,k),   v(i+2,j,k)/)
      V4(:,3) = (/w(i,j,k),       w(i,j,k),   w(i+1,j,k),   w(i+2,j,k)/)
      call calc_4points(1.d0,1.d0,1.d0,1.d0/3.d0,rho4,p4,V4,rho2,p2,V2)
      E(i,j-offset,k-offset,:) = SLAU(id_slau_wall,1,rho2,p2,V2,Normal)
    else
      rho4(:) = (/rho(i-1,j,k), rho(i,j,k), rho(i+1,j,k), rho(i+1,j,k)/)
      p4(:)   = (/p(i-1,j,k),     p(i,j,k),   p(i+1,j,k),   p(i+1,j,k)/)
      V4(:,1) = (/u(i-1,j,k),     u(i,j,k),   u(i+1,j,k),   u(i+1,j,k)/)
      V4(:,2) = (/v(i-1,j,k),     v(i,j,k),   v(i+1,j,k),   v(i+1,j,k)/)
      V4(:,3) = (/w(i-1,j,k),     w(i,j,k),   w(i+1,j,k),   w(i+1,j,k)/)
      call calc_4points(1.d0,1.d0,1.d0,1.d0/3.d0,rho4,p4,V4,rho2,p2,V2)
      E(i,j-offset,k-offset,:) = SLAU(id_slau_wall,1,rho2,p2,V2,Normal)
    endif
  end subroutine calc_E

  attributes(global) subroutine calc_F(nx, ny, nz, rho, u, v, w, p, fd, F)
    use mod_globals, only : id_tvd
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p, fd
    real(8), intent(out), device                      :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    integer i, j, k
    integer(kind=2) id_slau_wall
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
    if (3 <= j .and. j <= ny-3 .and. kind(id_tvd) == 8 .and. id_scheme <= 8) then
      rho6(:) = rho(i,j-2:j+3,k)
      p6(:)   =   p(i,j-2:j+3,k)
      V6(:,1) =   u(i,j-2:j+3,k)
      V6(:,2) =   v(i,j-2:j+3,k)
      V6(:,3) =   w(i,j-2:j+3,k)
      F(i-offset,j,k-offset,:) = Fmuscl6(2,rho6,p6,V6,Normal,fdy)
    elseif (2 <= j .and. j <= ny-2 .and. id_scheme <= 8) then
      rho4(:) = rho(i,j-1:j+2,k)
      p4(:)   =   p(i,j-1:j+2,k)
      V4(:,1) =   u(i,j-1:j+2,k)
      V4(:,2) =   v(i,j-1:j+2,k)
      V4(:,3) =   w(i,j-1:j+2,k)
      if (id_scheme == 1) then
        F(i-offset,j,k-offset,:) = KEEP4(2,rho4,p4,V4,Normal)
      else
        F(i-offset,j,k-offset,:) = Fmuscl4(2,rho4,p4,V4,Normal,fdy)
      endif
    elseif (id_scheme == 9) then
      rho2(:) = rho(i,j:j+1,k)
      p2(:)   =   p(i,j:j+1,k)
      V2(:,1) =   u(i,j:j+1,k)
      V2(:,2) =   v(i,j:j+1,k)
      V2(:,3) =   w(i,j:j+1,k)
      F(i-offset,j,k-offset,:) = KEEP2(2,rho2,p2,V2,Normal)
    ! use 3rd-order SLAU at wall
    elseif (j == 1) then
      rho4(:) = (/rho(i,j,k),   rho(i,j,k), rho(i,j+1,k), rho(i,j+2,k)/)
      p4(:)   = (/p(i,j,k),       p(i,j,k),   p(i,j+1,k),   p(i,j+2,k)/)
      V4(:,1) = (/u(i,j,k),       u(i,j,k),   u(i,j+1,k),   u(i,j+2,k)/)
      V4(:,2) = (/v(i,j,k),       v(i,j,k),   v(i,j+1,k),   v(i,j+2,k)/)
      V4(:,3) = (/w(i,j,k),       w(i,j,k),   w(i,j+1,k),   w(i,j+2,k)/)
      call calc_4points(1.d0,1.d0,1.d0,1.d0/3.d0,rho4,p4,V4,rho2,p2,V2)
      F(i-offset,j,k-offset,:) = SLAU(id_slau_wall,2,rho2,p2,V2,Normal)
    else
      rho4(:) = (/rho(i,j-1,k), rho(i,j,k), rho(i,j+1,k), rho(i,j+1,k)/)
      p4(:)   = (/p(i,j-1,k),     p(i,j,k),   p(i,j+1,k),   p(i,j+1,k)/)
      V4(:,1) = (/u(i,j-1,k),     u(i,j,k),   u(i,j+1,k),   u(i,j+1,k)/)
      V4(:,2) = (/v(i,j-1,k),     v(i,j,k),   v(i,j+1,k),   v(i,j+1,k)/)
      V4(:,3) = (/w(i,j-1,k),     w(i,j,k),   w(i,j+1,k),   w(i,j+1,k)/)
      call calc_4points(1.d0,1.d0,1.d0,1.d0/3.d0,rho4,p4,V4,rho2,p2,V2)
      F(i-offset,j,k-offset,:) = SLAU(id_slau_wall,2,rho2,p2,V2,Normal)
    endif
  end subroutine calc_F

  attributes(global) subroutine calc_G(nx, ny, nz, rho, u, v, w, p, fd, G)
    use mod_globals, only : id_tvd
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p, fd
    real(8), intent(out), device                      :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    integer i, j, k
    integer(kind=2) id_slau_wall
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
    if (3 <= k .and. k <= nz-3 .and. kind(id_tvd) == 8 .and. id_scheme <= 8) then
      rho6(:) = rho(i,j,k-2:k+3)
      p6(:)   =   p(i,j,k-2:k+3)
      V6(:,1) =   u(i,j,k-2:k+3)
      V6(:,2) =   v(i,j,k-2:k+3)
      V6(:,3) =   w(i,j,k-2:k+3)
      G(i-offset,j-offset,k,:) = Fmuscl6(3,rho6,p6,V6,Normal,fdz)
    elseif (2 <= k .and. k <= nz-2 .and. id_scheme <= 8) then
      rho4(:) = rho(i,j,k-1:k+2)
      p4(:)   =   p(i,j,k-1:k+2)
      V4(:,1) =   u(i,j,k-1:k+2)
      V4(:,2) =   v(i,j,k-1:k+2)
      V4(:,3) =   w(i,j,k-1:k+2)
      if (id_scheme == 1) then
        G(i-offset,j-offset,k,:) = KEEP4(3,rho4,p4,V4,Normal)
      else
        G(i-offset,j-offset,k,:) = Fmuscl4(3,rho4,p4,V4,Normal,fdz)
      endif
    elseif (id_scheme == 9) then
      rho2(:) = rho(i,j,k:k+1)
      p2(:)   =   p(i,j,k:k+1)
      V2(:,1) =   u(i,j,k:k+1)
      V2(:,2) =   v(i,j,k:k+1)
      V2(:,3) =   w(i,j,k:k+1)
      G(i-offset,j-offset,k,:) = KEEP2(3,rho2,p2,V2,Normal)
    ! use 3rd-order SLAU at wall
    elseif (k == 1) then
      rho4(:) = (/rho(i,j,k),   rho(i,j,k), rho(i,j,k+1), rho(i,j,k+2)/)
      p4(:)   = (/p(i,j,k),       p(i,j,k),   p(i,j,k+1),   p(i,j,k+2)/)
      V4(:,1) = (/u(i,j,k),       u(i,j,k),   u(i,j,k+1),   u(i,j,k+2)/)
      V4(:,2) = (/v(i,j,k),       v(i,j,k),   v(i,j,k+1),   v(i,j,k+2)/)
      V4(:,3) = (/w(i,j,k),       w(i,j,k),   w(i,j,k+1),   w(i,j,k+2)/)
      call calc_4points(1.d0,1.d0,1.d0,1.d0/3.d0,rho4,p4,V4,rho2,p2,V2)
      G(i-offset,j-offset,k,:) = SLAU(id_slau_wall,3,rho2,p2,V2,Normal)
    else
      rho4(:) = (/rho(i,j,k-1), rho(i,j,k), rho(i,j,k+1), rho(i,j,k+1)/)
      p4(:)   = (/p(i,j,k-1),     p(i,j,k),   p(i,j,k+1),   p(i,j,k+1)/)
      V4(:,1) = (/u(i,j,k-1),     u(i,j,k),   u(i,j,k+1),   u(i,j,k+1)/)
      V4(:,2) = (/v(i,j,k-1),     v(i,j,k),   v(i,j,k+1),   v(i,j,k+1)/)
      V4(:,3) = (/w(i,j,k-1),     w(i,j,k),   w(i,j,k+1),   w(i,j,k+1)/)
      call calc_4points(1.d0,1.d0,1.d0,1.d0/3.d0,rho4,p4,V4,rho2,p2,V2)
      G(i-offset,j-offset,k,:) = SLAU(id_slau_wall,3,rho2,p2,V2,Normal)
    endif
  end subroutine calc_G
end module calc_flux

