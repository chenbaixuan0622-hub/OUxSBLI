module calc_flux
  use mod_globals, only : id_scheme, id_sensor, id_muscl, accuracy, offset, gamma, threshold
  use calc_keep
  use calc_slau
  use calc_roe
  use calc_hybrid
  use calc_muscl
  implicit none
  interface flux6
    module procedure flux_KEEP6, flux_SLAU6, flux_Weighted6, flux_Threshold6
  end interface flux6

  interface flux4
    module procedure flux_KEEP4, flux_SLAU4, flux_Weighted4, flux_Threshold4
  end interface flux4

  interface flux2
    module procedure flux_KEEP2, flux_SLAU2, flux_Weighted2, flux_Threshold2
  end interface flux2
contains
  attributes(device) function flux_KEEP6(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    integer(kind=2), intent(in), value          :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(6), device   :: rho, p
    real(8), intent(in), dimension(6,3), device :: V
    real(8), intent(in), dimension(5), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(5)
    F = KEEP6(id,rho,p,V,Normal)
  end function flux_KEEP6

  attributes(device) function flux_KEEP4(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    integer(kind=2), intent(in), value          :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(4), device   :: rho, p
    real(8), intent(in), dimension(4,3), device :: V
    real(8), intent(in), dimension(5), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(5)
    F = KEEP4(id,rho,p,V,Normal)
  end function flux_KEEP4
  
  attributes(device) function flux_KEEP2(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    integer(kind=2), intent(in), value          :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(2), device   :: rho, p
    real(8), intent(in), dimension(2,3), device :: V
    real(8), intent(in), dimension(5), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(5)
    F = KEEP2(id,rho,p,V,Normal)
  end function flux_KEEP2

  attributes(device) function flux_SLAU6(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    use mod_globals, only : id_slau
    real(kind=2), intent(in), value             :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(6), device   :: rho, p
    real(8), intent(in), dimension(6,3), device :: V
    real(8), intent(in), dimension(5), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) wiggle, rho2(2), p2(2), V2(2,3), p4(4), F(5)
    p4(:)  = p(2:5)
    wiggle = wiggle_detector(p4)
    call calc_6points(sensor,rho,p,V,rho2,p2,V2)
    F = SLAU(id_slau,id,rho2,p2,V2,Normal,wiggle,sensor)
  end function flux_SLAU6

  attributes(device) function flux_SLAU4(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    use mod_globals, only : id_slau
    real(kind=2), intent(in), value             :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(4), device   :: rho, p
    real(8), intent(in), dimension(4,3), device :: V
    real(8), intent(in), dimension(5), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) wiggle, rho2(2), p2(2), V2(2,3), F(5)
    wiggle = wiggle_detector(p)
    call calc_4points(1.d0,1.d0,1.d0/3.d0,sensor,rho,p,V,rho2,p2,V2)
    F = SLAU(id_slau,id,rho2,p2,V2,Normal,wiggle,sensor)
  end function flux_SLAU4

  attributes(device) function flux_SLAU2(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    use mod_globals, only : id_slau
    real(kind=2), intent(in), value             :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(2), device   :: rho, p
    real(8), intent(in), dimension(2,3), device :: V
    real(8), intent(in), dimension(5), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(5)
    F = SLAU(id_slau,id,rho,p,V,Normal,1.d0,sensor)
  end function flux_SLAU2

  attributes(device) function flux_Weighted6(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    real(4), intent(in), value                  :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(6), device   :: rho, p
    real(8), intent(in), dimension(6,3), device :: V
    real(8), intent(in), dimension(5), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(5)
    real(2) slau
    F = (1.d0 - sensor) * KEEP6(id,rho,p,V,Normal) &
        + sensor * flux_SLAU6(slau,id,rho,p,V,Normal,sensor)
  end function flux_Weighted6

  attributes(device) function flux_Weighted4(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    real(4), intent(in), value                  :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(4), device   :: rho, p
    real(8), intent(in), dimension(4,3), device :: V
    real(8), intent(in), dimension(5), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(5)
    real(2) slau
    F = (1.d0 - sensor) * KEEP4(id,rho,p,V,Normal) &
        + sensor * flux_SLAU4(slau,id,rho,p,V,Normal,sensor)
  end function flux_Weighted4

  attributes(device) function flux_Weighted2(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    use mod_globals, only : id_slau
    real(4), intent(in), value                 :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(2), device   :: rho, p
    real(8), intent(in), dimension(2,3), device :: V
    real(8), intent(in), dimension(5), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(5)
    F = (1.d0 - sensor) * KEEP2(id,rho,p,V,Normal) &
        + sensor * SLAU(id_slau,id,rho,p,V,Normal,1.d0,sensor)
  end function flux_Weighted2

  attributes(device) function flux_Threshold6(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    real(8), intent(in), value                  :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(6), device   :: rho, p
    real(8), intent(in), dimension(6,3), device :: V
    real(8), intent(in), dimension(5), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(5)
    real(2) slau
    if (sensor < threshold) then
      F = KEEP6(id,rho,p,V,Normal)
    else
      F = flux_SLAU6(slau,id,rho,p,V,Normal,sensor)
    endif
  end function flux_Threshold6

  attributes(device) function flux_Threshold4(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    real(8), intent(in), value                  :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(4), device   :: rho, p
    real(8), intent(in), dimension(4,3), device :: V
    real(8), intent(in), dimension(5), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(5)
    real(2) slau
    if (sensor < threshold) then
      F = KEEP4(id,rho,p,V,Normal)
    else
      F = flux_SLAU4(slau,id,rho,p,V,Normal,sensor)
    endif
  end function flux_Threshold4

  attributes(device) function flux_Threshold2(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    use mod_globals, only : id_slau
    real(8), intent(in), value                  :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(2), device   :: rho, p
    real(8), intent(in), dimension(2,3), device :: V
    real(8), intent(in), dimension(5), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(5)
    if (sensor < threshold) then
      F = KEEP2(id,rho,p,V,Normal)
    else
      F = SLAU(id_slau,id,rho,p,V,Normal,1.d0,sensor)
    endif
  end function flux_Threshold2

  attributes(global) subroutine calc_E(nx, ny, nz, rho, u, v, w, p, E)
    use mod_globals, only : id_accuracy, id_scheme
    integer, intent(in), value                       :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w, p
    real(8), intent(out), device :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    integer i, j, k
    integer(kind=2) id_slau_wall
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4, e4
    real(8), dimension(4,3) :: V4
    real(8), dimension(6)   :: rho6, p6
    real(8), dimension(6,3) :: V6
    real(8), device         :: Normal(5) = (/0.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(8) fdx
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    if (3 <= i .and. i <= nx-3 .and. 8 <= kind(id_accuracy)) then
      !e4(:)   = QJ(i-1:i+2,j,k,5) !   e / Jacobian
      !rho4(:) = QJ(i-1:i+2,j,k,1) ! rho / Jacobian
      rho6(:) = rho(i-2:i+3,j,k)
      V6(:,1) =   u(i-2:i+3,j,k)
      V6(:,2) =   v(i-2:i+3,j,k)
      V6(:,3) =   w(i-2:i+3,j,k)
      p6(:)   =   p(i-2:i+3,j,k)
      !fdx     = Albada(e4, rho4)
      E(i,j-offset,k-offset,:) = flux6(id_scheme,1,rho6,p6,V6,Normal,fdx)
    elseif (2 <= i .and. i <= nx-2 .and. 4 <= kind(id_accuracy)) then
      !e4(:)   = QJ(i-1:i+2,j,k,5) * Jacobian(j)
      rho4(:) = rho(i-1:i+2,j,k)
      V4(:,1) =   u(i-1:i+2,j,k)
      V4(:,2) =   v(i-1:i+2,j,k)
      V4(:,3) =   w(i-1:i+2,j,k)
      p4(:)   =   p(i-1:i+2,j,k)
      !fdx     = Albada(e4, rho4)
      E(i,j-offset,k-offset,:) = flux4(id_scheme,1,rho4,p4,V4,Normal,fdx)
    else
      rho2(:) = rho(i:i+1,j,k)
      V2(:,1) =   u(i:i+1,j,k)
      V2(:,2) =   v(i:i+1,j,k)
      V2(:,3) =   w(i:i+1,j,k)
      p2(:)   =   p(i:i+1,j,k)
      !fdx     = 0.d0
      E(i,j-offset,k-offset,:) = flux2(id_scheme,1,rho2,p2,V2,Normal,fdx)
    endif
  end subroutine calc_E

  attributes(global) subroutine calc_F(nx, ny, nz, rho, u, v, w, p, F)
    use mod_globals, only : id_accuracy, id_scheme, slau_wall
    integer, intent(in), value                       :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w, p
    real(8), intent(out), device :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    integer i, j, k
    integer(kind=2) id_slau_wall
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4, e4
    real(8), dimension(4,3) :: V4
    real(8), dimension(6)   :: rho6, p6
    real(8), dimension(6,3) :: V6
    real(8), device         :: Normal(5) = (/0.d0, 0.d0, 1.d0, 0.d0, 0.d0/)
    real(8) fdy
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    if (3 <= j .and. j <= ny-3 .and. 8 <= kind(id_accuracy)) then
      !e4(:)   = QJ(i,j-1:j+2,k,5) !   e / Jacobian
      !rho4(:) = QJ(i,j-1:j+2,k,1) ! rho / Jacobian
      rho6(:) = rho(i,j-2:j+3,k)
      V6(:,1) =   u(i,j-2:j+3,k)
      V6(:,2) =   v(i,j-2:j+3,k)
      V6(:,3) =   w(i,j-2:j+3,k)
      p6(:)   =   p(i,j-2:j+3,k)
      !fdy     = Albada(e4, rho4)
      F(i-offset,j,k-offset,:) = flux6(id_scheme,2,rho6,p6,V6,Normal,fdy)
    elseif (2 <= j .and. j <= ny-2 .and. 4 <= kind(id_accuracy)) then
      !e4(:)   = QJ(i,j-1:j+2,k,5) * Jacobian(j-1:j+2)
      rho4(:) = rho(i,j-1:j+2,k)
      V4(:,1) =   u(i,j-1:j+2,k)
      V4(:,2) =   v(i,j-1:j+2,k)
      V4(:,3) =   w(i,j-1:j+2,k)
      p4(:)   =   p(i,j-1:j+2,k)
      !fdy     = Albada(e4, rho4)
      F(i-offset,j,k-offset,:) = flux4(id_scheme,2,rho4,p4,V4,Normal,fdy)
    else
      rho2(:) = rho(i,j:j+1,k)
      V2(:,1) =   u(i,j:j+1,k)
      V2(:,2) =   v(i,j:j+1,k)
      V2(:,3) =   w(i,j:j+1,k)
      p2(:)   =   p(i,j:j+1,k)
      !fdy     = 0.d0
      if (kind(slau_wall) /= 4) then
        F(i-offset,j,k-offset,:) = flux2(id_scheme,2,rho2,p2,V2,Normal,fdy)
      else
        F(i-offset,j,k-offset,:) = SLAU(id_slau_wall,2,rho2,p2,V2,Normal)
      endif
    endif
  end subroutine calc_F

  attributes(global) subroutine calc_G(nx, ny, nz, rho, u, v, w, p, G)
    use mod_globals, only : id_accuracy, id_scheme
    integer, intent(in), value                       :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device :: rho, u, v, w, p
    real(8), intent(out), device :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    integer i, j, k
    integer(kind=2) id_slau_wall
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4, e4
    real(8), dimension(4,3) :: V4
    real(8), dimension(6)   :: rho6, p6
    real(8), dimension(6,3) :: V6
    real(8), device         :: Normal(5) = (/0.d0, 0.d0, 0.d0, 1.d0, 0.d0/)
    real(8) :: fdz
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (3 <= k .and. k <= nz-3 .and. 8 <= kind(id_accuracy)) then
      !e4(:)   = QJ(i,j,k-1:k+2,5) !   e / Jacobian
      !rho4(:) = QJ(i,j,k-1:k+2,1) ! rho / Jacobian
      rho6(:) = rho(i,j,k-2:k+3)
      V6(:,1) =   u(i,j,k-2:k+3)
      V6(:,2) =   v(i,j,k-2:k+3)
      V6(:,3) =   w(i,j,k-2:k+3)
      p6(:)   =   p(i,j,k-2:k+3)
      !fdz     = Albada(e4, rho4)
      G(i-offset,j-offset,k,:) = flux6(id_scheme,3,rho6,p6,V6,Normal,fdz)
    elseif (2 <= k .and. k <= nz-2 .and. 4 <= kind(id_accuracy)) then
      !e4(:)   = QJ(i,j,k-1:k+2,5) * Jacobian(j)
      rho4(:) = rho(i,j,k-1:k+2)
      V4(:,1) =   u(i,j,k-1:k+2)
      V4(:,2) =   v(i,j,k-1:k+2)
      V4(:,3) =   w(i,j,k-1:k+2)
      p4(:)   =   p(i,j,k-1:k+2)
      !fdz     = Albada(e4, rho4)
      G(i-offset,j-offset,k,:) = flux4(id_scheme,3,rho4,p4,V4,Normal,fdz)
    else
      rho2(:) = rho(i,j,k:k+1)
      V2(:,1) =   u(i,j,k:k+1)
      V2(:,2) =   v(i,j,k:k+1)
      V2(:,3) =   w(i,j,k:k+1)
      p2(:)   =   p(i,j,k:k+1)
      !fdz     = 0.d0
      G(i-offset,j-offset,k,:) = flux2(id_scheme,3,rho2,p2,V2,Normal,fdz)
    endif
  end subroutine calc_G
end module calc_flux

