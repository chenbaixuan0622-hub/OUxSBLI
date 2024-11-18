module calc_flux
  use mod_globals, only : id_scheme, id_sensor, id_muscl, accuracy, offset, gamma, threshold
  use calc_keep
  use calc_slau
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
    real(8), intent(in), dimension(6,2), device :: V
    real(8), intent(in), dimension(4), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(4)
    F = KEEP6(id,rho,p,V,Normal)
  end function flux_KEEP6

  attributes(device) function flux_KEEP4(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    integer(kind=2), intent(in), value          :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(4), device   :: rho, p
    real(8), intent(in), dimension(4,2), device :: V
    real(8), intent(in), dimension(4), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(4)
    F = KEEP4(id,rho,p,V,Normal)
  end function flux_KEEP4
  
  attributes(device) function flux_KEEP2(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    integer(kind=2), intent(in), value          :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(2), device   :: rho, p
    real(8), intent(in), dimension(2,2), device :: V
    real(8), intent(in), dimension(4), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(4)
    F = KEEP2(id,rho,p,V,Normal)
  end function flux_KEEP2

  attributes(device) function flux_SLAU6(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    use mod_globals, only : id_slau
    real(kind=2), intent(in), value             :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(6), device   :: rho, p
    real(8), intent(in), dimension(6,2), device :: V
    real(8), intent(in), dimension(4), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) wiggle, rho2(2), p2(2), V2(2,2), p4(4), F(4)
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
    real(8), intent(in), dimension(4,2), device :: V
    real(8), intent(in), dimension(4), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) wiggle, rho2(2), p2(2), V2(2,2), F(4)
    wiggle = wiggle_detector(p)
    call calc_4points(1.d0,1.d0,1.d0/3.d0,sensor,rho,p,V,rho2,p2,V2)
    F = SLAU(id_slau,id,rho2,p2,V2,Normal,wiggle,sensor)
  end function flux_SLAU4

  attributes(device) function flux_SLAU2(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    use mod_globals, only : id_slau
    real(kind=2), intent(in), value             :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(2), device   :: rho, p
    real(8), intent(in), dimension(2,2), device :: V
    real(8), intent(in), dimension(4), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(4)
    F = SLAU(id_slau,id,rho,p,V,Normal,1.d0,sensor)
  end function flux_SLAU2

  attributes(device) function flux_Weighted6(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    real(4), intent(in), value                  :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(6), device   :: rho, p
    real(8), intent(in), dimension(6,2), device :: V
    real(8), intent(in), dimension(4), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(4)
    real(2) slau
    F = (1.d0 - sensor) * KEEP6(id,rho,p,V,Normal) &
        + sensor * flux_SLAU6(slau,id,rho,p,V,Normal,sensor)
  end function flux_Weighted6

  attributes(device) function flux_Weighted4(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    real(4), intent(in), value                  :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(4), device   :: rho, p
    real(8), intent(in), dimension(4,2), device :: V
    real(8), intent(in), dimension(4), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(4)
    real(2) slau
    F = (1.d0 - sensor) * KEEP4(id,rho,p,V,Normal) &
        + sensor * flux_SLAU4(slau,id,rho,p,V,Normal,sensor)
  end function flux_Weighted4

  attributes(device) function flux_Weighted2(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    use mod_globals, only : id_slau
    real(4), intent(in), value                 :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(2), device   :: rho, p
    real(8), intent(in), dimension(2,2), device :: V
    real(8), intent(in), dimension(4), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(4)
    F = (1.d0 - sensor) * KEEP2(id,rho,p,V,Normal) &
        + sensor * SLAU(id_slau,id,rho,p,V,Normal,1.d0,sensor)
  end function flux_Weighted2

  attributes(device) function flux_Threshold6(id_scheme,id,rho,p,V,Normal,sensor) result(F)
    real(8), intent(in), value                  :: id_scheme
    integer, intent(in), value                  :: id
    real(8), intent(in), dimension(6), device   :: rho, p
    real(8), intent(in), dimension(6,2), device :: V
    real(8), intent(in), dimension(4), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(4)
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
    real(8), intent(in), dimension(4,2), device :: V
    real(8), intent(in), dimension(4), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(4)
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
    real(8), intent(in), dimension(2,2), device :: V
    real(8), intent(in), dimension(4), device   :: Normal
    real(8), intent(in), value                  :: sensor
    real(8) F(4)
    if (sensor < threshold) then
      F = KEEP2(id,rho,p,V,Normal)
    else
      F = SLAU(id_slau,id,rho,p,V,Normal,1.d0,sensor)
    endif
  end function flux_Threshold2

  attributes(global) subroutine calc_E(nx, ny, rho, u, v, p, E)
    use mod_globals, only : id_accuracy, id_scheme
    integer, intent(in), value                    :: nx, ny
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(out), device :: E(nx-accuracy+1,ny-accuracy,4)
    integer i, j
    integer(kind=2) id_slau_wall
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,2) :: V2
    real(8), dimension(4)   :: rho4, p4, e4
    real(8), dimension(4,2) :: V4
    real(8), dimension(6)   :: rho6, p6
    real(8), dimension(6,2) :: V6
    real(8), device         :: Normal(4) = (/0.d0, 1.d0, 0.d0, 0.d0/)
    real(8) fdx
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    if (3 <= i .and. i <= nx-3 .and. 8 <= kind(id_accuracy)) then
      rho6(:) = rho(i-2:i+3,j)
      V6(:,1) =   u(i-2:i+3,j)
      V6(:,2) =   v(i-2:i+3,j)
      p6(:)   =   p(i-2:i+3,j)
      E(i,j-offset,:) = flux6(id_scheme,1,rho6,p6,V6,Normal,fdx)
    elseif (2 <= i .and. i <= nx-2 .and. 4 <= kind(id_accuracy)) then
      rho4(:) = rho(i-1:i+2,j)
      V4(:,1) =   u(i-1:i+2,j)
      V4(:,2) =   v(i-1:i+2,j)
      p4(:)   =   p(i-1:i+2,j)
      E(i,j-offset,:) = flux4(id_scheme,1,rho4,p4,V4,Normal,fdx)
    else
      rho2(:) = rho(i:i+1,j)
      V2(:,1) =   u(i:i+1,j)
      V2(:,2) =   v(i:i+1,j)
      p2(:)   =   p(i:i+1,j)
      E(i,j-offset,:) = flux2(id_scheme,1,rho2,p2,V2,Normal,fdx)
    endif
  end subroutine calc_E

  attributes(global) subroutine calc_F(nx, ny, rho, u, v, p, F)
    use mod_globals, only : id_accuracy, id_scheme, slau_wall
    integer, intent(in), value                    :: nx, ny
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(out), device :: F(nx-accuracy,ny-accuracy+1,4)
    integer i, j
    integer(kind=2) id_slau_wall
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,2) :: V2
    real(8), dimension(4)   :: rho4, p4, e4
    real(8), dimension(4,2) :: V4
    real(8), dimension(6)   :: rho6, p6
    real(8), dimension(6,2) :: V6
    real(8), device         :: Normal(4) = (/0.d0, 0.d0, 1.d0, 0.d0/)
    real(8) fdy
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    if (3 <= j .and. j <= ny-3 .and. 8 <= kind(id_accuracy)) then
      rho6(:) = rho(i,j-2:j+3)
      V6(:,1) =   u(i,j-2:j+3)
      V6(:,2) =   v(i,j-2:j+3)
      p6(:)   =   p(i,j-2:j+3)
      F(i-offset,j,:) = flux6(id_scheme,2,rho6,p6,V6,Normal,fdy)
    elseif (2 <= j .and. j <= ny-2 .and. 4 <= kind(id_accuracy)) then
      rho4(:) = rho(i,j-1:j+2)
      V4(:,1) =   u(i,j-1:j+2)
      V4(:,2) =   v(i,j-1:j+2)
      p4(:)   =   p(i,j-1:j+2)
      F(i-offset,j,:) = flux4(id_scheme,2,rho4,p4,V4,Normal,fdy)
    else
      rho2(:) = rho(i,j:j+1)
      V2(:,1) =   u(i,j:j+1)
      V2(:,2) =   v(i,j:j+1)
      p2(:)   =   p(i,j:j+1)
      if (kind(slau_wall) /= 4) then
        F(i-offset,j,:) = flux2(id_scheme,2,rho2,p2,V2,Normal,fdy)
      else
        F(i-offset,j,:) = SLAU(id_slau_wall,2,rho2,p2,V2,Normal)
      endif
    endif
  end subroutine calc_F
end module calc_flux

