module calc_flux
  use mod_globals, only : id_scheme, id_sensor, id_muscl, id_igr, gamma, threshold, threadsE, threadsF
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

  interface calc_E
    module procedure calc_E2, calc_E4, calc_E6
  end interface calc_E

  interface calc_F
    module procedure calc_F2, calc_F4, calc_F6
  end interface calc_F
contains
  attributes(device) subroutine flux_KEEP6(id_scheme, id, rho, u, v, uu, p, T, Normal, sensor, F)
    integer(kind=2), intent(in), value :: id_scheme
    integer, intent(in), value         :: id
    real(8), intent(in), contiguous    :: rho(6), u(6), v(6), uu(6), p(6), T(6)
    real(8), intent(in), contiguous    :: Normal(4)
    real(8), intent(in), value         :: sensor
    real(8), intent(out), contiguous   :: F(4)
    F = KEEP6(rho, u, v, uu, p, T, Normal)
  end subroutine flux_KEEP6


  attributes(device) subroutine flux_KEEP4(id_scheme, id, rho, u, v, uu, p, T, Normal, sensor, F)
    integer(kind=2), intent(in), value :: id_scheme
    integer, intent(in), value         :: id
    real(8), intent(in), contiguous    :: rho(4), u(4), v(4), uu(4), p(4), T(4)
    real(8), intent(in), contiguous    :: Normal(4)
    real(8), intent(in), value         :: sensor
    real(8), intent(out), contiguous   :: F(4)
    F = KEEP4(rho, u, v, uu, p, T, Normal)
  end subroutine flux_KEEP4


  attributes(device) subroutine flux_KEEP2(id_scheme, id, rho, u, v, uu, p, T, Normal, sensor, F)
    integer(kind=2), intent(in), value :: id_scheme
    integer, intent(in), value         :: id
    real(8), intent(in), contiguous    :: rho(2), u(2), v(2), uu(2), p(2), T(2)
    real(8), intent(in), contiguous    :: Normal(4)
    real(8), intent(in), value         :: sensor
    real(8), intent(out), contiguous   :: F(4)
    F = KEEP2(rho, u, v, uu, p, T, Normal)
  end subroutine flux_KEEP2


  attributes(device) subroutine flux_IGR6(id, rho, u, v, uu, p, T, sigma, Normal, sensor, F)
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(6), u(6), v(6), uu(6), p(6), T(6), sigma(6)
    real(8), intent(in), contiguous  :: Normal(4)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(4)
    F = KEEP_IGR6(rho, u, v, uu, p, T, sigma, Normal)
  end subroutine flux_IGR6


  attributes(device) subroutine flux_IGR4(id, rho, u, v, uu, p, T, sigma, Normal, sensor, F)
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(4), u(4), v(4), uu(4), p(4), T(4), sigma(4)
    real(8), intent(in), contiguous  :: Normal(4)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(4)
    F = KEEP_IGR4(rho, u, v, uu, p, T, sigma, Normal)
  end subroutine flux_IGR4


  attributes(device) subroutine flux_IGR2(id, rho, u, v, uu, p, T, sigma, Normal, sensor, F)
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(2), u(2), v(2), uu(2), p(2), T(2), sigma(2)
    real(8), intent(in), contiguous  :: Normal(4)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(4)
    F = KEEP_IGR2(rho, u, v, uu, p, T, sigma, Normal)
  end subroutine flux_IGR2


  attributes(device) subroutine flux_SLAU6(id_scheme, id, rho, u, v, uu, p, T, Normal, sensor, F)
    use mod_globals, only : id_slau
    real(kind=2), intent(in), value  :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(6), u(6), v(6), uu(6), p(6), T(6)
    real(8), intent(in), contiguous  :: Normal(4)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(4)
    real(8) wiggle, rho2(2), p2(2), V2(2,2)
    wiggle = wiggle_detector(p(2:5))
    call calc_6points(sensor, rho, u, v, p, rho2, p2, V2)
    F = SLAU(id_slau, id, rho2, p2, V2, Normal, wiggle)
  end subroutine flux_SLAU6


  attributes(device) subroutine flux_SLAU4(id_scheme, id, rho, u, v, uu, p, T, Normal, sensor, F)
    use mod_globals, only : id_slau
    real(kind=2), intent(in), value  :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(4), u(4), v(4), uu(4), p(4), T(4)
    real(8), intent(in), contiguous  :: Normal(4)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(4)
    real(8) wiggle, rho2(2), p2(2), V2(2,2)
    wiggle = wiggle_detector(p)
    call calc_4points(sensor, rho, u, v, p, rho2, p2, V2)
    F = SLAU(id_slau, id, rho2, p2, V2, Normal, wiggle)
  end subroutine flux_SLAU4


  attributes(device) subroutine flux_SLAU2(id_scheme, id, rho, u, v, uu, p, T, Normal, sensor, F)
    use mod_globals, only : id_slau
    real(kind=2), intent(in), value  :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(2), u(2), v(2), uu(2), p(2), T(2)
    real(8), intent(in), contiguous  :: Normal(4)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(4)
    real(8) V2(2,2)
    V2(:,1) = u
    V2(:,2) = v
    F = SLAU(id_slau, id, rho, p, V2, Normal, 1.d0)
  end subroutine flux_SLAU2


  attributes(device) subroutine flux_Weighted6(id_scheme, id, rho, u, v, uu, p, T, Normal, sensor, F)
    real(4), intent(in), value       :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(6), u(6), v(6), uu(6), p(6), T(6)
    real(8), intent(in), contiguous  :: Normal(4)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(4)
    real(2) slau
    call flux_SLAU6(slau, id, rho, u, v, uu, p, T, Normal, sensor, F)
    F = sensor * F + (1.d0 - sensor) * KEEP6(rho, u, v, uu, p, T, Normal) 
  end subroutine flux_Weighted6


  attributes(device) subroutine flux_Weighted4(id_scheme, id, rho, u, v, uu, p, T, Normal, sensor, F)
    real(4), intent(in), value       :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(4), u(4), v(4), uu(4), p(4), T(4)
    real(8), intent(in), contiguous  :: Normal(4)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(4)
    real(2) slau
    call flux_SLAU4(slau, id, rho, u, v, uu, p, T, Normal, sensor, F)
    F = sensor * F + (1.d0 - sensor) * KEEP4(rho, u, v, uu, p, T, Normal)
  end subroutine flux_Weighted4


  attributes(device) subroutine flux_Weighted2(id_scheme, id, rho, u, v, uu, p, T, Normal, sensor, F)
    use mod_globals, only : id_slau
    real(4), intent(in), value       :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(2), u(2), v(2), uu(2), p(2), T(2)
    real(8), intent(in), contiguous  :: Normal(4)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(4)
    real(8) V2(2,2)
    V2(:,1) = u
    V2(:,2) = v
    F = (1.d0 - sensor) * KEEP2(rho, u, v, uu, p, T, Normal)
    F = F + sensor * SLAU(id_slau, id, rho, p, V2, Normal, 1.d0)
  end subroutine flux_Weighted2


  attributes(device) subroutine flux_Threshold6(id_scheme, id, rho, u, v, uu, p, T, Normal, sensor, F)
    real(8), intent(in), value       :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(6), u(6), v(6), uu(6), p(6), T(6)
    real(8), intent(in), contiguous  :: Normal(4)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(4)
    real(2) slau
    if (sensor < threshold) then
      F = KEEP6(rho, u, v, uu, p, T, Normal)
    else
      call flux_SLAU6(slau, id, rho, u, v, uu, p, T, Normal, sensor, F)
    endif
  end subroutine flux_Threshold6


  attributes(device) subroutine flux_Threshold4(id_scheme, id, rho, u, v, uu, p, T, Normal, sensor, F)
    real(8), intent(in), value       :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(4), u(4), v(4), uu(4), p(4), T(4)
    real(8), intent(in), contiguous  :: Normal(4)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(4)
    real(2) slau
    if (sensor < threshold) then
      F = KEEP4(rho, u, v, uu, p, T, Normal)
    else
      call flux_SLAU4(slau, id, rho, u, v, uu, p, T, Normal, sensor, F)
    endif
  end subroutine flux_Threshold4


  attributes(device) subroutine flux_Threshold2(id_scheme, id, rho, u, v, uu, p, T, Normal, sensor, F)
    use mod_globals, only : id_slau
    real(8), intent(in), value       :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(2), u(2), v(2), uu(2), p(2), T(2)
    real(8), intent(in), contiguous  :: Normal(4)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(4)
    real(8) V2(2,2)
    V2(:,1) = u
    V2(:,2) = v
    if (sensor < threshold) then
      F = KEEP2(rho, u, v, uu, p, T, Normal)
    else
      F = SLAU(id_slau, id, rho, p, V2, Normal, 1.d0)
    endif
  end subroutine flux_Threshold2


  attributes(global) subroutine calc_E6(id_accuracy, nx, ny, Q, T, sensor, E, sigma)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_x
    integer(kind=8), intent(in), value              :: id_accuracy
    integer, intent(in), value                      :: nx, ny
    real(8), intent(in), dimension(4,nx,ny), device :: Q
    real(8), intent(in), dimension(nx,ny), device   :: T, sensor
    real(8), intent(out), device                    :: E(4,nx-1,ny-2)
    real(4), intent(in), device, optional           :: sigma(nx,ny)
    integer i, j, it, jt, ii, i_base
    real(8), dimension(-1:threadsE%x+3,threadsE%y), shared :: rho, u, v, p
    real(8) fdx, tmp(6)
    it = threadIdx%x
    jt = threadIdx%y
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    do ii = it-2, threadsE%x+3, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
        rho(ii,jt) = Q(1,i,j)
          u(ii,jt) = Q(2,i,j)
          v(ii,jt) = Q(3,i,j)
          p(ii,jt) = Q(4,i,j)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j) return
    fdx = 0.5d0 * (sensor(i,j) + sensor(i+1,j))
    associate(uu => u)
    if (3 <= i .and. i <= nx-3) then
      tmp(:) = T(i-2:i+3,j)
      if (kind(id_igr) == 4) then
        block
          real(8), device :: sig(6)
          sig = dble(sigma(i-2:i+3,j))
          call flux_IGR6(1,rho(it-2:it+3,jt),u(it-2:it+3,jt),v(it-2:it+3,jt),&
                         uu(it-2:it+3,jt),p(it-2:it+3,jt),tmp,sig,Normal_x(1:4),fdx,E(:,i,j-1))
        end block
      else
        call flux6(id_scheme,1,rho(it-2:it+3,jt),u(it-2:it+3,jt),v(it-2:it+3,jt),&
                   uu(it-2:it+3,jt),p(it-2:it+3,jt),tmp,Normal_x(1:4),fdx,E(:,i,j-1))
      endif
    elseif (2 <= i .and. i <= nx-2) then
      tmp(2:5) = T(i-1:i+2,j)
      if (kind(id_igr) == 4) then
        block
          real(8), device :: sig(4)
          sig = dble(sigma(i-1:i+2,j))
          call flux_IGR4(1,rho(it-1:it+2,jt),u(it-1:it+2,jt),v(it-1:it+2,jt),&
                         uu(it-1:it+2,jt),p(it-1:it+2,jt),tmp(2:5),sig,Normal_x(1:4),fdx,E(:,i,j-1))
        end block
      else
        call flux4(id_scheme,1,rho(it-1:it+2,jt),u(it-1:it+2,jt),v(it-1:it+2,jt),&
                   uu(it-1:it+2,jt),p(it-1:it+2,jt),tmp(2:5),Normal_x(1:4),fdx,E(:,i,j-1))
      endif
    else
      tmp(3:4) = T(i:i+1,j)
      if (kind(id_igr) == 4) then
        block
          real(8), device :: sig(2)
          sig = dble(sigma(i:i+1,j))
          call flux_IGR2(1,rho(it:it+1,jt),u(it:it+1,jt),v(it:it+1,jt),&
                         uu(it:it+1,jt),p(it:it+1,jt),tmp(3:4),sig,Normal_x(1:4),fdx,E(:,i,j-1))
        end block
      else
        call flux2(id_scheme,1,rho(it:it+1,jt),u(it:it+1,jt),v(it:it+1,jt),&
                   uu(it:it+1,jt),p(it:it+1,jt),tmp(3:4),Normal_x(1:4),fdx,E(:,i,j-1))
      endif
    endif
    end associate
  end subroutine calc_E6


  attributes(global) subroutine calc_F6(id_accuracy, nx, ny, Q, T, sensor, F, sigma)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_y
    integer(kind=8), intent(in), value              :: id_accuracy
    integer, intent(in), value                      :: nx, ny
    real(8), intent(in), dimension(4,nx,ny), device :: Q
    real(8), intent(in), dimension(nx,ny), device   :: T, sensor
    real(8), intent(out), device                    :: F(4,nx-2,ny-1)
    real(4), intent(in), device, optional           :: sigma(nx,ny)
    integer i, j, it, jt, jj, j_base
    real(8), dimension(-1:threadsF%y+3,threadsF%x), shared :: rho, u, v, p
    real(8) fdy, tmp(6)
    it = threadIdx%x
    jt = threadIdx%y
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j_base = (blockIdx%y-1)*blockDim%y
    do jj = jt-2, threadsF%y+3, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
        rho(jj,it) = Q(1,i,j)
          u(jj,it) = Q(2,i,j)
          v(jj,it) = Q(3,i,j)
          p(jj,it) = Q(4,i,j)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j) return
    fdy = 0.5d0 * (sensor(i,j) + sensor(i,j+1))
    associate(vv => v)
    if (3 <= j .and. j <= ny-3) then
      tmp(:) = T(i,j-2:j+3)
      if (kind(id_igr) == 4) then
        block
          real(8), device :: sig(6)
          sig = dble(sigma(i,j-2:j+3))
          call flux_IGR6(2,rho(jt-2:jt+3,it),u(jt-2:jt+3,it),v(jt-2:jt+3,it),&
                         vv(jt-2:jt+3,it),p(jt-2:jt+3,it),tmp,sig,Normal_y(1:4),fdy,F(:,i-1,j))
        end block
      else
        call flux6(id_scheme,2,rho(jt-2:jt+3,it),u(jt-2:jt+3,it),v(jt-2:jt+3,it),&
                   vv(jt-2:jt+3,it),p(jt-2:jt+3,it),tmp,Normal_y(1:4),fdy,F(:,i-1,j))
      endif
    elseif (2 <= j .and. j <= ny-2) then
      tmp(2:5) = T(i,j-1:j+2)
      if (kind(id_igr) == 4) then
        block
          real(8), device :: sig(4)
          sig = dble(sigma(i,j-1:j+2))
          call flux_IGR4(2,rho(jt-1:jt+2,it),u(jt-1:jt+2,it),v(jt-1:jt+2,it),&
                         vv(jt-1:jt+2,it),p(jt-1:jt+2,it),tmp(2:5),sig,Normal_y(1:4),fdy,F(:,i-1,j))
        end block
      else
        call flux4(id_scheme,2,rho(jt-1:jt+2,it),u(jt-1:jt+2,it),v(jt-1:jt+2,it),&
                   vv(jt-1:jt+2,it),p(jt-1:jt+2,it),tmp(2:5),Normal_y(1:4),fdy,F(:,i-1,j))
      endif
    else
      tmp(3:4) = T(i,j:j+1)
      if (kind(id_igr) == 4) then
        block
          real(8), device :: sig(2)
          sig = dble(sigma(i,j:j+1))
          call flux_IGR2(2,rho(jt:jt+1,it),u(jt:jt+1,it),v(jt:jt+1,it),&
                         vv(jt:jt+1,it),p(jt:jt+1,it),tmp(3:4),sig,Normal_y(1:4),fdy,F(:,i-1,j))
        end block
      else
        call flux2(id_scheme,2,rho(jt:jt+1,it),u(jt:jt+1,it),v(jt:jt+1,it),&
                   vv(jt:jt+1,it),p(jt:jt+1,it),tmp(3:4),Normal_y(1:4),fdy,F(:,i-1,j))
      endif
    endif
    end associate
  end subroutine calc_F6


  attributes(global) subroutine calc_E4(id_accuracy, nx, ny, Q, T, sensor, E, sigma)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_x
    integer(kind=4), intent(in), value              :: id_accuracy
    integer, intent(in), value                      :: nx, ny
    real(8), intent(in), dimension(4,nx,ny), device :: Q
    real(8), intent(in), dimension(nx,ny), device   :: T, sensor
    real(8), intent(out), device                    :: E(4,nx-1,ny-2)
    real(4), intent(in), device, optional           :: sigma(nx,ny)
    integer i, j, it, jt, ii, i_base
    real(8), dimension(0:threadsE%x+2,threadsE%y), shared :: rho, u, v, p
    real(8) fdx, tmp(4)
    it = threadIdx%x
    jt = threadIdx%y
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    do ii = it-1, threadsE%x+2, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
        rho(ii,jt) = Q(1,i,j)
          u(ii,jt) = Q(2,i,j)
          v(ii,jt) = Q(3,i,j)
          p(ii,jt) = Q(4,i,j)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j) return
    fdx = 0.5d0 * (sensor(i,j) + sensor(i+1,j))
    associate(uu => u)
    if (2 <= i .and. i <= nx-2) then
      tmp(:) = T(i-1:i+2,j)
      if (kind(id_igr) == 4) then
        block
          real(8), device :: sig(4)
          sig = dble(sigma(i-1:i+2,j))
          call flux_IGR4(1,rho(it-1:it+2,jt),u(it-1:it+2,jt),v(it-1:it+2,jt),&
                         uu(it-1:it+2,jt),p(it-1:it+2,jt),tmp,sig,Normal_x(1:4),fdx,E(:,i,j-1))
        end block
      else
        call flux4(id_scheme,1,rho(it-1:it+2,jt),u(it-1:it+2,jt),v(it-1:it+2,jt),&
                   uu(it-1:it+2,jt),p(it-1:it+2,jt),tmp,Normal_x(1:4),fdx,E(:,i,j-1))
      endif
    else
      tmp(2:3) = T(i:i+1,j)
      if (kind(id_igr) == 4) then
        block
          real(8), device :: sig(2)
          sig = dble(sigma(i:i+1,j))
          call flux_IGR2(1,rho(it:it+1,jt),u(it:it+1,jt),v(it:it+1,jt),&
                         uu(it:it+1,jt),p(it:it+1,jt),tmp(2:3),sig,Normal_x(1:4),fdx,E(:,i,j-1))
        end block
      else
        call flux2(id_scheme,1,rho(it:it+1,jt),u(it:it+1,jt),v(it:it+1,jt),&
                   uu(it:it+1,jt),p(it:it+1,jt),tmp(2:3),Normal_x(1:4),fdx,E(:,i,j-1))
      endif
    endif
    end associate
  end subroutine calc_E4


  attributes(global) subroutine calc_F4(id_accuracy, nx, ny, Q, T, sensor, F, sigma)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_y
    integer(kind=4), intent(in), value              :: id_accuracy
    integer, intent(in), value                      :: nx, ny
    real(8), intent(in), dimension(4,nx,ny), device :: Q
    real(8), intent(in), dimension(nx,ny), device   :: T, sensor
    real(8), intent(out), device                    :: F(4,nx-2,ny-1)
    real(4), intent(in), device, optional           :: sigma(nx,ny)
    integer i, j, it, jt, jj, j_base
    real(8), dimension(0:threadsF%y+2,threadsF%x), shared :: rho, u, v, p
    real(8) fdy, tmp(4)
    it = threadIdx%x
    jt = threadIdx%y
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j_base = (blockIdx%y-1)*blockDim%y
    do jj = jt-1, threadsF%y+2, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny) then
        rho(jj,it) = Q(1,i,j)
          u(jj,it) = Q(2,i,j)
          v(jj,it) = Q(3,i,j)
          p(jj,it) = Q(4,i,j)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j) return
    fdy = 0.5d0 * (sensor(i,j) + sensor(i,j+1))
    associate(vv => v)
    if (2 <= j .and. j <= ny-2) then
      tmp(:) = T(i,j-1:j+2)
      if (kind(id_igr) == 4) then
        block
          real(8), device :: sig(4)
          sig = dble(sigma(i,j-1:j+2))
          call flux_IGR4(2,rho(jt-1:jt+2,it),u(jt-1:jt+2,it),v(jt-1:jt+2,it),&
                         vv(jt-1:jt+2,it),p(jt-1:jt+2,it),tmp,sig,Normal_y(1:4),fdy,F(:,i-1,j))
        end block
      else
        call flux4(id_scheme,2,rho(jt-1:jt+2,it),u(jt-1:jt+2,it),v(jt-1:jt+2,it),&
                   vv(jt-1:jt+2,it),p(jt-1:jt+2,it),tmp,Normal_y(1:4),fdy,F(:,i-1,j))
      endif
    else
      tmp(2:3) = T(i,j:j+1)
      if (kind(id_igr) == 4) then
        block
          real(8), device :: sig(2)
          sig = dble(sigma(i,j:j+1))
          call flux_IGR2(2,rho(jt:jt+1,it),u(jt:jt+1,it),v(jt:jt+1,it),&
                         vv(jt:jt+1,it),p(jt:jt+1,it),tmp(2:3),sig,Normal_y(1:4),fdy,F(:,i-1,j))
        end block
      else
        call flux2(id_scheme,2,rho(jt:jt+1,it),u(jt:jt+1,it),v(jt:jt+1,it),&
                   vv(jt:jt+1,it),p(jt:jt+1,it),tmp(2:3),Normal_y(1:4),fdy,F(:,i-1,j))
      endif
    endif
    end associate
  end subroutine calc_F4


  attributes(global) subroutine calc_E2(id_accuracy, nx, ny, Q, T, sensor, E, sigma)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_x
    integer(kind=2), intent(in), value              :: id_accuracy
    integer, intent(in), value                      :: nx, ny
    real(8), intent(in), dimension(4,nx,ny), device :: Q
    real(8), intent(in), dimension(nx,ny), device   :: T, sensor
    real(8), intent(out), device                    :: E(4,nx-1,ny-2)
    real(4), intent(in), device, optional           :: sigma(nx,ny)
    integer i, j, it, jt
    real(8), dimension(2), device :: rho, u, v, p, tmp
    real(8) fdx
    it = threadIdx%x
    jt = threadIdx%y
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    if (nx-1 < i .or. ny-1 < j) return
    fdx = 0.5d0 * (sensor(i,j) + sensor(i+1,j))
    rho = Q(1,i:i+1,j)
    u   = Q(2,i:i+1,j)
    v   = Q(3,i:i+1,j)
    p   = Q(4,i:i+1,j)
    tmp = T(i:i+1,j)
    if (kind(id_igr) == 4) then
      block
        real(8), device :: sig(2)
        sig = dble(sigma(i:i+1,j))
        sig = sig * sensor(i:i+1,j)
        call flux_IGR2(1, rho, u, v, u, p, tmp, sig, Normal_x(1:4), fdx, E(:,i,j-1))
      end block
    else
      call flux2(id_scheme, 1, rho, u, v, u, p, tmp, Normal_x(1:4), fdx, E(:,i,j-1))
    endif
  end subroutine calc_E2


  attributes(global) subroutine calc_F2(id_accuracy, nx, ny, Q, T, sensor, F, sigma)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_y
    integer(kind=2), intent(in), value              :: id_accuracy
    integer, intent(in), value                      :: nx, ny
    real(8), intent(in), dimension(4,nx,ny), device :: Q
    real(8), intent(in), dimension(nx,ny), device   :: T, sensor
    real(8), intent(out), device                    :: F(4,nx-2,ny-1)
    real(4), intent(in), device, optional           :: sigma(nx,ny)
    integer i, j, it, jt
    real(8), dimension(2), device :: rho, u, v, p, tmp
    real(8) fdy
    it = threadIdx%x
    jt = threadIdx%y
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j) return
    fdy = 0.5d0 * (sensor(i,j) + sensor(i,j+1))
    rho = Q(1,i,j:j+1)
    u   = Q(2,i,j:j+1)
    v   = Q(3,i,j:j+1)
    p   = Q(4,i,j:j+1)
    tmp = T(i,j:j+1)
    if (kind(id_igr) == 4) then
      block
        real(8), device :: sig(2)
        sig = dble(sigma(i,j:j+1))
        sig = sig * sensor(i,j:j+1)
        call flux_IGR2(2, rho, u, v, v, p, tmp, sig, Normal_y(1:4), fdy, F(:,i-1,j))
      end block
    else
      call flux2(id_scheme, 2, rho, u, v, v, p, tmp, Normal_y(1:4), fdy, F(:,i-1,j))
    endif
  end subroutine calc_F2
end module calc_flux

