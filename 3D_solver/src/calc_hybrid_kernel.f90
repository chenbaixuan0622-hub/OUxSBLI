module calc_hybrid_kernel
  use mod_globals, only : id_slau, gamma, threshold, threadsE, threadsF, threadsG
  use mod_constant, only : over_gamma_1, R_over_gamma_1, one_third, one_sixth, one_twelfth, two_third
  use calc_muscl
  use calc_hybrid
  implicit none
  private
  public calc_hybrid_x, calc_hybrid_y, calc_hybrid_z
  real(8), parameter :: one_24        = 1.d0 / 24.d0
  real(8), parameter :: one_48        = 1.d0 / 48.d0
  real(8), parameter :: one_60        = 1.d0 / 60.d0
  real(8), parameter :: one_120       = 1.d0 / 120.d0
  real(8), parameter :: one_240       = 1.d0 / 240.d0
  real(8), parameter :: seven_twelfth = 7.d0 / 12.d0

  interface KEEP
    module procedure KEEP2, KEEP4, KEEP6
  end interface KEEP
  
  interface SLAU
    module procedure SLAU1, HRSLAU2
  end interface SLAU

  interface calc_hybrid_x
    module procedure calc_hybrid_x2, calc_hybrid_x4, calc_hybrid_x6
  end interface calc_hybrid_x

  interface calc_hybrid_y
    module procedure calc_hybrid_y2, calc_hybrid_y4, calc_hybrid_y6
  end interface calc_hybrid_y

  interface calc_hybrid_z
    module procedure calc_hybrid_z2, calc_hybrid_z4, calc_hybrid_z6
  end interface calc_hybrid_z
contains
  include 'calc_keep_3d.f90'
  include 'calc_slau_3d.f90'

  attributes(global) subroutine calc_hybrid_x6(id_accuracy, nx, ny, nz, Q, T, sensor, E)
    use mod_constant, only : Normal_x
    integer(kind=8), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: E(nx-1,5,ny-2,nz-2)
    integer i, j, k, it, jt, kt, ii, i_base
    real(8), dimension(-1:threadsE%x+3,threadsE%y,threadsE%z), shared :: rho,  u,  v,  w,  p
    real(8), dimension(   threadsE%x,  threadsE%y,threadsE%z), shared :: rhor, ur, vr, wr, pr
    real(8) fdx
    integer(kind=4) id_accuracy4
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    do ii = it-2, threadsE%x+3, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(ii,jt,kt) = Q(i,1,j,k)
          u(ii,jt,kt) = Q(i,2,j,k)
          v(ii,jt,kt) = Q(i,3,j,k)
          w(ii,jt,kt) = Q(i,4,j,k)
          p(ii,jt,kt) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    block
      real(8) rhol, ul, vl, wl, pl
      fdx = 0.5d0 * (sensor(i,j,k) + sensor(i+1,j,k))
      if (3 <= i .and. i <= nx-3) then
        if (fdx <= threshold) then
          block
            real(8) tmp(6)
            tmp = T(i-2:i+3,j,k)
            associate(uu => u)
              E(i,:,j-1,k-1) = KEEP(id_accuracy, &
                                    rho(it-2:it+3,jt,kt), u(it-2:it+3,jt,kt), &
                                      v(it-2:it+3,jt,kt), w(it-2:it+3,jt,kt), &
                                     uu(it-2:it+3,jt,kt), p(it-2:it+3,jt,kt), tmp, Normal_x)
            end associate
          end block
        else
          call delta6(fdx, rho(it-2:it+3,jt,kt), rhol, rhor(it,jt,kt))
          call delta6(fdx,   u(it-2:it+3,jt,kt),   ul,   ur(it,jt,kt))
          call delta6(fdx,   v(it-2:it+3,jt,kt),   vl,   vr(it,jt,kt))
          call delta6(fdx,   w(it-2:it+3,jt,kt),   wl,   wr(it,jt,kt))
          call delta6(fdx,   p(it-2:it+3,jt,kt),   pl,   pr(it,jt,kt))
        endif
      elseif (2 <= i .and. i <= nx-2) then
        if (fdx <= threshold) then
          block
            real(8) tmp(4)
            tmp = T(i-1:i+2,j,k)
            associate(uu => u)
              E(i,:,j-1,k-1) = KEEP(id_accuracy4, &
                                    rho(it-1:it+2,jt,kt), u(it-1:it+2,jt,kt), &
                                      v(it-1:it+2,jt,kt), w(it-1:it+2,jt,kt), &
                                     uu(it-1:it+2,jt,kt), p(it-1:it+2,jt,kt), tmp, Normal_x)
            end associate
          end block
        else
          call delta4(fdx, rho(it-1:it+2,jt,kt), rhol, rhor(it,jt,kt))
          call delta4(fdx,   u(it-1:it+2,jt,kt),   ul,   ur(it,jt,kt))
          call delta4(fdx,   v(it-1:it+2,jt,kt),   vl,   vr(it,jt,kt))
          call delta4(fdx,   w(it-1:it+2,jt,kt),   wl,   wr(it,jt,kt))
          call delta4(fdx,   p(it-1:it+2,jt,kt),   pl,   pr(it,jt,kt))
        endif
      else
        if (fdx <= threshold) then
          block
            real(8) tmp(2)
            tmp = T(i:i+1,j,k)
            associate(uu => u)
              E(i,:,j-1,k-1) = KEEP(id_accuracy2, &
                                    rho(it:it+1,jt,kt), u(it:it+1,jt,kt), &
                                      v(it:it+1,jt,kt), w(it:it+1,jt,kt), &
                                     uu(it:it+1,jt,kt), p(it:it+1,jt,kt), tmp, Normal_x)
            end associate
          end block
        else
          rhol = rho(it,jt,kt); rhor(it,jt,kt) = rho(it+1,jt,kt)
            ul =   u(it,jt,kt);   ur(it,jt,kt) =   u(it+1,jt,kt)
            vl =   v(it,jt,kt);   vr(it,jt,kt) =   v(it+1,jt,kt)
            wl =   w(it,jt,kt);   wr(it,jt,kt) =   w(it+1,jt,kt)
            pl =   p(it,jt,kt);   pr(it,jt,kt) =   p(it+1,jt,kt)
        endif
      endif
      call syncthreads()
      rho(it,jt,kt) = rhol
        u(it,jt,kt) = ul
        v(it,jt,kt) = vl
        w(it,jt,kt) = wl
        p(it,jt,kt) = pl
    end block
    if (fdx > threshold) then
      associate(un1 => u(it,jt,kt), un2 => ur(it,jt,kt))
        call SLAU(id_slau, rho(it,jt,kt), rhor(it,jt,kt), u(it,jt,kt), ur(it,jt,kt), v(it,jt,kt), vr(it,jt,kt), &
                  w(it,jt,kt), wr(it,jt,kt), un1, un2, p(it,jt,kt), pr(it,jt,kt), Normal_x, 1.d0, &
                  E(i,1,j-1,k-1), E(i,2,j-1,k-1), E(i,3,j-1,k-1), E(i,4,j-1,k-1), E(i,5,j-1,k-1))
      end associate
    endif
  end subroutine calc_hybrid_x6


  attributes(global) subroutine calc_hybrid_y6(id_accuracy, nx, ny, nz, Q, T, sensor, F)
    use mod_constant, only : Normal_y
    integer(kind=8), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: F(nx-2,5,ny-1,nz-2)
    integer i, j, k, it, jt, kt, jj, j_base
    real(8), dimension(-1:threadsF%y+3,threadsF%x,threadsF%z), shared :: rho,  u,  v,  w,  p
    real(8), dimension(   threadsF%y,  threadsF%x,threadsF%z), shared :: rhor, ur, vr, wr, pr
    real(8) fdy
    integer(kind=4) id_accuracy4
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    do jj = jt-2, threadsF%y+3, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(jj,it,kt) = Q(i,1,j,k)
          u(jj,it,kt) = Q(i,2,j,k)
          v(jj,it,kt) = Q(i,3,j,k)
          w(jj,it,kt) = Q(i,4,j,k)
          p(jj,it,kt) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    block
      real(8) rhol, ul, vl, wl, pl
      fdy = 0.5d0 * (sensor(i,j,k) + sensor(i,j+1,k))
      if (3 <= j .and. j <= ny-3 .and. 8 <= kind(id_accuracy)) then
        if (fdy <= threshold) then
          block
            real(8) tmp(6)
            tmp = T(i,j-2:j+3,k)
            associate(vv => v)
              F(i-1,:,j,k-1) = KEEP(id_accuracy, &
                                    rho(jt-2:jt+3,it,kt), u(jt-2:jt+3,it,kt), &
                                      v(jt-2:jt+3,it,kt), w(jt-2:jt+3,it,kt), &
                                     vv(jt-2:jt+3,it,kt), p(jt-2:jt+3,it,kt), tmp, Normal_y)
            end associate
          end block
        else
          call delta6(fdy, rho(jt-2:jt+3,it,kt), rhol, rhor(jt,it,kt))
          call delta6(fdy,   u(jt-2:jt+3,it,kt),   ul,   ur(jt,it,kt))
          call delta6(fdy,   v(jt-2:jt+3,it,kt),   vl,   vr(jt,it,kt))
          call delta6(fdy,   w(jt-2:jt+3,it,kt),   wl,   wr(jt,it,kt))
          call delta6(fdy,   p(jt-2:jt+3,it,kt),   pl,   pr(jt,it,kt))
        endif
      elseif (2 <= j .and. j <= ny-2) then
        if (fdy <= threshold) then
          block
            real(8) tmp(4)
            tmp = T(i,j-1:j+2,k)
            associate(vv => v)
              F(i-1,:,j,k-1) = KEEP(id_accuracy4, &
                                    rho(jt-1:jt+2,it,kt), u(jt-1:jt+2,it,kt), &
                                      v(jt-1:jt+2,it,kt), w(jt-1:jt+2,it,kt), &
                                     vv(jt-1:jt+2,it,kt), p(jt-1:jt+2,it,kt), tmp, Normal_y)
            end associate
          end block
        else
          call delta4(fdy, rho(jt-1:jt+2,it,kt), rhol, rhor(jt,it,kt))
          call delta4(fdy,   u(jt-1:jt+2,it,kt),   ul,   ur(jt,it,kt))
          call delta4(fdy,   v(jt-1:jt+2,it,kt),   vl,   vr(jt,it,kt))
          call delta4(fdy,   w(jt-1:jt+2,it,kt),   wl,   wr(jt,it,kt))
          call delta4(fdy,   p(jt-1:jt+2,it,kt),   pl,   pr(jt,it,kt))
        endif
      else
        if (fdy <= threshold) then
          block
            real(8) tmp(2)
            tmp = T(i,j:j+1,k)
            associate(vv => v)
              F(i-1,:,j,k-1) = KEEP(id_accuracy2, &
                                    rho(jt:jt+1,it,kt), u(jt:jt+1,it,kt), &
                                      v(jt:jt+1,it,kt), w(jt:jt+1,it,kt), &
                                     vv(jt:jt+1,it,kt), p(jt:jt+1,it,kt), tmp, Normal_y)
            end associate
          end block
        else
          rhol = rho(jt,it,kt); rhor(jt,it,kt) = rho(jt+1,it,kt)
            ul =   u(jt,it,kt);   ur(jt,it,kt) =   u(jt+1,it,kt)
            vl =   v(jt,it,kt);   vr(jt,it,kt) =   v(jt+1,it,kt)
            wl =   w(jt,it,kt);   wr(jt,it,kt) =   w(jt+1,it,kt)
            pl =   p(jt,it,kt);   pr(jt,it,kt) =   p(jt+1,it,kt)
        endif
      endif
      call syncthreads()
      rho(jt,it,kt) = rhol
        u(jt,it,kt) = ul
        v(jt,it,kt) = vl
        w(jt,it,kt) = wl
        p(jt,it,kt) = pl
    end block
    if (fdy > threshold) then
      associate(un1 => v(jt,it,kt), un2 => vr(jt,it,kt))
        call SLAU(id_slau, rho(jt,it,kt), rhor(jt,it,kt), u(jt,it,kt), ur(jt,it,kt), v(jt,it,kt), vr(jt,it,kt), &
                  w(jt,it,kt), wr(jt,it,kt), un1, un2, p(jt,it,kt), pr(jt,it,kt), Normal_y, 1.d0, &
                  F(i-1,1,j,k-1), F(i-1,2,j,k-1), F(i-1,3,j,k-1), F(i-1,4,j,k-1), F(i-1,5,j,k-1))
      end associate
    endif
  end subroutine calc_hybrid_y6


  attributes(global) subroutine calc_hybrid_z6(id_accuracy, nx, ny, nz, Q, T, sensor, G)
    use mod_constant, only : Normal_z
    integer(kind=8), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: G(nx-2,5,ny-2,nz-1)
    integer i, j, k, it, jt, kt, kk, k_base
    real(8), dimension(-1:threadsG%z+3,threadsG%y,threadsG%x), shared :: rho,  u,  v,  w,  p
    real(8), dimension(   threadsG%z,  threadsG%y,threadsG%x), shared :: rhor, ur, vr, wr, pr
    real(8) fdz
    integer(kind=4) id_accuracy4
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    do kk = kt-2, threadsG%z+3, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(kk,jt,it) = Q(i,1,j,k)
          u(kk,jt,it) = Q(i,2,j,k)
          v(kk,jt,it) = Q(i,3,j,k)
          w(kk,jt,it) = Q(i,4,j,k)
          p(kk,jt,it) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
      block
      real(8) rhol, ul, vl, wl, pl
      fdz = 0.5d0 * (sensor(i,j,k) + sensor(i,j,k+1))
      if (3 <= k .and. k <= nz-3 .and. 8 <= kind(id_accuracy)) then
        if (fdz <= threshold) then
          block
            real(8) tmp(6)
            tmp = T(i,j,k-2:k+3)
            associate(ww => w)
              G(i-1,:,j-1,k) = KEEP(id_accuracy, &
                                    rho(kt-2:kt+3,jt,it), u(kt-2:kt+3,jt,it), &
                                      v(kt-2:kt+3,jt,it), w(kt-2:kt+3,jt,it), &
                                     ww(kt-2:kt+3,jt,it), p(kt-2:kt+3,jt,it), tmp, Normal_z)
            end associate
          end block
        else
          call delta6(fdz, rho(kt-2:kt+3,jt,it), rhol, rhor(kt,jt,it))
          call delta6(fdz,   u(kt-2:kt+3,jt,it),   ul,   ur(kt,jt,it))
          call delta6(fdz,   v(kt-2:kt+3,jt,it),   vl,   vr(kt,jt,it))
          call delta6(fdz,   w(kt-2:kt+3,jt,it),   wl,   wr(kt,jt,it))
          call delta6(fdz,   p(kt-2:kt+3,jt,it),   pl,   pr(kt,jt,it))
        endif
      elseif (2 <= k .and. k <= nz-2) then
        if (fdz <= threshold) then
          block
            real(8) tmp(4)
            tmp = T(i,j,k-1:k+2)
            associate(ww => w)
              G(i-1,:,j-1,k) = KEEP(id_accuracy4, &
                                    rho(kt-1:kt+2,jt,it), u(kt-1:kt+2,jt,it), &
                                      v(kt-1:kt+2,jt,it), w(kt-1:kt+2,jt,it), &
                                     ww(kt-1:kt+2,jt,it), p(kt-1:kt+2,jt,it), tmp, Normal_z)
            end associate
          end block
        else
          call delta4(fdz, rho(kt-1:kt+2,jt,it), rhol, rhor(kt,jt,it))
          call delta4(fdz,   u(kt-1:kt+2,jt,it),   ul,   ur(kt,jt,it))
          call delta4(fdz,   v(kt-1:kt+2,jt,it),   vl,   vr(kt,jt,it))
          call delta4(fdz,   w(kt-1:kt+2,jt,it),   wl,   wr(kt,jt,it))
          call delta4(fdz,   p(kt-1:kt+2,jt,it),   pl,   pr(kt,jt,it))
        endif
      else
        if (fdz <= threshold) then
          block
            real(8) tmp(2)
            tmp = T(i,j,k:k+1)
            associate(ww => w)
              G(i-1,:,j-1,k) = KEEP(id_accuracy2, &
                                    rho(kt:kt+1,jt,it), u(kt:kt+1,jt,it), &
                                      v(kt:kt+1,jt,it), w(kt:kt+1,jt,it), &
                                     ww(kt:kt+1,jt,it), p(kt:kt+1,jt,it), tmp, Normal_z)
            end associate
          end block
        else
          rhol = rho(kt,jt,it); rhor(kt,jt,it) = rho(kt+1,jt,it)
            ul =   u(kt,jt,it);   ur(kt,jt,it) =   u(kt+1,jt,it)
            vl =   v(kt,jt,it);   vr(kt,jt,it) =   v(kt+1,jt,it)
            wl =   w(kt,jt,it);   wr(kt,jt,it) =   w(kt+1,jt,it)
            pl =   p(kt,jt,it);   pr(kt,jt,it) =   p(kt+1,jt,it)
        endif
      endif
      call syncthreads()
      rho(kt,jt,it) = rhol
        u(kt,jt,it) = ul
        v(kt,jt,it) = vl
        w(kt,jt,it) = wl
        p(kt,jt,it) = pl
      end block
    if (fdz > threshold) then
      associate(un1 => w(kt,jt,it), un2 => wr(kt,jt,it))
        call SLAU(id_slau, rho(kt,jt,it), rhor(kt,jt,it), u(kt,jt,it), ur(kt,jt,it), v(kt,jt,it), vr(kt,jt,it), &
                  w(kt,jt,it), wr(kt,jt,it), un1, un2, p(kt,jt,it), pr(kt,jt,it), Normal_z, 1.d0, &
                  G(i-1,1,j-1,k), G(i-1,2,j-1,k), G(i-1,3,j-1,k), G(i-1,4,j-1,k), G(i-1,5,j-1,k))
      end associate
    endif
  end subroutine calc_hybrid_z6


  attributes(global) subroutine calc_hybrid_x4(id_accuracy, nx, ny, nz, Q, T, sensor, E)
    use mod_constant, only : Normal_x
    integer(kind=4), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: E(nx-1,5,ny-2,nz-2)
    integer i, j, k, it, jt, kt, ii, i_base
    real(8), dimension(0:threadsE%x+2,threadsE%y,threadsE%z), shared :: rho,  u,  v,  w,  p
    real(8), dimension(  threadsE%x,  threadsE%y,threadsE%z), shared :: rhor, ur, vr, wr, pr
    real(8) fdx
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    do ii = it-1, threadsE%x+2, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(ii,jt,kt) = Q(i,1,j,k)
          u(ii,jt,kt) = Q(i,2,j,k)
          v(ii,jt,kt) = Q(i,3,j,k)
          w(ii,jt,kt) = Q(i,4,j,k)
          p(ii,jt,kt) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    block
      real(8) rhol, ul, vl, wl, pl
      fdx = 0.5d0 * (sensor(i,j,k) + sensor(i+1,j,k))
      if (2 <= i .and. i <= nx-2) then
        if (fdx <= threshold) then
          block
            real(8) tmp(4)
            tmp = T(i-1:i+2,j,k)
            associate(uu => u)
              E(i,:,j-1,k-1) = KEEP(id_accuracy, &
                                    rho(it-1:it+2,jt,kt), u(it-1:it+2,jt,kt), &
                                      v(it-1:it+2,jt,kt), w(it-1:it+2,jt,kt), &
                                     uu(it-1:it+2,jt,kt), p(it-1:it+2,jt,kt), tmp, Normal_x)
            end associate
          end block
        else
          call delta4(fdx, rho(it-1:it+2,jt,kt), rhol, rhor(it,jt,kt))
          call delta4(fdx,   u(it-1:it+2,jt,kt),   ul,   ur(it,jt,kt))
          call delta4(fdx,   v(it-1:it+2,jt,kt),   vl,   vr(it,jt,kt))
          call delta4(fdx,   w(it-1:it+2,jt,kt),   wl,   wr(it,jt,kt))
          call delta4(fdx,   p(it-1:it+2,jt,kt),   pl,   pr(it,jt,kt))
        endif
      else
        if (fdx <= threshold) then
          block
            real(8) tmp(2)
            tmp = T(i:i+1,j,k)
            associate(uu => u)
              E(i,:,j-1,k-1) = KEEP(id_accuracy2, &
                                    rho(it:it+1,jt,kt), u(it:it+1,jt,kt), &
                                      v(it:it+1,jt,kt), w(it:it+1,jt,kt), &
                                     uu(it:it+1,jt,kt), p(it:it+1,jt,kt), tmp, Normal_x)
            end associate
          end block
        else
          rhol = rho(it,jt,kt); rhor(it,jt,kt) = rho(it+1,jt,kt)
            ul =   u(it,jt,kt);   ur(it,jt,kt) =   u(it+1,jt,kt)
            vl =   v(it,jt,kt);   vr(it,jt,kt) =   v(it+1,jt,kt)
            wl =   w(it,jt,kt);   wr(it,jt,kt) =   w(it+1,jt,kt)
            pl =   p(it,jt,kt);   pr(it,jt,kt) =   p(it+1,jt,kt)
        endif
      endif
      call syncthreads()
      rho(it,jt,kt) = rhol
        u(it,jt,kt) = ul
        v(it,jt,kt) = vl
        w(it,jt,kt) = wl
        p(it,jt,kt) = pl
    end block
    if (fdx > threshold) then
      associate(un1 => u(it,jt,kt), un2 => ur(it,jt,kt))
        call SLAU(id_slau, rho(it,jt,kt), rhor(it,jt,kt), u(it,jt,kt), ur(it,jt,kt), v(it,jt,kt), vr(it,jt,kt), &
                  w(it,jt,kt), wr(it,jt,kt), un1, un2, p(it,jt,kt), pr(it,jt,kt), Normal_x, 1.d0, &
                  E(i,1,j-1,k-1), E(i,2,j-1,k-1), E(i,3,j-1,k-1), E(i,4,j-1,k-1), E(i,5,j-1,k-1))
      end associate
    endif
  end subroutine calc_hybrid_x4


  attributes(global) subroutine calc_hybrid_y4(id_accuracy, nx, ny, nz, Q, T, sensor, F)
    use mod_constant, only : Normal_y
    integer(kind=4), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: F(nx-2,5,ny-1,nz-2)
    integer i, j, k, it, jt, kt, jj, j_base
    real(8), dimension(0:threadsF%y+2,threadsF%x,threadsF%z), shared :: rho,  u,  v,  w,  p
    real(8), dimension(  threadsF%y+2,threadsF%x,threadsF%z), shared :: rhor, ur, vr, wr, pr
    real(8) fdy
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    do jj = jt-1, threadsF%y+2, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(jj,it,kt) = Q(i,1,j,k)
          u(jj,it,kt) = Q(i,2,j,k)
          v(jj,it,kt) = Q(i,3,j,k)
          w(jj,it,kt) = Q(i,4,j,k)
          p(jj,it,kt) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    block
      real(8) rhol, ul, vl, wl, pl
      fdy = 0.5d0 * (sensor(i,j,k) + sensor(i,j+1,k))
      if (2 <= j .and. j <= ny-2) then
        if (fdy <= threshold) then
          block
            real(8) tmp(4)
            tmp = T(i,j-1:j+2,k)
            associate(vv => v)
              F(i-1,:,j,k-1) = KEEP(id_accuracy, &
                                    rho(jt-1:jt+2,it,kt), u(jt-1:jt+2,it,kt), &
                                      v(jt-1:jt+2,it,kt), w(jt-1:jt+2,it,kt), &
                                     vv(jt-1:jt+2,it,kt), p(jt-1:jt+2,it,kt), tmp, Normal_y)
            end associate
          end block
        else
          call delta4(fdy, rho(jt-1:jt+2,it,kt), rhol, rhor(jt,it,kt))
          call delta4(fdy,   u(jt-1:jt+2,it,kt),   ul,   ur(jt,it,kt))
          call delta4(fdy,   v(jt-1:jt+2,it,kt),   vl,   vr(jt,it,kt))
          call delta4(fdy,   w(jt-1:jt+2,it,kt),   wl,   wr(jt,it,kt))
          call delta4(fdy,   p(jt-1:jt+2,it,kt),   pl,   pr(jt,it,kt))
        endif
      else
        if (fdy <= threshold) then
          block
            real(8) tmp(2)
            tmp = T(i,j:j+1,k)
            associate(vv => v)
              F(i-1,:,j,k-1) = KEEP(id_accuracy2, &
                                    rho(jt:jt+1,it,kt), u(jt:jt+1,it,kt), &
                                      v(jt:jt+1,it,kt), w(jt:jt+1,it,kt), &
                                     vv(jt:jt+1,it,kt), p(jt:jt+1,it,kt), tmp, Normal_y)
            end associate
          end block
        else
          rhol = rho(jt,it,kt); rhor(jt,it,kt) = rho(jt+1,it,kt)
            ul =   u(jt,it,kt);   ur(jt,it,kt) =   u(jt+1,it,kt)
            vl =   v(jt,it,kt);   vr(jt,it,kt) =   v(jt+1,it,kt)
            wl =   w(jt,it,kt);   wr(jt,it,kt) =   w(jt+1,it,kt)
            pl =   p(jt,it,kt);   pr(jt,it,kt) =   p(jt+1,it,kt)
        endif
      endif
      call syncthreads()
      rho(jt,it,kt) = rhol
        u(jt,it,kt) = ul
        v(jt,it,kt) = vl
        w(jt,it,kt) = wl
        p(jt,it,kt) = pl
    end block
    if (fdy > threshold) then
      associate(un1 => v(jt,it,kt), un2 => vr(jt,it,kt))
        call SLAU(id_slau, rho(jt,it,kt), rhor(jt,it,kt), u(jt,it,kt), ur(jt,it,kt), v(jt,it,kt), vr(jt,it,kt), &
                  w(jt,it,kt), wr(jt,it,kt), un1, un2, p(jt,it,kt), pr(jt,it,kt), Normal_y, 1.d0, &
                  F(i-1,1,j,k-1), F(i-1,2,j,k-1), F(i-1,3,j,k-1), F(i-1,4,j,k-1), F(i-1,5,j,k-1))
      end associate
    endif
  end subroutine calc_hybrid_y4


  attributes(global) subroutine calc_hybrid_z4(id_accuracy, nx, ny, nz, Q, T, sensor, G)
    use mod_constant, only : Normal_z
    integer(kind=4), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: G(nx-2,5,ny-2,nz-1)
    integer i, j, k, it, jt, kt, kk, k_base
    real(8), dimension(0:threadsG%z+2,threadsG%y,threadsG%x), shared :: rho,  u,  v,  w,  p
    real(8), dimension(  threadsG%z+2,threadsG%y,threadsG%x), shared :: rhor, ur, vr, wr, pr
    real(8) fdz
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    do kk = kt-1, threadsG%z+2, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(kk,jt,it) = Q(i,1,j,k)
          u(kk,jt,it) = Q(i,2,j,k)
          v(kk,jt,it) = Q(i,3,j,k)
          w(kk,jt,it) = Q(i,4,j,k)
          p(kk,jt,it) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    block
      real(8) rhol, ul, vl, wl, pl
      fdz = 0.5d0 * (sensor(i,j,k) + sensor(i,j,k+1))
      if (2 <= k .and. k <= nz-2) then
        if (fdz <= threshold) then
          block
            real(8) tmp(4)
            tmp = T(i,j,k-1:k+2)
            associate(ww => w)
              G(i-1,:,j-1,k) = KEEP(id_accuracy, &
                                    rho(kt-1:kt+2,jt,it), u(kt-1:kt+2,jt,it), &
                                      v(kt-1:kt+2,jt,it), w(kt-1:kt+2,jt,it), &
                                     ww(kt-1:kt+2,jt,it), p(kt-1:kt+2,jt,it), tmp, Normal_z)
            end associate
          end block
        else
          call delta4(fdz, rho(kt-1:kt+2,jt,it), rhol, rhor(kt,jt,it))
          call delta4(fdz,   u(kt-1:kt+2,jt,it),   ul,   ur(kt,jt,it))
          call delta4(fdz,   v(kt-1:kt+2,jt,it),   vl,   vr(kt,jt,it))
          call delta4(fdz,   w(kt-1:kt+2,jt,it),   wl,   wr(kt,jt,it))
          call delta4(fdz,   p(kt-1:kt+2,jt,it),   pl,   pr(kt,jt,it))
        endif
      else
        if (fdz <= threshold) then
          block
            real(8) tmp(2)
            tmp = T(i,j,k:k+1)
            associate(ww => w)
              G(i-1,:,j-1,k) = KEEP(id_accuracy2, &
                                    rho(kt:kt+1,jt,it), u(kt:kt+1,jt,it), &
                                      v(kt:kt+1,jt,it), w(kt:kt+1,jt,it), &
                                     ww(kt:kt+1,jt,it), p(kt:kt+1,jt,it), tmp, Normal_z)
            end associate
          end block
        else
          rhol = rho(kt,jt,it); rhor(kt,jt,it) = rho(kt+1,jt,it)
            ul =   u(kt,jt,it);   ur(kt,jt,it) =   u(kt+1,jt,it)
            vl =   v(kt,jt,it);   vr(kt,jt,it) =   v(kt+1,jt,it)
            wl =   w(kt,jt,it);   wr(kt,jt,it) =   w(kt+1,jt,it)
            pl =   p(kt,jt,it);   pr(kt,jt,it) =   p(kt+1,jt,it)
        endif
      endif
      call syncthreads()
      rho(kt,jt,it) = rhol
        u(kt,jt,it) = ul
        v(kt,jt,it) = vl
        w(kt,jt,it) = wl
        p(kt,jt,it) = pl
    end block
    if (fdz > threshold) then
      associate(un1 => w(kt,jt,it), un2 => wr(kt,jt,it))
        call SLAU(id_slau, rho(kt,jt,it), rhor(kt,jt,it), u(kt,jt,it), ur(kt,jt,it), v(kt,jt,it), vr(kt,jt,it), &
                  w(kt,jt,it), wr(kt,jt,it), un1, un2, p(kt,jt,it), pr(kt,jt,it), Normal_z, 1.d0, &
                  G(i-1,1,j-1,k), G(i-1,2,j-1,k), G(i-1,3,j-1,k), G(i-1,4,j-1,k), G(i-1,5,j-1,k))
      end associate
    endif
  end subroutine calc_hybrid_z4


  attributes(global) subroutine calc_hybrid_x2(id_accuracy, nx, ny, nz, Q, T, sensor, E)
    use mod_constant, only : Normal_x
    integer(kind=2), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: E(nx-1,5,ny-2,nz-2)
    integer i, j, k, it, jt, kt, ii, i_base
    real(8), dimension(threadsE%x+1,threadsE%y,threadsE%z), shared :: rho, u, v, w, p
    real(8) fdx
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    do ii = it, threadsE%x+1, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(ii,jt,kt) = Q(i,1,j,k)
          u(ii,jt,kt) = Q(i,2,j,k)
          v(ii,jt,kt) = Q(i,3,j,k)
          w(ii,jt,kt) = Q(i,4,j,k)
          p(ii,jt,kt) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    fdx = 0.5d0 * (sensor(i,j,k) + sensor(i+1,j,k))
    if (fdx <= threshold) then
      block
        real(8) tmp(2)
        tmp = T(i:i+1,j,k)
        associate(uu => u)
          E(i,:,j-1,k-1) = KEEP(id_accuracy, &
                                rho(it:it+1,jt,kt), u(it:it+1,jt,kt), &
                                  v(it:it+1,jt,kt), w(it:it+1,jt,kt), &
                                 uu(it:it+1,jt,kt), p(it:it+1,jt,kt), tmp, Normal_x)
        end associate
      end block
    else
      associate(un1 => u(it,jt,kt), un2 => u(it+1,jt,kt))
        call SLAU(id_slau, rho(it,jt,kt), rho(it+1,jt,kt), u(it,jt,kt), u(it+1,jt,kt), v(it,jt,kt), v(it+1,jt,kt), &
                  w(it,jt,kt), w(it+1,jt,kt), un1, un2, p(it,jt,kt), p(it+1,jt,kt), Normal_x, 1.d0, &
                  E(i,1,j-1,k-1), E(i,2,j-1,k-1), E(i,3,j-1,k-1), E(i,4,j-1,k-1), E(i,5,j-1,k-1))
      end associate
    endif
  end subroutine calc_hybrid_x2


  attributes(global) subroutine calc_hybrid_y2(id_accuracy, nx, ny, nz, Q, T, sensor, F)
    use mod_constant, only : Normal_y
    integer(kind=2), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: F(nx-2,5,ny-1,nz-2)
    integer i, j, k, it, jt, kt, jj, j_base
    real(8), dimension(threadsF%y+1,threadsF%x,threadsF%z), shared :: rho, u, v, w, p
    real(8) fdy
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    do jj = jt, threadsF%y+1, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(jj,it,kt) = Q(i,1,j,k)
          u(jj,it,kt) = Q(i,2,j,k)
          v(jj,it,kt) = Q(i,3,j,k)
          w(jj,it,kt) = Q(i,4,j,k)
          p(jj,it,kt) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    fdy = 0.5d0 * (sensor(i,j,k) + sensor(i,j+1,k))
    if (fdy <= threshold) then
      block
        real(8) tmp(2)
        tmp = T(i,j:j+1,k)
        associate(vv => v)
          F(i-1,:,j,k-1) = KEEP(id_accuracy, &
                                rho(jt:jt+1,it,kt), u(jt:jt+1,it,kt), &
                                  v(jt:jt+1,it,kt), w(jt:jt+1,it,kt), &
                                 vv(jt:jt+1,it,kt), p(jt:jt+1,it,kt), tmp, Normal_y)
        end associate
      end block
    else
      associate(un1 => v(jt,it,kt), un2 => v(jt+1,it,kt))
        call SLAU(id_slau, rho(jt,it,kt), rho(jt+1,it,kt), u(jt,it,kt), u(jt+1,it,kt), v(jt,it,kt), v(jt+1,it,kt), &
                  w(jt,it,kt), w(jt+1,it,kt), un1, un2, p(jt,it,kt), p(jt+1,it,kt), Normal_y, 1.d0, &
                  F(i-1,1,j,k-1), F(i-1,2,j,k-1), F(i-1,3,j,k-1), F(i-1,4,j,k-1), F(i-1,5,j,k-1))
      end associate
    endif
  end subroutine calc_hybrid_y2


  attributes(global) subroutine calc_hybrid_z2(id_accuracy, nx, ny, nz, Q, T, sensor, G)
    use mod_constant, only : Normal_z
    integer(kind=2), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), T(nx,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: G(nx-2,5,ny-2,nz-1)
    integer i, j, k, it, jt, kt, kk, k_base
    real(8), dimension(threadsG%z+1,threadsG%y,threadsG%x), shared :: rho, u, v, w, p
    real(8) fdz
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    do kk = kt, threadsG%z+1, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(kk,jt,it) = Q(i,1,j,k)
          u(kk,jt,it) = Q(i,2,j,k)
          v(kk,jt,it) = Q(i,3,j,k)
          w(kk,jt,it) = Q(i,4,j,k)
          p(kk,jt,it) = Q(i,5,j,k)
      endif
    enddo
    call syncthreads()
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    fdz = 0.5d0 * (sensor(i,j,k) + sensor(i,j,k+1))
    if (fdz <= threshold) then
      block
        real(8) tmp(2)
        tmp = T(i,j,k:k+1)
        associate(ww => w)
          G(i-1,:,j-1,k) = KEEP(id_accuracy, &
                                rho(kt:kt+1,jt,it), u(kt:kt+1,jt,it), &
                                  v(kt:kt+1,jt,it), w(kt:kt+1,jt,it), &
                                 ww(kt:kt+1,jt,it), p(kt:kt+1,jt,it), tmp, Normal_z)
        end associate
      end block
    else
      associate(un1 => w(kt,jt,it), un2 => w(kt+1,jt,it))
        call SLAU(id_slau, rho(kt,jt,it), rho(kt+1,jt,it), u(kt,jt,it), u(kt+1,jt,it), v(kt,jt,it), v(kt+1,jt,it), &
                  w(kt,jt,it), w(kt+1,jt,it), un1, un2, p(kt,jt,it), p(kt+1,jt,it), Normal_z, 1.d0, &
                  G(i-1,1,j-1,k), G(i-1,2,j-1,k), G(i-1,3,j-1,k), G(i-1,4,j-1,k), G(i-1,5,j-1,k))
      end associate
    endif
  end subroutine calc_hybrid_z2
end module calc_hybrid_kernel

