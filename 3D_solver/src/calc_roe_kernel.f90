module calc_roe_kernel
  use mod_globals, only : gamma, threadsE, threadsF, threadsG
  use mod_constant, only : gamma_1, over_gamma_1
  use calc_muscl
  use calc_hybrid
  implicit none
  private
  public calc_roe_x, calc_roe_y, calc_roe_z
  
  interface calc_roe_x
    module procedure calc_roe_x2, calc_roe_x4, calc_roe_x6
  end interface calc_roe_x

  interface calc_roe_y
    module procedure calc_roe_y2, calc_roe_y4, calc_roe_y6
  end interface calc_roe_y

  interface calc_roe_z
    module procedure calc_roe_z2, calc_roe_z4, calc_roe_z6
  end interface calc_roe_z
contains
  include 'calc_roe_3d.f90'

  attributes(global) subroutine calc_roe_x6(id_accuracy, nx, ny, nz, Q, sensor, E)
    integer(kind=8), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: E(nx-1,5,ny-2,nz-2)
    integer i, j, k, it, jt, kt, ii, i_base
    real(8), dimension(-1:threadsE%x+3,threadsE%y,threadsE%z), shared :: rho,  u,  v,  w,  p
    real(8), dimension(   threadsE%x,  threadsE%y,threadsE%z), shared :: rhor, ur, vr, wr, pr
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
      real(8) rhol, ul, vl, wl, pl, fdx
      fdx = 0.5d0 * (sensor(i,j,k) + sensor(i+1,j,k))
      if (3 <= i .and. i <= nx-3) then
        call delta6(fdx, rho(it-2:it+3,jt,kt), rhol, rhor(it,jt,kt))
        call delta6(fdx,   u(it-2:it+3,jt,kt),   ul,   ur(it,jt,kt))
        call delta6(fdx,   v(it-2:it+3,jt,kt),   vl,   vr(it,jt,kt))
        call delta6(fdx,   w(it-2:it+3,jt,kt),   wl,   wr(it,jt,kt))
        call delta6(fdx,   p(it-2:it+3,jt,kt),   pl,   pr(it,jt,kt))
      elseif (2 <= i .and. i <= nx-2) then
        call delta4(fdx, rho(it-1:it+2,jt,kt), rhol, rhor(it,jt,kt))
        call delta4(fdx,   u(it-1:it+2,jt,kt),   ul,   ur(it,jt,kt))
        call delta4(fdx,   v(it-1:it+2,jt,kt),   vl,   vr(it,jt,kt))
        call delta4(fdx,   w(it-1:it+2,jt,kt),   wl,   wr(it,jt,kt))
        call delta4(fdx,   p(it-1:it+2,jt,kt),   pl,   pr(it,jt,kt))
      else
        rhol = rho(it,jt,kt); rhor(it,jt,kt) = rho(it+1,jt,kt)
          ul =   u(it,jt,kt);   ur(it,jt,kt) =   u(it+1,jt,kt)
          vl =   v(it,jt,kt);   vr(it,jt,kt) =   v(it+1,jt,kt)
          wl =   w(it,jt,kt);   wr(it,jt,kt) =   w(it+1,jt,kt)
          pl =   p(it,jt,kt);   pr(it,jt,kt) =   p(it+1,jt,kt)
      endif
      call syncthreads()
      rho(it,jt,kt) = rhol
        u(it,jt,kt) = ul
        v(it,jt,kt) = vl
        w(it,jt,kt) = wl
        p(it,jt,kt) = pl
    end block
    call Roe(rho(it,jt,kt), rhor(it,jt,kt), u(it,jt,kt), ur(it,jt,kt), v(it,jt,kt), vr(it,jt,kt), &
             w(it,jt,kt), wr(it,jt,kt), p(it,jt,kt), pr(it,jt,kt), &
             E(i,1,j-1,k-1), E(i,2,j-1,k-1), E(i,3,j-1,k-1), E(i,4,j-1,k-1), E(i,5,j-1,k-1))
  end subroutine calc_roe_x6


  attributes(global) subroutine calc_roe_y6(id_accuracy, nx, ny, nz, Q, sensor, F)
    integer(kind=8), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: F(nx-2,5,ny-1,nz-2)
    integer i, j, k, it, jt, kt, jj, j_base
    real(8), dimension(-1:threadsF%y+3,threadsF%x,threadsF%z), shared :: rho,  u,  v,  w,  p
    real(8), dimension(   threadsF%y,  threadsF%x,threadsF%z), shared :: rhor, ur, vr, wr, pr
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
      real(8) rhol, ul, vl, wl, pl, fdy
      fdy = 0.5d0 * (sensor(i,j,k) + sensor(i,j+1,k))
      if (3 <= j .and. j <= ny-3 .and. 8 <= kind(id_accuracy)) then
        call delta6(fdy, rho(jt-2:jt+3,it,kt), rhol, rhor(jt,it,kt))
        call delta6(fdy,   u(jt-2:jt+3,it,kt),   ul,   ur(jt,it,kt))
        call delta6(fdy,   v(jt-2:jt+3,it,kt),   vl,   vr(jt,it,kt))
        call delta6(fdy,   w(jt-2:jt+3,it,kt),   wl,   wr(jt,it,kt))
        call delta6(fdy,   p(jt-2:jt+3,it,kt),   pl,   pr(jt,it,kt))
      elseif (2 <= j .and. j <= ny-2) then
        call delta4(fdy, rho(jt-1:jt+2,it,kt), rhol, rhor(jt,it,kt))
        call delta4(fdy,   u(jt-1:jt+2,it,kt),   ul,   ur(jt,it,kt))
        call delta4(fdy,   v(jt-1:jt+2,it,kt),   vl,   vr(jt,it,kt))
        call delta4(fdy,   w(jt-1:jt+2,it,kt),   wl,   wr(jt,it,kt))
        call delta4(fdy,   p(jt-1:jt+2,it,kt),   pl,   pr(jt,it,kt))
      else
        rhol = rho(jt,it,kt); rhor(jt,it,kt) = rho(jt+1,it,kt)
          ul =   u(jt,it,kt);   ur(jt,it,kt) =   u(jt+1,it,kt)
          vl =   v(jt,it,kt);   vr(jt,it,kt) =   v(jt+1,it,kt)
          wl =   w(jt,it,kt);   wr(jt,it,kt) =   w(jt+1,it,kt)
          pl =   p(jt,it,kt);   pr(jt,it,kt) =   p(jt+1,it,kt)
      endif
      call syncthreads()
      rho(jt,it,kt) = rhol
        u(jt,it,kt) = ul
        v(jt,it,kt) = vl
        w(jt,it,kt) = wl
        p(jt,it,kt) = pl
    end block
    call Roe(rho(jt,it,kt), rhor(jt,it,kt), v(jt,it,kt), vr(jt,it,kt), w(jt,it,kt), wr(jt,it,kt), &
             u(jt,it,kt), ur(jt,it,kt), p(jt,it,kt), pr(jt,it,kt), &
             F(i-1,1,j,k-1), F(i-1,3,j,k-1), F(i-1,4,j,k-1), F(i-1,2,j,k-1), F(i-1,5,j,k-1))
  end subroutine calc_roe_y6


  attributes(global) subroutine calc_roe_z6(id_accuracy, nx, ny, nz, Q, sensor, G)
    integer(kind=8), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: G(nx-2,5,ny-2,nz-1)
    integer i, j, k, it, jt, kt, kk, k_base
    real(8), dimension(-1:threadsG%z+3,threadsG%y,threadsG%x), shared :: rho,  u,  v,  w,  p
    real(8), dimension(   threadsG%z,  threadsG%y,threadsG%x), shared :: rhor, ur, vr, wr, pr
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
      real(8) rhol, ul, vl, wl, pl, fdz
      fdz = 0.5d0 * (sensor(i,j,k) + sensor(i,j,k+1))
      if (3 <= k .and. k <= nz-3 .and. 8 <= kind(id_accuracy)) then
        call delta6(fdz, rho(kt-2:kt+3,jt,it), rhol, rhor(kt,jt,it))
        call delta6(fdz,   u(kt-2:kt+3,jt,it),   ul,   ur(kt,jt,it))
        call delta6(fdz,   v(kt-2:kt+3,jt,it),   vl,   vr(kt,jt,it))
        call delta6(fdz,   w(kt-2:kt+3,jt,it),   wl,   wr(kt,jt,it))
        call delta6(fdz,   p(kt-2:kt+3,jt,it),   pl,   pr(kt,jt,it))
      elseif (2 <= k .and. k <= nz-2) then
        call delta4(fdz, rho(kt-1:kt+2,jt,it), rhol, rhor(kt,jt,it))
        call delta4(fdz,   u(kt-1:kt+2,jt,it),   ul,   ur(kt,jt,it))
        call delta4(fdz,   v(kt-1:kt+2,jt,it),   vl,   vr(kt,jt,it))
        call delta4(fdz,   w(kt-1:kt+2,jt,it),   wl,   wr(kt,jt,it))
        call delta4(fdz,   p(kt-1:kt+2,jt,it),   pl,   pr(kt,jt,it))
      else
        rhol = rho(kt,jt,it); rhor(kt,jt,it) = rho(kt+1,jt,it)
          ul =   u(kt,jt,it);   ur(kt,jt,it) =   u(kt+1,jt,it)
          vl =   v(kt,jt,it);   vr(kt,jt,it) =   v(kt+1,jt,it)
          wl =   w(kt,jt,it);   wr(kt,jt,it) =   w(kt+1,jt,it)
          pl =   p(kt,jt,it);   pr(kt,jt,it) =   p(kt+1,jt,it)
      endif
      call syncthreads()
      rho(kt,jt,it) = rhol
        u(kt,jt,it) = ul
        v(kt,jt,it) = vl
        w(kt,jt,it) = wl
        p(kt,jt,it) = pl
    end block
    call Roe(rho(kt,jt,it), rhor(kt,jt,it), w(kt,jt,it), wr(kt,jt,it), u(kt,jt,it), ur(kt,jt,it), &
             v(kt,jt,it), vr(kt,jt,it), p(kt,jt,it), pr(kt,jt,it), &
             G(i-1,1,j-1,k), G(i-1,4,j-1,k), G(i-1,2,j-1,k), G(i-1,3,j-1,k), G(i-1,5,j-1,k))
  end subroutine calc_roe_z6


  attributes(global) subroutine calc_roe_x4(id_accuracy, nx, ny, nz, Q, sensor, E)
    integer(kind=4), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: E(nx-1,5,ny-2,nz-2)
    integer i, j, k, it, jt, kt, ii, i_base
    real(8), dimension(0:threadsE%x+2,threadsE%y,threadsE%z), shared :: rho,  u,  v,  w,  p
    real(8), dimension(  threadsE%x,  threadsE%y,threadsE%z), shared :: rhor, ur, vr, wr, pr
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
      real(8) rhol, ul, vl, wl, pl, fdx
      fdx = 0.5d0 * (sensor(i,j,k) + sensor(i+1,j,k))
      if (2 <= i .and. i <= nx-2) then
        call delta4(fdx, rho(it-1:it+2,jt,kt), rhol, rhor(it,jt,kt))
        call delta4(fdx,   u(it-1:it+2,jt,kt),   ul,   ur(it,jt,kt))
        call delta4(fdx,   v(it-1:it+2,jt,kt),   vl,   vr(it,jt,kt))
        call delta4(fdx,   w(it-1:it+2,jt,kt),   wl,   wr(it,jt,kt))
        call delta4(fdx,   p(it-1:it+2,jt,kt),   pl,   pr(it,jt,kt))
      else
        rhol = rho(it,jt,kt); rhor(it,jt,kt) = rho(it+1,jt,kt)
          ul =   u(it,jt,kt);   ur(it,jt,kt) =   u(it+1,jt,kt)
          vl =   v(it,jt,kt);   vr(it,jt,kt) =   v(it+1,jt,kt)
          wl =   w(it,jt,kt);   wr(it,jt,kt) =   w(it+1,jt,kt)
          pl =   p(it,jt,kt);   pr(it,jt,kt) =   p(it+1,jt,kt)
      endif
      call syncthreads()
      rho(it,jt,kt) = rhol
        u(it,jt,kt) = ul
        v(it,jt,kt) = vl
        w(it,jt,kt) = wl
        p(it,jt,kt) = pl
    end block
    call Roe(rho(it,jt,kt), rhor(it,jt,kt), u(it,jt,kt), ur(it,jt,kt), v(it,jt,kt), vr(it,jt,kt), &
             w(it,jt,kt), wr(it,jt,kt), p(it,jt,kt), pr(it,jt,kt), &
             E(i,1,j-1,k-1), E(i,2,j-1,k-1), E(i,3,j-1,k-1), E(i,4,j-1,k-1), E(i,5,j-1,k-1))
  end subroutine calc_roe_x4


  attributes(global) subroutine calc_roe_y4(id_accuracy, nx, ny, nz, Q, sensor, F)
    integer(kind=4), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: F(nx-2,5,ny-1,nz-2)
    integer i, j, k, it, jt, kt, jj, j_base
    real(8), dimension(0:threadsF%y+2,threadsF%x,threadsF%z), shared :: rho,  u,  v,  w,  p
    real(8), dimension(  threadsF%y+2,threadsF%x,threadsF%z), shared :: rhor, ur, vr, wr, pr
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
      real(8) rhol, ul, vl, wl, pl, fdy
      fdy = 0.5d0 * (sensor(i,j,k) + sensor(i,j+1,k))
      if (2 <= j .and. j <= ny-2) then
        call delta4(fdy, rho(jt-1:jt+2,it,kt), rhol, rhor(jt,it,kt))
        call delta4(fdy,   u(jt-1:jt+2,it,kt),   ul,   ur(jt,it,kt))
        call delta4(fdy,   v(jt-1:jt+2,it,kt),   vl,   vr(jt,it,kt))
        call delta4(fdy,   w(jt-1:jt+2,it,kt),   wl,   wr(jt,it,kt))
        call delta4(fdy,   p(jt-1:jt+2,it,kt),   pl,   pr(jt,it,kt))
      else
        rhol = rho(jt,it,kt); rhor(jt,it,kt) = rho(jt+1,it,kt)
          ul =   u(jt,it,kt);   ur(jt,it,kt) =   u(jt+1,it,kt)
          vl =   v(jt,it,kt);   vr(jt,it,kt) =   v(jt+1,it,kt)
          wl =   w(jt,it,kt);   wr(jt,it,kt) =   w(jt+1,it,kt)
          pl =   p(jt,it,kt);   pr(jt,it,kt) =   p(jt+1,it,kt)
      endif
      call syncthreads()
      rho(jt,it,kt) = rhol
        u(jt,it,kt) = ul
        v(jt,it,kt) = vl
        w(jt,it,kt) = wl
        p(jt,it,kt) = pl
    end block
    call Roe(rho(jt,it,kt), rhor(jt,it,kt), v(jt,it,kt), vr(jt,it,kt), w(jt,it,kt), wr(jt,it,kt), &
             u(jt,it,kt), ur(jt,it,kt), p(jt,it,kt), pr(jt,it,kt), &
             F(i-1,1,j,k-1), F(i-1,3,j,k-1), F(i-1,4,j,k-1), F(i-1,2,j,k-1), F(i-1,5,j,k-1))
  end subroutine calc_roe_y4


  attributes(global) subroutine calc_roe_z4(id_accuracy, nx, ny, nz, Q, sensor, G)
    integer(kind=4), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: G(nx-2,5,ny-2,nz-1)
    integer i, j, k, it, jt, kt, kk, k_base
    real(8), dimension(0:threadsG%z+2,threadsG%y,threadsG%x), shared :: rho,  u,  v,  w,  p
    real(8), dimension(  threadsG%z+2,threadsG%y,threadsG%x), shared :: rhor, ur, vr, wr, pr
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
      real(8) rhol, ul, vl, wl, pl, fdz
      if (2 <= k .and. k <= nz-2) then
        call delta4(fdz, rho(kt-1:kt+2,jt,it), rhol, rhor(kt,jt,it))
        call delta4(fdz,   u(kt-1:kt+2,jt,it),   ul,   ur(kt,jt,it))
        call delta4(fdz,   v(kt-1:kt+2,jt,it),   vl,   vr(kt,jt,it))
        call delta4(fdz,   w(kt-1:kt+2,jt,it),   wl,   wr(kt,jt,it))
        call delta4(fdz,   p(kt-1:kt+2,jt,it),   pl,   pr(kt,jt,it))
      else
        rhol = rho(kt,jt,it); rhor(kt,jt,it) = rho(kt+1,jt,it)
          ul =   u(kt,jt,it);   ur(kt,jt,it) =   u(kt+1,jt,it)
          vl =   v(kt,jt,it);   vr(kt,jt,it) =   v(kt+1,jt,it)
          wl =   w(kt,jt,it);   wr(kt,jt,it) =   w(kt+1,jt,it)
          pl =   p(kt,jt,it);   pr(kt,jt,it) =   p(kt+1,jt,it)
      endif
      call syncthreads()
      rho(kt,jt,it) = rhol
        u(kt,jt,it) = ul
        v(kt,jt,it) = vl
        w(kt,jt,it) = wl
        p(kt,jt,it) = pl
    end block
    call Roe(rho(kt,jt,it), rhor(kt,jt,it), w(kt,jt,it), wr(kt,jt,it), u(kt,jt,it), ur(kt,jt,it), &
             v(kt,jt,it), vr(kt,jt,it), p(kt,jt,it), pr(kt,jt,it), &
             G(i-1,1,j-1,k), G(i-1,4,j-1,k), G(i-1,2,j-1,k), G(i-1,3,j-1,k), G(i-1,5,j-1,k))
  end subroutine calc_roe_z4


  attributes(global) subroutine calc_roe_x2(id_accuracy, nx, ny, nz, Q, sensor, E)
    integer(kind=2), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: E(nx-1,5,ny-2,nz-2)
    integer i, j, k, it, jt, kt, ii, i_base
    real(8), dimension(threadsE%x+1,threadsE%y,threadsE%z), shared :: rho, u, v, w, p
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
    call Roe(rho(it,jt,kt), rho(it+1,jt,kt), u(it,jt,kt), u(it+1,jt,kt), v(it,jt,kt), v(it+1,jt,kt), &
             w(it,jt,kt), w(it+1,jt,kt), p(it,jt,kt), p(it+1,jt,kt), &
             E(i,1,j-1,k-1), E(i,2,j-1,k-1), E(i,3,j-1,k-1), E(i,4,j-1,k-1), E(i,5,j-1,k-1))
  end subroutine calc_roe_x2


  attributes(global) subroutine calc_roe_y2(id_accuracy, nx, ny, nz, Q, sensor, F)
    integer(kind=2), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: F(nx-2,5,ny-1,nz-2)
    integer i, j, k, it, jt, kt, jj, j_base
    real(8), dimension(threadsF%y+1,threadsF%x,threadsF%z), shared :: rho, u, v, w, p
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
    call Roe(rho(jt,it,kt), rho(jt+1,it,kt), v(jt,it,kt), v(jt+1,it,kt), w(jt,it,kt), w(jt+1,it,kt), &
             u(jt,it,kt), u(jt+1,it,kt), p(jt,it,kt), p(jt+1,it,kt), &
             F(i-1,1,j,k-1), F(i-1,3,j,k-1), F(i-1,4,j,k-1), F(i-1,2,j,k-1), F(i-1,5,j,k-1))
  end subroutine calc_roe_y2


  attributes(global) subroutine calc_roe_z2(id_accuracy, nx, ny, nz, Q, sensor, G)
    integer(kind=2), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), device, contiguous   :: Q(nx,5,ny,nz), sensor(nx,ny,nz)
    real(8), intent(out), device, contiguous  :: G(nx-2,5,ny-2,nz-1)
    integer i, j, k, it, jt, kt, kk, k_base
    real(8), dimension(threadsG%z+1,threadsG%y,threadsG%x), shared :: rho, u, v, w, p
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
    call Roe(rho(kt,jt,it), rho(kt+1,jt,it), w(kt,jt,it), w(kt+1,jt,it), u(kt,jt,it), u(kt+1,jt,it), &
             v(kt,jt,it), v(kt+1,jt,it), p(kt,jt,it), p(kt+1,jt,it), &
             G(i-1,1,j-1,k), G(i-1,4,j-1,k), G(i-1,2,j-1,k), G(i-1,3,j-1,k), G(i-1,5,j-1,k))
  end subroutine calc_roe_z2
end module calc_roe_kernel

