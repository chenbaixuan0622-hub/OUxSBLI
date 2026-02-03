module calc_keep_kernel
  use mod_globals, only : id_igr, threadsE, threadsF, threadsG
  use mod_constant, only : R_over_gamma_1, one_third, one_sixth, one_twelfth, two_third
  implicit none
  private
  public calc_keep_x, calc_keep_y, calc_keep_z
  real(8), parameter :: one_24        = 1.d0 / 24.d0
  real(8), parameter :: one_48        = 1.d0 / 48.d0
  real(8), parameter :: one_60        = 1.d0 / 60.d0
  real(8), parameter :: one_120       = 1.d0 / 120.d0
  real(8), parameter :: one_240       = 1.d0 / 240.d0
  real(8), parameter :: seven_twelfth = 7.d0 / 12.d0
  
  interface calc_keep_x
    module procedure calc_keep_x2, calc_keep_x4, calc_keep_x6
  end interface calc_keep_x

  interface calc_keep_y
    module procedure calc_keep_y2, calc_keep_y4, calc_keep_y6
  end interface calc_keep_y

  interface calc_keep_z
    module procedure calc_keep_z2, calc_keep_z4, calc_keep_z6
  end interface calc_keep_z
contains
  include 'calc_keep_3d.f90'
 
  attributes(global) subroutine calc_keep_x6(id_accuracy, nx, ny, nz, Q, T, E, sigma)
    use mod_constant, only : Normal_x
    integer(kind=8), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: E(5,nx-1,ny-2,nz-2)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt, ii, i_base
    real(8), dimension(-1:threadsE%x+3,threadsE%y,threadsE%z), shared :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    do ii = it-2, threadsE%x+3, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(ii,jt,kt) = Q(1,i,j,k)
          u(ii,jt,kt) = Q(2,i,j,k)
          v(ii,jt,kt) = Q(3,i,j,k)
          w(ii,jt,kt) = Q(4,i,j,k)
          p(ii,jt,kt) = Q(5,i,j,k)
        tmp(ii,jt,kt) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    associate(uu => u)
    if (3 <= i .and. i <= nx-3) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(6)
        sig = dble(sigma(i-2:i+3,j,k))
        E(:,i,j-1,k-1) = KEEP_IGR6(rho(it-2:it+3,jt,kt), u(it-2:it+3,jt,kt), &
                                     v(it-2:it+3,jt,kt), w(it-2:it+3,jt,kt), &
                                    uu(it-2:it+3,jt,kt), p(it-2:it+3,jt,kt), &
                                   tmp(it-2:it+3,jt,kt), sig, Normal_x)
        end block
      else
        E(:,i,j-1,k-1) = KEEP6(rho(it-2:it+3,jt,kt), u(it-2:it+3,jt,kt), &
                                 v(it-2:it+3,jt,kt), w(it-2:it+3,jt,kt), &
                                uu(it-2:it+3,jt,kt), p(it-2:it+3,jt,kt), &
                               tmp(it-2:it+3,jt,kt), Normal_x)
      endif
    elseif (2 <= i .and. i <= nx-2) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(4)
        sig = dble(sigma(i-1:i+2,j,k))
        E(:,i,j-1,k-1) = KEEP_IGR4(rho(it-1:it+2,jt,kt), u(it-1:it+2,jt,kt), &
                                     v(it-1:it+2,jt,kt), w(it-1:it+2,jt,kt), &
                                    uu(it-1:it+2,jt,kt), p(it-1:it+2,jt,kt), &
                                   tmp(it-1:it+2,jt,kt), sig, Normal_x)
        end block
      else
        E(:,i,j-1,k-1) = KEEP4(rho(it-1:it+2,jt,kt), u(it-1:it+2,jt,kt), &
                                 v(it-1:it+2,jt,kt), w(it-1:it+2,jt,kt), &
                                uu(it-1:it+2,jt,kt), p(it-1:it+2,jt,kt), &
                               tmp(it-1:it+2,jt,kt), Normal_x)
      endif
    else
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(2)
        sig = dble(sigma(i:i+1,j,k))
        E(:,i,j-1,k-1) = KEEP_IGR2(rho(it:it+1,jt,kt), u(it:it+1,jt,kt), &
                                     v(it:it+1,jt,kt), w(it:it+1,jt,kt), &
                                    uu(it:it+1,jt,kt), p(it:it+1,jt,kt), &
                                   tmp(it:it+1,jt,kt), sig, Normal_x)
        end block
      else
        E(:,i,j-1,k-1) = KEEP2(rho(it:it+1,jt,kt), u(it:it+1,jt,kt), &
                                 v(it:it+1,jt,kt), w(it:it+1,jt,kt), &
                                uu(it:it+1,jt,kt), p(it:it+1,jt,kt), &
                               tmp(it:it+1,jt,kt), Normal_x)
      endif
    endif
    end associate
  end subroutine calc_keep_x6


  attributes(global) subroutine calc_keep_y6(id_accuracy, nx, ny, nz, Q, T, F, sigma)
    use mod_constant, only : Normal_y
    integer(kind=8), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: F(5,nx-2,ny-1,nz-2)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt, jj, j_base
    real(8), dimension(-1:threadsF%y+3,threadsF%x,threadsF%z), shared :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    do jj = jt-2, threadsF%y+3, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(jj,it,kt) = Q(1,i,j,k)
          u(jj,it,kt) = Q(2,i,j,k)
          v(jj,it,kt) = Q(3,i,j,k)
          w(jj,it,kt) = Q(4,i,j,k)
          p(jj,it,kt) = Q(5,i,j,k)
        tmp(jj,it,kt) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    associate(vv => v)
    if (3 <= j .and. j <= ny-3) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(6)
        sig = dble(sigma(i,j-2:j+3,k))
        F(:,i-1,j,k-1) = KEEP_IGR6(rho(jt-2:jt+3,it,kt), u(jt-2:jt+3,it,kt), &
                                     v(jt-2:jt+3,it,kt), w(jt-2:jt+3,it,kt), &
                                    vv(jt-2:jt+3,it,kt), p(jt-2:jt+3,it,kt), &
                                   tmp(jt-2:jt+3,it,kt), sig, Normal_y)
        end block
      else
        F(:,i-1,j,k-1) = KEEP6(rho(jt-2:jt+3,it,kt), u(jt-2:jt+3,it,kt), &
                                 v(jt-2:jt+3,it,kt), w(jt-2:jt+3,it,kt), &
                                vv(jt-2:jt+3,it,kt), p(jt-2:jt+3,it,kt), &
                               tmp(jt-2:jt+3,it,kt), Normal_y)
      endif
    elseif (2 <= j .and. j <= ny-2) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(4)
        sig = dble(sigma(i,j-1:j+2,k))
        F(:,i-1,j,k-1) = KEEP_IGR4(rho(jt-1:jt+2,it,kt), u(jt-1:jt+2,it,kt), &
                                     v(jt-1:jt+2,it,kt), w(jt-1:jt+2,it,kt), &
                                    vv(jt-1:jt+2,it,kt), p(jt-1:jt+2,it,kt), &
                                   tmp(jt-1:jt+2,it,kt), sig, Normal_y)
        end block
      else
        F(:,i-1,j,k-1) = KEEP4(rho(jt-1:jt+2,it,kt), u(jt-1:jt+2,it,kt), &
                                 v(jt-1:jt+2,it,kt), w(jt-1:jt+2,it,kt), &
                                vv(jt-1:jt+2,it,kt), p(jt-1:jt+2,it,kt), &
                               tmp(jt-1:jt+2,it,kt), Normal_y)
      endif
    else
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(2)
        sig = dble(sigma(i,j:j+1,k))
        F(:,i-1,j,k-1) = KEEP_IGR2(rho(jt:jt+1,it,kt), u(jt:jt+1,it,kt), &
                                     v(jt:jt+1,it,kt), w(jt:jt+1,it,kt), &
                                    vv(jt:jt+1,it,kt), p(jt:jt+1,it,kt), &
                                   tmp(jt:jt+1,it,kt), sig, Normal_y)
        end block
      else
        F(:,i-1,j,k-1) = KEEP2(rho(jt:jt+1,it,kt), u(jt:jt+1,it,kt), &
                                 v(jt:jt+1,it,kt), w(jt:jt+1,it,kt), &
                                vv(jt:jt+1,it,kt), p(jt:jt+1,it,kt), &
                               tmp(jt:jt+1,it,kt), Normal_y)
      endif
    endif
    end associate
  end subroutine calc_keep_y6


  attributes(global) subroutine calc_keep_z6(id_accuracy, nx, ny, nz, Q, T, G, sigma)
    use mod_constant, only : Normal_z
    integer(kind=8), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: G(5,nx-2,ny-2,nz-1)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt, kk, k_base
    real(8), dimension(-1:threadsG%z+3,threadsG%y,threadsG%x), shared :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    do kk = kt-2, threadsG%z+3, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(kk,jt,it) = Q(1,i,j,k)
          u(kk,jt,it) = Q(2,i,j,k)
          v(kk,jt,it) = Q(3,i,j,k)
          w(kk,jt,it) = Q(4,i,j,k)
          p(kk,jt,it) = Q(5,i,j,k)
        tmp(kk,jt,it) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    associate(ww => w)
    if (3 <= k .and. k <= nz-3) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(6)
        sig = dble(sigma(i,j,k-2:k+3))
        G(:,i-1,j-1,k) = KEEP_IGR6(rho(kt-2:kt+3,jt,it), u(kt-2:kt+3,jt,it), &
                                     v(kt-2:kt+3,jt,it), w(kt-2:kt+3,jt,it), &
                                    ww(kt-2:kt+3,jt,it), p(kt-2:kt+3,jt,it), &
                                   tmp(kt-2:kt+3,jt,it), sig, Normal_z)
        end block
      else
        G(:,i-1,j-1,k) = KEEP6(rho(kt-2:kt+3,jt,it), u(kt-2:kt+3,jt,it), &
                                 v(kt-2:kt+3,jt,it), w(kt-2:kt+3,jt,it), &
                                ww(kt-2:kt+3,jt,it), p(kt-2:kt+3,jt,it), &
                               tmp(kt-2:kt+3,jt,it), Normal_z)
      endif
    elseif (2 <= k .and. k <= nz-2) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(4)
        sig = dble(sigma(i,j,k-1:k+2))
        G(:,i-1,j-1,kt) = KEEP_IGR4(rho(kt-1:kt+2,jt,it), u(kt-1:kt+2,jt,it), &
                                      v(kt-1:kt+2,jt,it), w(kt-1:kt+2,jt,it), &
                                     ww(kt-1:kt+2,jt,it), p(kt-1:kt+2,jt,it), &
                                    tmp(kt-1:kt+2,jt,it), sig, Normal_z)
        end block
      else
        G(:,i-1,j-1,k) = KEEP4(rho(kt-1:kt+2,jt,it), u(kt-1:kt+2,jt,it), &
                                 v(kt-1:kt+2,jt,it), w(kt-1:kt+2,jt,it), &
                                ww(kt-1:kt+2,jt,it), p(kt-1:kt+2,jt,it), &
                               tmp(kt-1:kt+2,jt,it), Normal_z)
      endif
    else
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(2)
        sig = dble(sigma(i,j,k:k+1))
        G(:,i-1,j-1,k) = KEEP_IGR2(rho(kt:kt+1,jt,it), u(kt:kt+1,jt,it), &
                                     v(kt:kt+1,jt,it), w(kt:kt+1,jt,it), &
                                    ww(kt:kt+1,jt,it), p(kt:kt+1,jt,it), &
                                   tmp(kt:kt+1,jt,it), sig, Normal_z)
        end block
      else
        G(:,i-1,j-1,k) = KEEP2(rho(kt:kt+1,jt,it), u(kt:kt+1,jt,it), &
                                 v(kt:kt+1,jt,it), w(kt:kt+1,jt,it), &
                                ww(kt:kt+1,jt,it), p(kt:kt+1,jt,it), &
                               tmp(kt:kt+1,jt,it), Normal_z)
      endif
    endif
    end associate
  end subroutine calc_keep_z6


  attributes(global) subroutine calc_keep_x4(id_accuracy, nx, ny, nz, Q, T, E, sigma)
    use mod_constant, only : Normal_x
    integer(kind=4), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: E(5,nx-1,ny-2,nz-2)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt, ii, i_base
    real(8), dimension(0:threadsE%x+2,threadsE%y,threadsE%z), shared :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    do ii = it-1, threadsE%x+2, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(ii,jt,kt) = Q(1,i,j,k)
          u(ii,jt,kt) = Q(2,i,j,k)
          v(ii,jt,kt) = Q(3,i,j,k)
          w(ii,jt,kt) = Q(4,i,j,k)
          p(ii,jt,kt) = Q(5,i,j,k)
        tmp(ii,jt,kt) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    associate(uu => u)
    if (2 <= i .and. i <= nx-2) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(4)
        sig = dble(sigma(i-1:i+2,j,k))
        E(:,i,j-1,k-1) = KEEP_IGR4(rho(it-1:it+2,jt,kt), u(it-1:it+2,jt,kt), &
                                     v(it-1:it+2,jt,kt), w(it-1:it+2,jt,kt), &
                                    uu(it-1:it+2,jt,kt), p(it-1:it+2,jt,kt), &
                                   tmp(it-1:it+2,jt,kt), sig, Normal_x)
        end block
      else
        E(:,i,j-1,k-1) = KEEP4(rho(it-1:it+2,jt,kt), u(it-1:it+2,jt,kt), &
                                 v(it-1:it+2,jt,kt), w(it-1:it+2,jt,kt), &
                                uu(it-1:it+2,jt,kt), p(it-1:it+2,jt,kt), &
                               tmp(it-1:it+2,jt,kt), Normal_x)
      endif
    else
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(2)
        sig = dble(sigma(i:i+1,j,k))
        E(:,i,j-1,k-1) = KEEP_IGR2(rho(it:it+1,jt,kt), u(it:it+1,jt,kt), &
                                     v(it:it+1,jt,kt), w(it:it+1,jt,kt), &
                                    uu(it:it+1,jt,kt), p(it:it+1,jt,kt), &
                                   tmp(it:it+1,jt,kt), sig, Normal_x)
        end block
      else
        E(:,i,j-1,k-1) = KEEP2(rho(it:it+1,jt,kt), u(it:it+1,jt,kt), &
                                 v(it:it+1,jt,kt), w(it:it+1,jt,kt), &
                                uu(it:it+1,jt,kt), p(it:it+1,jt,kt), &
                               tmp(it:it+1,jt,kt), Normal_x)
      endif
    endif
    end associate
  end subroutine calc_keep_x4


  attributes(global) subroutine calc_keep_y4(id_accuracy, nx, ny, nz, Q, T, F, sigma)
    use mod_constant, only : Normal_y
    integer(kind=4), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: F(5,nx-2,ny-1,nz-2)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt, jj, j_base
    real(8), dimension(0:threadsF%y+2,threadsF%x,threadsF%z), shared :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    do jj = jt-1, threadsF%y+2, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(jj,it,kt) = Q(1,i,j,k)
          u(jj,it,kt) = Q(2,i,j,k)
          v(jj,it,kt) = Q(3,i,j,k)
          w(jj,it,kt) = Q(4,i,j,k)
          p(jj,it,kt) = Q(5,i,j,k)
        tmp(jj,it,kt) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    associate(vv => v)
    if (2 <= j .and. j <= ny-2) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(4)
        sig = dble(sigma(i,j-1:j+2,k))
        F(:,i-1,j,k-1) = KEEP_IGR4(rho(jt-1:jt+2,it,kt), u(jt-1:jt+2,it,kt), &
                                     v(jt-1:jt+2,it,kt), w(jt-1:jt+2,it,kt), &
                                    vv(jt-1:jt+2,it,kt), p(jt-1:jt+2,it,kt), &
                                   tmp(jt-1:jt+2,it,kt), sig, Normal_y)
        end block
      else
        F(:,i-1,j,k-1) = KEEP4(rho(jt-1:jt+2,it,kt), u(jt-1:jt+2,it,kt), &
                                 v(jt-1:jt+2,it,kt), w(jt-1:jt+2,it,kt), &
                                vv(jt-1:jt+2,it,kt), p(jt-1:jt+2,it,kt), &
                               tmp(jt-1:jt+2,it,kt), Normal_y)
      endif
    else
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(2)
        sig = dble(sigma(i,j:j+1,k))
        F(:,i-1,j,k-1) = KEEP_IGR2(rho(jt:jt+1,it,kt), u(jt:jt+1,it,kt), &
                                     v(jt:jt+1,it,kt), w(jt:jt+1,it,kt), &
                                    vv(jt:jt+1,it,kt), p(jt:jt+1,it,kt), &
                                   tmp(jt:jt+1,it,kt), sig, Normal_y)
        end block
      else
        F(:,i-1,j,k-1) = KEEP2(rho(jt:jt+1,it,kt), u(jt:jt+1,it,kt), &
                                 v(jt:jt+1,it,kt), w(jt:jt+1,it,kt), &
                                vv(jt:jt+1,it,kt), p(jt:jt+1,it,kt), &
                               tmp(jt:jt+1,it,kt), Normal_y)
      endif
    endif
    end associate
  end subroutine calc_keep_y4


  attributes(global) subroutine calc_keep_z4(id_accuracy, nx, ny, nz, Q, T, G, sigma)
    use mod_constant, only : Normal_z
    integer(kind=4), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: G(5,nx-2,ny-2,nz-1)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt, kk, k_base
    real(8), dimension(0:threadsG%z+2,threadsG%y,threadsG%x), shared :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    do kk = kt-1, threadsG%z+2, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(kk,jt,it) = Q(1,i,j,k)
          u(kk,jt,it) = Q(2,i,j,k)
          v(kk,jt,it) = Q(3,i,j,k)
          w(kk,jt,it) = Q(4,i,j,k)
          p(kk,jt,it) = Q(5,i,j,k)
        tmp(kk,jt,it) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    associate(ww => w)
    if (2 <= k .and. k <= nz-2) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(4)
        sig = dble(sigma(i,j,k-1:k+2))
        G(:,i-1,j-1,k) = KEEP_IGR4(rho(kt-1:kt+2,jt,it), u(kt-1:kt+2,jt,it), &
                                     v(kt-1:kt+2,jt,it), w(kt-1:kt+2,jt,it), &
                                    ww(kt-1:kt+2,jt,it), p(kt-1:kt+2,jt,it), &
                                   tmp(kt-1:kt+2,jt,it), sig, Normal_z)
        end block
      else
        G(:,i-1,j-1,k) = KEEP4(rho(kt-1:kt+2,jt,it), u(kt-1:kt+2,jt,it), &
                                 v(kt-1:kt+2,jt,it), w(kt-1:kt+2,jt,it), &
                                ww(kt-1:kt+2,jt,it), p(kt-1:kt+2,jt,it), &
                               tmp(kt-1:kt+2,jt,it), Normal_z)
      endif
    else
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(2)
        sig = dble(sigma(i,j,k:k+1))
        G(:,i-1,j-1,k) = KEEP_IGR2(rho(kt:kt+1,jt,it), u(kt:kt+1,jt,it), &
                                     v(kt:kt+1,jt,it), w(kt:kt+1,jt,it), &
                                    ww(kt:kt+1,jt,it), p(kt:kt+1,jt,it), &
                                   tmp(kt:kt+1,jt,it), sig, Normal_z)
        end block
      else
        G(:,i-1,j-1,k) = KEEP2(rho(kt:kt+1,jt,it), u(kt:kt+1,jt,it), &
                                 v(kt:kt+1,jt,it), w(kt:kt+1,jt,it), &
                                ww(kt:kt+1,jt,it), p(kt:kt+1,jt,it), &
                               tmp(kt:kt+1,jt,it), Normal_z)
      endif
    endif
    end associate
  end subroutine calc_keep_z4


  attributes(global) subroutine calc_keep_x2(id_accuracy, nx, ny, nz, Q, T, E, sigma)
    use mod_constant, only : Normal_x
    integer(kind=2), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: E(5,nx-1,ny-2,nz-2)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt
    real(8), dimension(2), device :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho = Q(1,i:i+1,j,k)
    u   = Q(2,i:i+1,j,k)
    v   = Q(3,i:i+1,j,k)
    w   = Q(4,i:i+1,j,k)
    p   = Q(5,i:i+1,j,k)
    tmp = T(i:i+1,j,k)
    if (kind(id_igr) == 4) then
      block
      real(8), device :: sig(2)
      sig = dble(sigma(i:i+1,j,k))
      E(:,i,j-1,k-1) = KEEP_IGR2(rho, u, v, w, u, p, tmp, sig, Normal_x)
      end block
    else
      E(:,i,j-1,k-1) = KEEP2(rho, u, v, w, u, p, tmp, Normal_x)
    endif
  end subroutine calc_keep_x2


  attributes(global) subroutine calc_keep_y2(id_accuracy, nx, ny, nz, Q, T, F, sigma)
    use mod_constant, only : Normal_y
    integer(kind=2), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: F(5,nx-2,ny-1,nz-2)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt
    real(8), dimension(2), device :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho = Q(1,i,j:j+1,k)
    u   = Q(2,i,j:j+1,k)
    v   = Q(3,i,j:j+1,k)
    w   = Q(4,i,j:j+1,k)
    p   = Q(5,i,j:j+1,k)
    tmp = T(i,j:j+1,k)
    if (kind(id_igr) == 4) then
      block
      real(8), device :: sig(2)
      sig = dble(sigma(i,j:j+1,k))
      F(:,i-1,j,k-1) = KEEP_IGR2(rho, u, v, w, v, p, tmp, sig, Normal_y)
      end block
    else
      F(:,i-1,j,k-1) = KEEP2(rho, u, v, w, v, p, tmp, Normal_y)
    endif
  end subroutine calc_keep_y2


  attributes(global) subroutine calc_keep_z2(id_accuracy, nx, ny, nz, Q, T, G, sigma)
    use mod_constant, only : Normal_z
    integer(kind=2), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: G(5,nx-2,ny-2,nz-1)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt
    real(8), dimension(2), device :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho = Q(1,i,j,k:k+1)
    u   = Q(2,i,j,k:k+1)
    v   = Q(3,i,j,k:k+1)
    w   = Q(4,i,j,k:k+1)
    p   = Q(5,i,j,k:k+1)
    tmp = T(i,j,k:k+1)
    if (kind(id_igr) == 4) then
      block
      real(8), device :: sig(2)
      sig = dble(sigma(i,j,k:k+1))
      G(:,i-1,j-1,k) = KEEP_IGR2(rho, u, v, w, w, p, tmp, sig, Normal_z)
      end block
    else
      G(:,i-1,j-1,k) = KEEP2(rho, u, v, w, w, p, tmp, Normal_z)
    endif
  end subroutine calc_keep_z2
end module calc_keep_kernel

