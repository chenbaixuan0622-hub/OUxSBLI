module calc_hybrid_kernel
  use mod_globals, only : id_slau, gamma, threshold, threadsE, threadsF, threadsG
  use mod_constant, only : over_gamma_1, R_over_gamma_1, one_third, one_sixth, one_twelfth, two_third
  use calc_muscl
  use calc_hybrid
  implicit none
  private
  public calc_hybrid_x, calc_hybrid_y !, calc_hybrid_z
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

  !interface calc_hybrid_z
  !  module procedure calc_hybrid_z2, calc_hybrid_z4, calc_hybrid_z6
 ! end interface calc_hybrid_z
contains
  include 'calc_keep_3d.f90'
  include 'calc_slau_3d.f90'

 
    attributes(global) subroutine calc_hybrid_x6(id_accuracy, nx, ny, Q, T, sensor, E)
    use mod_constant, only : Normal_x
    integer(kind=8), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny), T(nx,ny)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)
    real(8), intent(out), device, contiguous  :: E(4,nx-1,ny-2)
    integer i, j, it, jt, ii, i_base
    real(8), dimension(-1:threadsE%x+3,threadsE%y), shared :: rho, u, v, p
    real(sp) fdx
    real(8) rho_r, u_r, v_r, p_r
    integer(kind=4) id_accuracy4
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    !kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    !k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    do ii = it-2, threadsE%x+3, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny ) then
        rho(ii,jt) = Q(i,1,j)
          u(ii,jt) = Q(i,2,j)
          v(ii,jt) = Q(i,3,j)
          !w(ii,jt) = Q(i,4,j)
          p(ii,jt) = Q(i,4,j)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j ) return
    block
      real(8) rhol, ul, vl, wl, pl
      fdx = 0.5_sp * (sensor(i,j) + sensor(i+1,j))
      if (3 <= i .and. i <= nx-3) then
        if (fdx <= threshold) then
          block
            real(8) tmp(6)
            tmp = T(i-2:i+3,j)
            associate(uu => u)
              E(:,i,j-1) = KEEP(id_accuracy, &
                                    rho(it-2:it+3,jt), u(it-2:it+3,jt), &
                                      v(it-2:it+3,jt), &
                                     uu(it-2:it+3,jt), p(it-2:it+3,jt), tmp, Normal_x)
            end associate
          end block
        else
          call delta6(fdx, rho(it-2:it+3,jt), rhol, rho_r)
          call delta6(fdx,   u(it-2:it+3,jt),   ul,   u_r)
          call delta6(fdx,   v(it-2:it+3,jt),   vl,   v_r)
          !call delta6(fdx,   w(it-2:it+3,jt,kt),   wl,   w_r)
          call delta6(fdx,   p(it-2:it+3,jt),   pl,   p_r)
        endif
      elseif (2 <= i .and. i <= nx-2) then
        if (fdx <= threshold) then
          block
            real(8) tmp(4)
            tmp = T(i-1:i+2,j)
            associate(uu => u)
              E(:,i,j-1) = KEEP(id_accuracy4, &
                                    rho(it-1:it+2,jt), u(it-1:it+2,jt), &
                                      v(it-1:it+2,jt), &
                                     uu(it-1:it+2,jt), p(it-1:it+2,jt), tmp, Normal_x)
            end associate
          end block
        else
          call delta4(fdx, rho(it-1:it+2,jt), rhol, rho_r)
          call delta4(fdx,   u(it-1:it+2,jt),   ul,   u_r)
          call delta4(fdx,   v(it-1:it+2,jt),   vl,   v_r)
          !call delta4(fdx,   w(it-1:it+2,jt,kt),   wl,   w_r)
          call delta4(fdx,   p(it-1:it+2,jt),   pl,   p_r)
        endif
      else
        if (fdx <= threshold) then
          block
            real(8) tmp(2)
            tmp = T(i:i+1,j)
            associate(uu => u)
              E(:,i,j-1) = KEEP(id_accuracy2, &
                                    rho(it:it+1,jt), u(it:it+1,jt), &
                                      v(it:it+1,jt),  &
                                     uu(it:it+1,jt), p(it:it+1,jt), tmp, Normal_x)
            end associate
          end block
        else
          rhol = rho(it,jt); rho_r = rho(it+1,jt)
            ul =   u(it,jt);   u_r =   u(it+1,jt)
            vl =   v(it,jt);   v_r =   v(it+1,jt)
            !wl =   w(it,jt,kt);   w_r =   w(it+1,jt,kt)
            pl =   p(it,jt);   p_r =   p(it+1,jt)
        endif
      endif
      call syncthreads()
      rho(it,jt) = rhol
        u(it,jt) = ul
        v(it,jt) = vl
        !w(it,jt,kt) = wl
        p(it,jt) = pl
    end block
    if (fdx > threshold) then
      associate(un1 => u(it,jt), un2 => u_r)
        call SLAU(id_slau, rho(it,jt), rho_r, u(it,jt), u_r, v(it,jt), v_r, &
                   un1, un2, p(it,jt), p_r, Normal_x, 1.0_sp, &
                  E(1,i,j-1), E(2,i,j-1), E(3,i,j-1), E(4,i,j-1))
      end associate
    endif
  end subroutine calc_hybrid_x6

    

  attributes(global) subroutine calc_hybrid_y6(id_accuracy, nx, ny,  Q, T, sensor, F)
    use mod_constant, only : Normal_y
    integer(kind=8), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny), T(nx,ny)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)
    real(8), intent(out), device, contiguous  :: F(4,nx-2,ny-1)
    integer i, j, it, jt, jj, j_base
    real(8), dimension(-1:threadsF%y+3,threadsF%x), shared :: rho, u, v, p
    real(sp) fdy
    real(8) rho_r, u_r, v_r, p_r
    integer(kind=4) id_accuracy4
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    !kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    !k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    do jj = jt-2, threadsF%y+3, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny ) then
        rho(jj,it) = Q(i,1,j)
          u(jj,it) = Q(i,2,j)
          v(jj,it) = Q(i,3,j)
          !w(jj,it) = Q(i,4,j)
          p(jj,it) = Q(i,4,j)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j) return
    block
      real(8) rhol, ul, vl, pl
      fdy = 0.5_sp * (sensor(i,j) + sensor(i,j+1))
      if (3 <= j .and. j <= ny-3) then
        if (fdy <= threshold) then
          block
            real(8) tmp(6)
            tmp = T(i,j-2:j+3)
            associate(vv => v)
              F(:,i-1,j) = KEEP(id_accuracy, &
                                    rho(jt-2:jt+3,it), u(jt-2:jt+3,it), &
                                      v(jt-2:jt+3,it),  &
                                     vv(jt-2:jt+3,it), p(jt-2:jt+3,it), tmp, Normal_y)
            end associate
          end block
        else
          call delta6(fdy, rho(jt-2:jt+3,it), rhol, rho_r)
          call delta6(fdy,   u(jt-2:jt+3,it),   ul,   u_r)
          call delta6(fdy,   v(jt-2:jt+3,it),   vl,   v_r)
          !call delta6(fdy,   w(jt-2:jt+3,it),   wl,   w_r)
          call delta6(fdy,   p(jt-2:jt+3,it),   pl,   p_r)
        endif
      elseif (2 <= j .and. j <= ny-2) then
        if (fdy <= threshold) then
          block
            real(8) tmp(4)
            tmp = T(i,j-1:j+2)
            associate(vv => v)
              F(:,i-1,j) = KEEP(id_accuracy4, &
                                    rho(jt-1:jt+2,it), u(jt-1:jt+2,it), &
                                      v(jt-1:jt+2,it),  &
                                     vv(jt-1:jt+2,it), p(jt-1:jt+2,it), tmp, Normal_y)
            end associate
          end block
        else
          call delta4(fdy, rho(jt-1:jt+2,it), rhol, rho_r)
          call delta4(fdy,   u(jt-1:jt+2,it),   ul,   u_r)
          call delta4(fdy,   v(jt-1:jt+2,it),   vl,   v_r)
         ! call delta4(fdy,   w(jt-1:jt+2,it,kt),   wl,   w_r)
          call delta4(fdy,   p(jt-1:jt+2,it),   pl,   p_r)
        endif
      else
        if (fdy <= threshold) then
          block
            real(8) tmp(2)
            tmp = T(i,j:j+1)
            associate(vv => v)
              F(:,i-1,j) = KEEP(id_accuracy2, &
                                    rho(jt:jt+1,it), u(jt:jt+1,it), &
                                      v(jt:jt+1,it),  &
                                     vv(jt:jt+1,it), p(jt:jt+1,it), tmp, Normal_y)
            end associate
          end block
        else
          rhol = rho(jt,it); rho_r = rho(jt+1,it)
            ul =   u(jt,it);   u_r =   u(jt+1,it)
            vl =   v(jt,it);   v_r =   v(jt+1,it)
            !wl =   w(jt,it,kt);   w_r =   w(jt+1,it,kt)
            pl =   p(jt,it);   p_r =   p(jt+1,it)
        endif
      endif
      call syncthreads()
      rho(jt,it) = rhol
        u(jt,it) = ul
        v(jt,it) = vl
        !w(jt,it,kt) = wl
        p(jt,it) = pl
    end block
    if (fdy > threshold) then
      associate(un1 => v(jt,it), un2 => v_r)
        call SLAU(id_slau, rho(jt,it), rho_r, u(jt,it), u_r, v(jt,it), v_r, &
                   un1, un2, p(jt,it), p_r, Normal_y, 1.0_sp, &
                  F(1,i-1,j), F(2,i-1,j), F(3,i-1,j), F(4,i-1,j))
      end associate
    endif
  end subroutine calc_hybrid_y6




  attributes(global) subroutine calc_hybrid_x4(id_accuracy, nx, ny,  Q, T, sensor, E)
    use mod_constant, only : Normal_x
    integer(kind=4), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny), T(nx,ny)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)
    real(8), intent(out), device, contiguous  :: E(4,nx-1,ny-2)
    integer i, j,  it, jt,  ii, i_base
    real(8), dimension(0:threadsE%x+2,threadsE%y), shared :: rho, u, v, p
    real(sp) fdx
    real(8) rho_r, u_r, v_r, p_r
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    !kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    !k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    do ii = it-1, threadsE%x+2, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny ) then
        rho(ii,jt) = Q(i,1,j)
          u(ii,jt) = Q(i,2,j)
          v(ii,jt) = Q(i,3,j)
        !  w(ii,jt,kt) = Q(i,4,j,k)
          p(ii,jt) = Q(i,4,j)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j ) return
    block
      real(8) rhol, ul, vl, pl
      fdx = 0.5_sp * (sensor(i,j) + sensor(i+1,j))
      if (2 <= i .and. i <= nx-2) then
        if (fdx <= threshold) then
          block
            real(8) tmp(4)
            tmp = T(i-1:i+2,j)
            associate(uu => u)
              E(:,i,j-1) = KEEP(id_accuracy, &
                                    rho(it-1:it+2,jt), u(it-1:it+2,jt), &
                                      v(it-1:it+2,jt), &
                                     uu(it-1:it+2,jt), p(it-1:it+2,jt), tmp, Normal_x)
            end associate
          end block
        else
          call delta4(fdx, rho(it-1:it+2,jt), rhol, rho_r)
          call delta4(fdx,   u(it-1:it+2,jt),   ul,   u_r)
          call delta4(fdx,   v(it-1:it+2,jt),   vl,   v_r)
         ! call delta4(fdx,   w(it-1:it+2,jt,kt),   wl,   w_r)
          call delta4(fdx,   p(it-1:it+2,jt),   pl,   p_r)
        endif
      else
        if (fdx <= threshold) then
          block
            real(8) tmp(2)
            tmp = T(i:i+1,j)
            associate(uu => u)
              E(:,i,j-1) = KEEP(id_accuracy2, &
                                    rho(it:it+1,jt), u(it:it+1,jt), &
                                      v(it:it+1,jt),  &
                                     uu(it:it+1,jt), p(it:it+1,jt), tmp, Normal_x)
            end associate
          end block
        else
          rhol = rho(it,jt); rho_r = rho(it+1,jt)
            ul =   u(it,jt);   u_r =   u(it+1,jt)
            vl =   v(it,jt);   v_r =   v(it+1,jt)
           ! wl =   w(it,jt);   w_r =   w(it+1,jt,kt)
            pl =   p(it,jt);   p_r =   p(it+1,jt)
        endif
      endif
      call syncthreads()
      rho(it,jt) = rhol
        u(it,jt) = ul
        v(it,jt) = vl
       !w(it,jt,kt) = wl
        p(it,jt) = pl
    end block
    if (fdx > threshold) then
      associate(un1 => u(it,jt), un2 => u_r)
        call SLAU(id_slau, rho(it,jt), rho_r, u(it,jt), u_r, v(it,jt), v_r, &
                   un1, un2, p(it,jt), p_r, Normal_x, 1.0_sp, &
                  E(1,i,j-1), E(2,i,j-1), E(3,i,j-1), E(4,i,j-1))
      end associate
    endif
  end subroutine calc_hybrid_x4


  attributes(global) subroutine calc_hybrid_y4(id_accuracy, nx, ny, Q, T, sensor, F)
    use mod_constant, only : Normal_y
    integer(kind=4), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny), T(nx,ny)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)
    real(8), intent(out), device, contiguous  :: F(4,nx-2,ny-1)
    integer i, j, it, jt, jj, j_base
    real(8), dimension(0:threadsF%y+2,threadsF%x), shared :: rho, u, v, p
    real(sp) fdy
    real(8) rho_r, u_r, v_r, p_r
    integer(kind=2) id_accuracy2
    it = threadIdx%x
    jt = threadIdx%y
    !kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
   !k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    do jj = jt-1, threadsF%y+2, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny ) then
        rho(jj,it) = Q(i,1,j)
          u(jj,it) = Q(i,2,j)
          v(jj,it) = Q(i,3,j)
         ! w(jj,it) = Q(i,4,j)
          p(jj,it) = Q(i,4,j)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j ) return
    block
      real(8) rhol, ul, vl,  pl
      fdy = 0.5_sp * (sensor(i,j) + sensor(i,j+1))
      if (2 <= j .and. j <= ny-2) then
        if (fdy <= threshold) then
          block
            real(8) tmp(4)
            tmp = T(i,j-1:j+2)
            associate(vv => v)
              F(:,i-1,j) = KEEP(id_accuracy, &
                                    rho(jt-1:jt+2,it), u(jt-1:jt+2,it), &
                                      v(jt-1:jt+2,it),  &
                                     vv(jt-1:jt+2,it), p(jt-1:jt+2,it), tmp, Normal_y)
            end associate
          end block
        else
          call delta4(fdy, rho(jt-1:jt+2,it), rhol, rho_r)
          call delta4(fdy,   u(jt-1:jt+2,it),   ul,   u_r)
          call delta4(fdy,   v(jt-1:jt+2,it),   vl,   v_r)
          !call delta4(fdy,   w(jt-1:jt+2,it,kt),   wl,   w_r)
          call delta4(fdy,   p(jt-1:jt+2,it),   pl,   p_r)
        endif
      else
        if (fdy <= threshold) then
          block
            real(8) tmp(2)
            tmp = T(i,j:j+1)
            associate(vv => v)
              F(:,i-1,j) = KEEP(id_accuracy2, &
                                    rho(jt:jt+1,it), u(jt:jt+1,it), &
                                      v(jt:jt+1,it),  &
                                     vv(jt:jt+1,it), p(jt:jt+1,it), tmp, Normal_y)
            end associate
          end block
        else
          rhol = rho(jt,it); rho_r = rho(jt+1,it)
            ul =   u(jt,it);   u_r =   u(jt+1,it)
            vl =   v(jt,it);   v_r =   v(jt+1,it)
           ! wl =   w(jt,it);   w_r =   w(jt+1,it,kt)
            pl =   p(jt,it);   p_r =   p(jt+1,it)
        endif
      endif
      call syncthreads()
      rho(jt,it) = rhol
        u(jt,it) = ul
        v(jt,it) = vl
       ! w(jt,it,kt) = wl
        p(jt,it) = pl
    end block
    if (fdy > threshold) then
      associate(un1 => v(jt,it), un2 => v_r)
        call SLAU(id_slau, rho(jt,it), rho_r, u(jt,it), u_r, v(jt,it), v_r, &
                   un1, un2, p(jt,it), p_r, Normal_y, 1.0_sp, &
                  F(1,i-1,j), F(2,i-1,j), F(3,i-1,j), F(4,i-1,j))
      end associate
    endif
  end subroutine calc_hybrid_y4



  attributes(global) subroutine calc_hybrid_x2(id_accuracy, nx, ny, Q, T, sensor, E)
    use mod_constant, only : Normal_x
    integer(kind=2), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny), T(nx,ny)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)
    real(8), intent(out), device, contiguous  :: E(4,nx-1,ny-2)
    integer i, j, it, jt, ii, i_base
    real(8), dimension(threadsE%x+1,threadsE%y), shared :: rho, u, v, p
    real(sp) fdx
    it = threadIdx%x
    jt = threadIdx%y
    !kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    !k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    do ii = it, threadsE%x+1, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny ) then
        rho(ii,jt) = Q(i,1,j)
          u(ii,jt) = Q(i,2,j)
          v(ii,jt) = Q(i,3,j)
          !w(ii,jt) = Q(i,4,j)
          p(ii,jt) = Q(i,5,j)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j ) return
    fdx = 0.5_sp * (sensor(i,j) + sensor(i+1,j))
    if (fdx <= threshold) then
      block
        real(8) tmp(2)
        tmp = T(i:i+1,j)
        associate(uu => u)
          E(:,i,j-1) = KEEP(id_accuracy, &
                                rho(it:it+1,jt), u(it:it+1,jt), &
                                  v(it:it+1,jt),&
                                 uu(it:it+1,jt), p(it:it+1,jt), tmp, Normal_x)
        end associate
      end block
    else
      associate(un1 => u(it,jt), un2 => u(it+1,jt))
        call SLAU(id_slau, rho(it,jt), rho(it+1,jt), u(it,jt), u(it+1,jt), v(it,jt), v(it+1,jt), &
                  un1, un2, p(it,jt), p(it+1,jt), Normal_x, 1.0_sp, &
                  E(1,i,j-1), E(2,i,j-1), E(3,i,j-1), E(4,i,j-1))
      end associate
    endif
  end subroutine calc_hybrid_x2


  attributes(global) subroutine calc_hybrid_y2(id_accuracy, nx, ny,  Q, T, sensor, F)
    use mod_constant, only : Normal_y
    integer(kind=2), intent(in), value        :: id_accuracy
    integer, intent(in), value                :: nx, ny
    real(8), intent(in), device, contiguous   :: Q(nx,4,ny), T(nx,ny)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny)
    real(8), intent(out), device, contiguous  :: F(4,nx-2,ny-1)
    integer i, j, it, jt, jj, j_base
    real(8), dimension(threadsF%y+1,threadsF%x), shared :: rho, u, v, p
    real(sp) fdy
    it = threadIdx%x
    jt = threadIdx%y
    !kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    !k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    do jj = jt, threadsF%y+1, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny ) then
        rho(jj,it) = Q(i,1,j)
          u(jj,it) = Q(i,2,j)
          v(jj,it) = Q(i,3,j)
          !w(jj,it,kt) = Q(i,4,j,k)
          p(jj,it) = Q(i,4,j)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j ) return
    fdy = 0.5_sp * (sensor(i,j) + sensor(i,j+1))
    if (fdy <= threshold) then
      block
        real(8) tmp(2)
        tmp = T(i,j:j+1)
        associate(vv => v)
          F(:,i-1,j) = KEEP(id_accuracy, &
                                rho(jt:jt+1,it), u(jt:jt+1,it), &
                                  v(jt:jt+1,it),  &
                                 vv(jt:jt+1,it), p(jt:jt+1,it), tmp, Normal_y)
        end associate
      end block
    else
      associate(un1 => v(jt,it), un2 => v(jt+1,it))
        call SLAU(id_slau, rho(jt,it), rho(jt+1,it), u(jt,it), u(jt+1,it), v(jt,it), v(jt+1,it), &
                   un1, un2, p(jt,it), p(jt+1,it), Normal_y, 1.0_sp, &
                  F(1,i-1,j), F(2,i-1,j), F(3,i-1,j), F(4,i-1,j))
      end associate
    endif
  end subroutine calc_hybrid_y2


 
end module calc_hybrid_kernel

