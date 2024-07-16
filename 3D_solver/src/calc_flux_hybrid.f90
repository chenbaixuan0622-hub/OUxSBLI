module calc_flux_hybrid
  use mod_globals, only : id_accuracy, accuracy, offset, gamma
  use calc_qlr
  use calc_keep
  use calc_slau
  implicit none
contains
  attributes(global) subroutine calc_E_hybrid(nx, ny, nz, rho, u, v, w, p, fd, E)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p, fd
    real(8), intent(out), device                      :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    integer i, j, k
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4
    real(8), dimension(4,3) :: V4
    real(8), device         :: Normal3(3) = (/1.d0, 0.d0, 0.d0/)
    real(8)                 :: Normal5(5) = (/0.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(8)                 :: zero(5) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    real(8), dimension(5)   :: Q1, Q2, Q3, Q4, Q5, Q6, Ql, Qr
    real(8) fdx, M, c
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    fdx = max(fd(i,j,k), fd(i+1,j,k))
    c = 0.5d0 * (sqrt(gamma * p(i,j,k) / rho(i,j,k)) + sqrt(gamma * p(i+1,j,k) / rho(i+1,j,k)))
    M = sqrt(0.5d0 * (u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2 + u(i+1,j,k)**2 + v(i+1,j,k)**2 + w(i+1,j,k)**2)) / c
    if (M < 1.d0) then
      ! KEEP
      if (2 <= i .and. i <= nx-2) then
        rho4(:) = rho(i-1:i+2,j,k)
        p4(:) = p(i-1:i+2,j,k)
        V4(:,1) = u(i-1:i+2,j,k)
        V4(:,2) = v(i-1:i+2,j,k)
        V4(:,3) = w(i-1:i+2,j,k)
        E(i,j-offset,k-offset,:) = KEEP4(1,rho4,p4,V4,Normal3)
        !E(i,j-offset,k-offset,:) = KEEPRho(1,rho4,p4,V4,Normal3,fdx)
      else
        ! calc SLAU at wall
        Q3 = (/rho(i,j,k),   u(i,j,k),   v(i,j,k),   w(i,j,k),   p(i,j,k)/) 
        Q4 = (/rho(i+1,j,k), u(i+1,j,k), v(i+1,j,k), w(i+1,j,k), p(i+1,j,k)/) 
        if (i == 1) then
          Q5 = (/rho(i+2,j,k), u(i+2,j,k), v(i+2,j,k), w(i+2,j,k), p(i+2,j,k)/)
          call Qlr_left(Q3,Q4,Q5,Ql,Qr)
        else
          Q2 = (/rho(i-1,j,k), u(i-1,j,k), v(i-1,j,k), w(i-1,j,k), p(i-1,j,k)/)
          call Qlr_right(Q2,Q3,Q4,Ql,Qr)
        endif
        E(i,j-offset,k-offset,:) = SLAU(1,Ql,Qr,Normal5,fdx)
      endif
    else
      ! SLAU
      Q3 = (/rho(i,j,k),   u(i,j,k),   v(i,j,k),   w(i,j,k),   p(i,j,k)/) 
      Q4 = (/rho(i+1,j,k), u(i+1,j,k), v(i+1,j,k), w(i+1,j,k), p(i+1,j,k)/) 
      if (3 <= i .and. i <= nx-3) then
        ! 4th-order MUSCL
        Q1 = (/rho(i-2,j,k), u(i-2,j,k), v(i-2,j,k), w(i-2,j,k), p(i-2,j,k)/)
        Q2 = (/rho(i-1,j,k), u(i-1,j,k), v(i-1,j,k), w(i-1,j,k), p(i-1,j,k)/)
        Q5 = (/rho(i+2,j,k), u(i+2,j,k), v(i+2,j,k), w(i+2,j,k), p(i+2,j,k)/)
        Q6 = (/rho(i+3,j,k), u(i+3,j,k), v(i+3,j,k), w(i+3,j,k), p(i+3,j,k)/)
        call Qlr_mid(1.d0,1.d0,Q1,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
      elseif (i == 2) then
        ! left 3rd-order MUSCL
        ! right 4th-order MUSCL
        Q2 = (/rho(i-1,j,k), u(i-1,j,k), v(i-1,j,k), w(i-1,j,k), p(i-1,j,k)/)
        Q5 = (/rho(i+2,j,k), u(i+2,j,k), v(i+2,j,k), w(i+2,j,k), p(i+2,j,k)/)
        Q6 = (/rho(i+3,j,k), u(i+3,j,k), v(i+3,j,k), w(i+3,j,k), p(i+3,j,k)/)
        call Qlr_mid(0.d0,1.d0,zero,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
      elseif (i == nx-2) then
        ! left 4th-order MUSCL
        ! right 3rd-order MUSCL
        Q1 = (/rho(i-2,j,k), u(i-2,j,k), v(i-2,j,k), w(i-2,j,k), p(i-2,j,k)/)
        Q2 = (/rho(i-1,j,k), u(i-1,j,k), v(i-1,j,k), w(i-1,j,k), p(i-1,j,k)/)
        Q5 = (/rho(i+2,j,k), u(i+2,j,k), v(i+2,j,k), w(i+2,j,k), p(i+2,j,k)/)
        call Qlr_mid(1.d0,0.d0,Q1,Q2,Q3,Q4,Q5,zero,Ql,Qr)
      elseif (i == 1) then
        Q5 = (/rho(i+2,j,k), u(i+2,j,k), v(i+2,j,k), w(i+2,j,k), p(i+2,j,k)/)
        call Qlr_left(Q3,Q4,Q5,Ql,Qr)
      else
        Q2 = (/rho(i-1,j,k), u(i-1,j,k), v(i-1,j,k), w(i-1,j,k), p(i-1,j,k)/)
        call Qlr_right(Q2,Q3,Q4,Ql,Qr)
      endif
      !E(i,j-offset,k-offset,:) = simpleSLAU(1,Ql,Qr,Normal5,fdx)
      E(i,j-offset,k-offset,:) = SLAU(1,Ql,Qr,Normal5,fdx)
    endif
  end subroutine calc_E_hybrid

  attributes(global) subroutine calc_F_hybrid(nx, ny, nz, rho, u, v, w, p, fd, F)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p, fd
    real(8), intent(out), device                      :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    integer i, j, k
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4
    real(8), dimension(4,3) :: V4
    real(8), device         :: Normal3(3) = (/0.d0, 1.d0, 0.d0/)
    real(8)                 :: Normal5(5) = (/0.d0, 0.d0, 1.d0, 0.d0, 0.d0/)
    real(8)                 :: zero(5) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    real(8), dimension(5)   :: Q1, Q2, Q3, Q4, Q5, Q6, Ql, Qr
    real(8) fdy, M, c
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    fdy = max(fd(i,j,k), fd(i,j+1,k))
    c = 0.5d0 * (sqrt(gamma * p(i,j,k) / rho(i,j,k)) + sqrt(gamma * p(i,j+1,k) / rho(i,j+1,k)))
    M = sqrt(0.5d0 * (u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2 + u(i,j+1,k)**2 + v(i,j+1,k)**2 + w(i,j+1,k)**2)) / c
    if (M < 1.d0) then
      ! KEEP
      if (2 <= j .and. j <= ny-2) then
        rho4(:) = rho(i,j-1:j+2,k)
        p4(:) = p(i,j-1:j+2,k)
        V4(:,1) = u(i,j-1:j+2,k)
        V4(:,2) = v(i,j-1:j+2,k)
        V4(:,3) = w(i,j-1:j+2,k)
        F(i-offset,j,k-offset,:) = KEEP4(2,rho4,p4,V4,Normal3)
        !F(i-offset,j,k-offset,:) = KEEPRho(2,rho4,p4,V4,Normal3,fdy)
      else
        ! calc SLAU at wall
        Q3 = (/rho(i,j,k),   u(i,j,k),   v(i,j,k),   w(i,j,k),   p(i,j,k)/) 
        Q4 = (/rho(i,j+1,k), u(i,j+1,k), v(i,j+1,k), w(i,j+1,k), p(i,j+1,k)/) 
        if (j == 1) then
          Q5 = (/rho(i,j+2,k), u(i,j+2,k), v(i,j+2,k), w(i,j+2,k), p(i,j+2,k)/)
          call Qlr_left(Q3,Q4,Q5,Ql,Qr)
        else
          Q2 = (/rho(i,j-1,k), u(i,j-1,k), v(i,j-1,k), w(i,j-1,k), p(i,j-1,k)/)
          call Qlr_right(Q2,Q3,Q4,Ql,Qr)
        endif
        F(i-offset,j,k-offset,:) = SLAU(2,Ql,Qr,Normal5,fdy)
      endif
    else
      ! SLAU
      Q3 = (/rho(i,j,k),   u(i,j,k),   v(i,j,k),   w(i,j,k),   p(i,j,k)/) 
      Q4 = (/rho(i,j+1,k), u(i,j+1,k), v(i,j+1,k), w(i,j+1,k), p(i,j+1,k)/) 
      if (3 <= j .and. j <= ny-3) then
        ! 4th-order MUSCL
        Q1 = (/rho(i,j-2,k), u(i,j-2,k), v(i,j-2,k), w(i,j-2,k), p(i,j-2,k)/)
        Q2 = (/rho(i,j-1,k), u(i,j-1,k), v(i,j-1,k), w(i,j-1,k), p(i,j-1,k)/)
        Q5 = (/rho(i,j+2,k), u(i,j+2,k), v(i,j+2,k), w(i,j+2,k), p(i,j+2,k)/)
        Q6 = (/rho(i,j+3,k), u(i,j+3,k), v(i,j+3,k), w(i,j+3,k), p(i,j+3,k)/)
        call Qlr_mid(1.d0,1.d0,Q1,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
      elseif (j == 2) then
        ! left 3rd-order MUSCL
        ! right 4th-order MUSCL
        Q2 = (/rho(i,j-1,k), u(i,j-1,k), v(i,j-1,k), w(i,j-1,k), p(i,j-1,k)/)
        Q5 = (/rho(i,j+2,k), u(i,j+2,k), v(i,j+2,k), w(i,j+2,k), p(i,j+2,k)/)
        Q6 = (/rho(i,j+3,k), u(i,j+3,k), v(i,j+3,k), w(i,j+3,k), p(i,j+3,k)/)
        call Qlr_mid(0.d0,1.d0,zero,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
      elseif (j == ny-2) then
        ! left 4th-order MUSCL
        ! right 3rd-order MUSCL
        Q1 = (/rho(i,j-2,k), u(i,j-2,k), v(i,j-2,k), w(i,j-2,k), p(i,j-2,k)/)
        Q2 = (/rho(i,j-1,k), u(i,j-1,k), v(i,j-1,k), w(i,j-1,k), p(i,j-1,k)/)
        Q5 = (/rho(i,j+2,k), u(i,j+2,k), v(i,j+2,k), w(i,j+2,k), p(i,j+2,k)/)
        call Qlr_mid(1.d0,0.d0,Q1,Q2,Q3,Q4,Q5,zero,Ql,Qr)
      elseif (j == 1) then
        Q5 = (/rho(i,j+2,k), u(i,j+2,k), v(i,j+2,k), w(i,j+2,k), p(i,j+2,k)/)
        call Qlr_left(Q3,Q4,Q5,Ql,Qr)
      else
        Q2 = (/rho(i,j-1,k), u(i,j-1,k), v(i,j-1,k), w(i,j-1,k), p(i,j-1,k)/)
        call Qlr_right(Q2,Q3,Q4,Ql,Qr)
      endif
      F(i-offset,j,k-offset,:) = SLAU(2,Ql,Qr,Normal5,fdy)
    endif
  end subroutine calc_F_hybrid

  attributes(global) subroutine calc_G_hybrid(nx, ny, nz, rho, u, v, w, p, fd, G)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p, fd
    real(8), intent(out), device                      :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    integer i, j, k
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4
    real(8), dimension(4,3) :: V4
    real(8), device         :: Normal3(3) = (/0.d0, 0.d0, 1.d0/)
    real(8)                 :: Normal5(5) = (/0.d0, 0.d0, 0.d0, 1.d0, 0.d0/)
    real(8)                 :: zero(5) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    real(8), dimension(5)   :: Q1, Q2, Q3, Q4, Q5, Q6, Ql, Qr
    real(8) fdz, M, c
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    fdz = max(fd(i,j,k), fd(i,j,k+1))
    c = 0.5d0 * (sqrt(gamma * p(i,j,k) / rho(i,j,k)) + sqrt(gamma * p(i,j,k+1) / rho(i,j,k+1)))
    M = sqrt(0.5d0 * (u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2 + u(i,j,k+1)**2 + v(i,j,k+1)**2 + w(i,j,k+1)**2)) / c
    if (M < 1.d0) then
      ! KEEP
      if (2 <= k .and. k <= nz-2) then
        rho4(:) = rho(i,j,k-1:k+2)
        p4(:) = p(i,j,k-1:k+2)
        V4(:,1) = u(i,j,k-1:k+2)
        V4(:,2) = v(i,j,k-1:k+2)
        V4(:,3) = w(i,j,k-1:k+2)
        G(i-offset,j-offset,k,:) = KEEP4(3,rho4,p4,V4,Normal3)
        !G(i-offset,j-offset,k,:) = KEEPRho(3,rho4,p4,V4,Normal3,fdz)
      else
        ! calc SLAU at wall
        Q3 = (/rho(i,j,k),   u(i,j,k),   v(i,j,k),   w(i,j,k),   p(i,j,k)/) 
        Q4 = (/rho(i,j,k+1), u(i,j,k+1), v(i,j,k+1), w(i,j,k+1), p(i,j,k+1)/) 
        if (k == 1) then
          Q5 = (/rho(i,j,k+2), u(i,j,k+2), v(i,j,k+2), w(i,j,k+2), p(i,j,k+2)/)
          call Qlr_left(Q3,Q4,Q5,Ql,Qr)
        else
          Q2 = (/rho(i,j,k-1), u(i,j,k-1), v(i,j,k-1), w(i,j,k-1), p(i,j,k-1)/)
          call Qlr_right(Q2,Q3,Q4,Ql,Qr)
        endif
        G(i-offset,j-offset,k,:) = SLAU(3,Ql,Qr,Normal5,fdz)
      endif
    else
      ! SLAU
      Q3 = (/rho(i,j,k),   u(i,j,k),   v(i,j,k),   w(i,j,k),   p(i,j,k)/) 
      Q4 = (/rho(i,j,k+1), u(i,j,k+1), v(i,j,k+1), w(i,j,k+1), p(i,j,k+1)/) 
      if (3 <= k .and. k <= nz-3) then
        ! 4th-order MUSCL
        Q1 = (/rho(i,j,k-2), u(i,j,k-2), v(i,j,k-2), w(i,j,k-2), p(i,j,k-2)/)
        Q2 = (/rho(i,j,k-1), u(i,j,k-1), v(i,j,k-1), w(i,j,k-1), p(i,j,k-1)/)
        Q5 = (/rho(i,j,k+2), u(i,j,k+2), v(i,j,k+2), w(i,j,k+2), p(i,j,k+2)/)
        Q6 = (/rho(i,j,k+3), u(i,j,k+3), v(i,j,k+3), w(i,j,k+3), p(i,j,k+3)/)
        call Qlr_mid(1.d0,1.d0,Q1,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
      elseif (k == 2) then
        ! left 3rd-order MUSCL
        ! right 4th-order MUSCL
        Q2 = (/rho(i,j,k-1), u(i,j,k-1), v(i,j,k-1), w(i,j,k-1), p(i,j,k-1)/)
        Q5 = (/rho(i,j,k+2), u(i,j,k+2), v(i,j,k+2), w(i,j,k+2), p(i,j,k+2)/)
        Q6 = (/rho(i,j,k+3), u(i,j,k+3), v(i,j,k+3), w(i,j,k+3), p(i,j,k+3)/)
        call Qlr_mid(0.d0,1.d0,zero,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
      elseif (k == nz-2) then
        ! left 4th-order MUSCL
        ! right 3rd-order MUSCL
        Q1 = (/rho(i,j,k-2), u(i,j,k-2), v(i,j,k-2), w(i,j,k-2), p(i,j,k-2)/)
        Q2 = (/rho(i,j,k-1), u(i,j,k-1), v(i,j,k-1), w(i,j,k-1), p(i,j,k-1)/)
        Q5 = (/rho(i,j,k+2), u(i,j,k+2), v(i,j,k+2), w(i,j,k+2), p(i,j,k+2)/)
        call Qlr_mid(1.d0,0.d0,Q1,Q2,Q3,Q4,Q5,zero,Ql,Qr)
      elseif (k == 1) then
        Q5 = (/rho(i,j,k+2), u(i,j,k+2), v(i,j,k+2), w(i,j,k+2), p(i,j,k+2)/)
        call Qlr_left(Q3,Q4,Q5,Ql,Qr)
      else
        Q2 = (/rho(i,j,k-1), u(i,j,k-1), v(i,j,k-1), w(i,j,k-1), p(i,j,k-1)/)
        call Qlr_right(Q2,Q3,Q4,Ql,Qr)
      endif
      G(i-offset,j-offset,k,:) = SLAU(3,Ql,Qr,Normal5,fdz)
    endif
  end subroutine calc_G_hybrid

  !weighted hybrid scheme!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(global) subroutine calc_E_weight(nx, ny, nz, rho, u, v, w, p, fd, E)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p, fd
    real(8), intent(out), device                      :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    integer i, j, k
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4
    real(8), dimension(4,3) :: V4
    real(8), device         :: Normal3(3) = (/1.d0, 0.d0, 0.d0/)
    real(8)                 :: Normal5(5) = (/0.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(8)                 :: zero(5) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    real(8), dimension(5)   :: Q1, Q2, Q3, Q4, Q5, Q6, Ql, Qr
    real(8) fdx, Fkeep(5), Fslau(5)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    fdx = 0.5d0 * (fd(i,j,k) + fd(i+1,j,k))
    if (fdx > 0.8) then
      fdx = 1.d0
    else
      fdx = 0.d0
    endif
    ! KEEP
    if (2 <= i .and. i <= nx-2) then
      rho4(:) = rho(i-1:i+2,j,k)
      p4(:) = p(i-1:i+2,j,k)
      V4(:,1) = u(i-1:i+2,j,k)
      V4(:,2) = v(i-1:i+2,j,k)
      V4(:,3) = w(i-1:i+2,j,k)
      Fkeep(:) = KEEP4(1,rho4,p4,V4,Normal3)
    else
      ! calc SLAU at wall
      Q3 = (/rho(i,j,k),   u(i,j,k),   v(i,j,k),   w(i,j,k),   p(i,j,k)/) 
      Q4 = (/rho(i+1,j,k), u(i+1,j,k), v(i+1,j,k), w(i+1,j,k), p(i+1,j,k)/) 
      if (i == 1) then
        Q5 = (/rho(i+2,j,k), u(i+2,j,k), v(i+2,j,k), w(i+2,j,k), p(i+2,j,k)/)
        call Qlr_left(Q3,Q4,Q5,Ql,Qr)
      else
        Q2 = (/rho(i-1,j,k), u(i-1,j,k), v(i-1,j,k), w(i-1,j,k), p(i-1,j,k)/)
        call Qlr_right(Q2,Q3,Q4,Ql,Qr)
      endif
      Fkeep(:) = SLAU(1,Ql,Qr,Normal5,fdx)
    endif
    ! SLAU
    Q3 = (/rho(i,j,k),   u(i,j,k),   v(i,j,k),   w(i,j,k),   p(i,j,k)/) 
    Q4 = (/rho(i+1,j,k), u(i+1,j,k), v(i+1,j,k), w(i+1,j,k), p(i+1,j,k)/) 
    if (3 <= i .and. i <= nx-3) then
      ! 4th-order MUSCL
      Q1 = (/rho(i-2,j,k), u(i-2,j,k), v(i-2,j,k), w(i-2,j,k), p(i-2,j,k)/)
      Q2 = (/rho(i-1,j,k), u(i-1,j,k), v(i-1,j,k), w(i-1,j,k), p(i-1,j,k)/)
      Q5 = (/rho(i+2,j,k), u(i+2,j,k), v(i+2,j,k), w(i+2,j,k), p(i+2,j,k)/)
      Q6 = (/rho(i+3,j,k), u(i+3,j,k), v(i+3,j,k), w(i+3,j,k), p(i+3,j,k)/)
      call Qlr_mid(1.d0,1.d0,Q1,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
    elseif (i == 2) then
      ! left 3rd-order MUSCL
      ! right 4th-order MUSCL
      Q2 = (/rho(i-1,j,k), u(i-1,j,k), v(i-1,j,k), w(i-1,j,k), p(i-1,j,k)/)
      Q5 = (/rho(i+2,j,k), u(i+2,j,k), v(i+2,j,k), w(i+2,j,k), p(i+2,j,k)/)
      Q6 = (/rho(i+3,j,k), u(i+3,j,k), v(i+3,j,k), w(i+3,j,k), p(i+3,j,k)/)
      call Qlr_mid(0.d0,1.d0,zero,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
    elseif (i == nx-2) then
      ! left 4th-order MUSCL
      ! right 3rd-order MUSCL
      Q1 = (/rho(i-2,j,k), u(i-2,j,k), v(i-2,j,k), w(i-2,j,k), p(i-2,j,k)/)
      Q2 = (/rho(i-1,j,k), u(i-1,j,k), v(i-1,j,k), w(i-1,j,k), p(i-1,j,k)/)
      Q5 = (/rho(i+2,j,k), u(i+2,j,k), v(i+2,j,k), w(i+2,j,k), p(i+2,j,k)/)
      call Qlr_mid(1.d0,0.d0,Q1,Q2,Q3,Q4,Q5,zero,Ql,Qr)
    elseif (i == 1) then
      Q5 = (/rho(i+2,j,k), u(i+2,j,k), v(i+2,j,k), w(i+2,j,k), p(i+2,j,k)/)
      call Qlr_left(Q3,Q4,Q5,Ql,Qr)
    else
      Q2 = (/rho(i-1,j,k), u(i-1,j,k), v(i-1,j,k), w(i-1,j,k), p(i-1,j,k)/)
      call Qlr_right(Q2,Q3,Q4,Ql,Qr)
    endif
    Fslau(:) = SLAU(1,Ql,Qr,Normal5,fdx)
    E(i,j-offset,k-offset,:) = (1.d0 - fdx**2) * Fkeep(:) + fdx**2 * Fslau(:)
  end subroutine calc_E_weight

  attributes(global) subroutine calc_F_weight(nx, ny, nz, rho, u, v, w, p, fd, F)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p, fd
    real(8), intent(out), device                      :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    integer i, j, k
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4
    real(8), dimension(4,3) :: V4
    real(8), device         :: Normal3(3) = (/0.d0, 1.d0, 0.d0/)
    real(8)                 :: Normal5(5) = (/0.d0, 0.d0, 1.d0, 0.d0, 0.d0/)
    real(8)                 :: zero(5) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    real(8), dimension(5)   :: Q1, Q2, Q3, Q4, Q5, Q6, Ql, Qr
    real(8) fdy, Fkeep(5), Fslau(5)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    fdy = 0.5d0 * (fd(i,j,k) + fd(i,j+1,k))
    if (fdy > 0.8) then
      fdy = 1.d0
    else
      fdy = 0.d0
    endif
    ! KEEP
    if (2 <= j .and. j <= ny-2) then
      rho4(:) = rho(i,j-1:j+2,k)
      p4(:) = p(i,j-1:j+2,k)
      V4(:,1) = u(i,j-1:j+2,k)
      V4(:,2) = v(i,j-1:j+2,k)
      V4(:,3) = w(i,j-1:j+2,k)
      Fkeep(:) = KEEP4(2,rho4,p4,V4,Normal3)
    else
      ! calc SLAU at wall
      Q3 = (/rho(i,j,k),   u(i,j,k),   v(i,j,k),   w(i,j,k),   p(i,j,k)/) 
      Q4 = (/rho(i,j+1,k), u(i,j+1,k), v(i,j+1,k), w(i,j+1,k), p(i,j+1,k)/) 
      if (j == 1) then
        Q5 = (/rho(i,j+2,k), u(i,j+2,k), v(i,j+2,k), w(i,j+2,k), p(i,j+2,k)/)
        call Qlr_left(Q3,Q4,Q5,Ql,Qr)
      else
        Q2 = (/rho(i,j-1,k), u(i,j-1,k), v(i,j-1,k), w(i,j-1,k), p(i,j-1,k)/)
        call Qlr_right(Q2,Q3,Q4,Ql,Qr)
      endif
      Fkeep(:) = SLAU(2,Ql,Qr,Normal5,fdy)
    endif
    ! SLAU
    Q3 = (/rho(i,j,k),   u(i,j,k),   v(i,j,k),   w(i,j,k),   p(i,j,k)/) 
    Q4 = (/rho(i,j+1,k), u(i,j+1,k), v(i,j+1,k), w(i,j+1,k), p(i,j+1,k)/) 
    if (3 <= j .and. j <= ny-3) then
      ! 4th-order MUSCL
      Q1 = (/rho(i,j-2,k), u(i,j-2,k), v(i,j-2,k), w(i,j-2,k), p(i,j-2,k)/)
      Q2 = (/rho(i,j-1,k), u(i,j-1,k), v(i,j-1,k), w(i,j-1,k), p(i,j-1,k)/)
      Q5 = (/rho(i,j+2,k), u(i,j+2,k), v(i,j+2,k), w(i,j+2,k), p(i,j+2,k)/)
      Q6 = (/rho(i,j+3,k), u(i,j+3,k), v(i,j+3,k), w(i,j+3,k), p(i,j+3,k)/)
      call Qlr_mid(1.d0,1.d0,Q1,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
    elseif (j == 2) then
      ! left 3rd-order MUSCL
      ! right 4th-order MUSCL
      Q2 = (/rho(i,j-1,k), u(i,j-1,k), v(i,j-1,k), w(i,j-1,k), p(i,j-1,k)/)
      Q5 = (/rho(i,j+2,k), u(i,j+2,k), v(i,j+2,k), w(i,j+2,k), p(i,j+2,k)/)
      Q6 = (/rho(i,j+3,k), u(i,j+3,k), v(i,j+3,k), w(i,j+3,k), p(i,j+3,k)/)
      call Qlr_mid(0.d0,1.d0,zero,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
    elseif (j == ny-2) then
      ! left 4th-order MUSCL
      ! right 3rd-order MUSCL
      Q1 = (/rho(i,j-2,k), u(i,j-2,k), v(i,j-2,k), w(i,j-2,k), p(i,j-2,k)/)
      Q2 = (/rho(i,j-1,k), u(i,j-1,k), v(i,j-1,k), w(i,j-1,k), p(i,j-1,k)/)
      Q5 = (/rho(i,j+2,k), u(i,j+2,k), v(i,j+2,k), w(i,j+2,k), p(i,j+2,k)/)
      call Qlr_mid(1.d0,0.d0,Q1,Q2,Q3,Q4,Q5,zero,Ql,Qr)
    elseif (j == 1) then
      Q5 = (/rho(i,j+2,k), u(i,j+2,k), v(i,j+2,k), w(i,j+2,k), p(i,j+2,k)/)
      call Qlr_left(Q3,Q4,Q5,Ql,Qr)
    else
      Q2 = (/rho(i,j-1,k), u(i,j-1,k), v(i,j-1,k), w(i,j-1,k), p(i,j-1,k)/)
      call Qlr_right(Q2,Q3,Q4,Ql,Qr)
    endif
    Fslau(:) = SLAU(2,Ql,Qr,Normal5,fdy)
    F(i-offset,j,k-offset,:) = (1.d0 - fdy**2) * Fkeep(:) + fdy**2 * Fslau(:)
  end subroutine calc_F_weight

  attributes(global) subroutine calc_G_weight(nx, ny, nz, rho, u, v, w, p, fd, G)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p, fd
    real(8), intent(out), device                      :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    integer i, j, k
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4
    real(8), dimension(4,3) :: V4
    real(8), device         :: Normal3(3) = (/0.d0, 0.d0, 1.d0/)
    real(8)                 :: Normal5(5) = (/0.d0, 0.d0, 0.d0, 1.d0, 0.d0/)
    real(8)                 :: zero(5) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    real(8), dimension(5)   :: Q1, Q2, Q3, Q4, Q5, Q6, Ql, Qr
    real(8) fdz, Fkeep(5), Fslau(5)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    fdz = 0.5d0 * (fd(i,j,k) + fd(i,j,k+1))
    if (fdz > 0.8) then
      fdz = 1.d0
    else
      fdz = 0.d0
    endif
    ! KEEP
    if (2 <= k .and. k <= nz-2) then
      rho4(:) = rho(i,j,k-1:k+2)
      p4(:) = p(i,j,k-1:k+2)
      V4(:,1) = u(i,j,k-1:k+2)
      V4(:,2) = v(i,j,k-1:k+2)
      V4(:,3) = w(i,j,k-1:k+2)
      Fkeep(:) = KEEP4(3,rho4,p4,V4,Normal3)
    else
      ! calc SLAU at wall
      Q3 = (/rho(i,j,k),   u(i,j,k),   v(i,j,k),   w(i,j,k),   p(i,j,k)/) 
      Q4 = (/rho(i,j,k+1), u(i,j,k+1), v(i,j,k+1), w(i,j,k+1), p(i,j,k+1)/) 
      if (k == 1) then
        Q5 = (/rho(i,j,k+2), u(i,j,k+2), v(i,j,k+2), w(i,j,k+2), p(i,j,k+2)/)
        call Qlr_left(Q3,Q4,Q5,Ql,Qr)
      else
        Q2 = (/rho(i,j,k-1), u(i,j,k-1), v(i,j,k-1), w(i,j,k-1), p(i,j,k-1)/)
        call Qlr_right(Q2,Q3,Q4,Ql,Qr)
      endif
      Fkeep(:) = SLAU(3,Ql,Qr,Normal5,fdz)
    endif
    ! SLAU
    Q3 = (/rho(i,j,k),   u(i,j,k),   v(i,j,k),   w(i,j,k),   p(i,j,k)/) 
    Q4 = (/rho(i,j,k+1), u(i,j,k+1), v(i,j,k+1), w(i,j,k+1), p(i,j,k+1)/) 
    if (3 <= k .and. k <= nz-3) then
      ! 4th-order MUSCL
      Q1 = (/rho(i,j,k-2), u(i,j,k-2), v(i,j,k-2), w(i,j,k-2), p(i,j,k-2)/)
      Q2 = (/rho(i,j,k-1), u(i,j,k-1), v(i,j,k-1), w(i,j,k-1), p(i,j,k-1)/)
      Q5 = (/rho(i,j,k+2), u(i,j,k+2), v(i,j,k+2), w(i,j,k+2), p(i,j,k+2)/)
      Q6 = (/rho(i,j,k+3), u(i,j,k+3), v(i,j,k+3), w(i,j,k+3), p(i,j,k+3)/)
      call Qlr_mid(1.d0,1.d0,Q1,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
    elseif (k == 2) then
      ! left 3rd-order MUSCL
      ! right 4th-order MUSCL
      Q2 = (/rho(i,j,k-1), u(i,j,k-1), v(i,j,k-1), w(i,j,k-1), p(i,j,k-1)/)
      Q5 = (/rho(i,j,k+2), u(i,j,k+2), v(i,j,k+2), w(i,j,k+2), p(i,j,k+2)/)
      Q6 = (/rho(i,j,k+3), u(i,j,k+3), v(i,j,k+3), w(i,j,k+3), p(i,j,k+3)/)
      call Qlr_mid(0.d0,1.d0,zero,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
    elseif (k == nz-2) then
      ! left 4th-order MUSCL
      ! right 3rd-order MUSCL
      Q1 = (/rho(i,j,k-2), u(i,j,k-2), v(i,j,k-2), w(i,j,k-2), p(i,j,k-2)/)
      Q2 = (/rho(i,j,k-1), u(i,j,k-1), v(i,j,k-1), w(i,j,k-1), p(i,j,k-1)/)
      Q5 = (/rho(i,j,k+2), u(i,j,k+2), v(i,j,k+2), w(i,j,k+2), p(i,j,k+2)/)
      call Qlr_mid(1.d0,0.d0,Q1,Q2,Q3,Q4,Q5,zero,Ql,Qr)
    elseif (k == 1) then
      Q5 = (/rho(i,j,k+2), u(i,j,k+2), v(i,j,k+2), w(i,j,k+2), p(i,j,k+2)/)
      call Qlr_left(Q3,Q4,Q5,Ql,Qr)
    else
      Q2 = (/rho(i,j,k-1), u(i,j,k-1), v(i,j,k-1), w(i,j,k-1), p(i,j,k-1)/)
      call Qlr_right(Q2,Q3,Q4,Ql,Qr)
    endif
    Fslau(:) = SLAU(3,Ql,Qr,Normal5,fdz)
    G(i-offset,j-offset,k,:) = (1.d0 - fdz**2) * Fkeep(:) + fdz**2 * Fslau(:)
  end subroutine calc_G_weight
end module calc_flux_hybrid

