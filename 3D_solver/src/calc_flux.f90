module calc_flux
  use mod_globals, only : id_scheme, id_accuracy, accuracy, offset, gamma
  use calc_qlr
  use calc_keep
  use calc_slau
  use calc_roe
  implicit none
  interface calc_E
    module procedure calc_E_NoMUSCL, calc_E_MUSCL, calc_E_MUSCL_4th
  end interface
  interface calc_F
    module procedure calc_F_NoMUSCL, calc_F_MUSCL, calc_F_MUSCL_4th
  end interface
  interface calc_G
    module procedure calc_G_NoMUSCL, calc_G_MUSCL, calc_G_MUSCL_4th
  end interface
contains
!NoMUSCL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(global) subroutine calc_E_NoMUSCL(id_muscl, nx, ny, nz, rho, u, v, w, p, E)
    integer(kind=2), intent(in), value                :: id_muscl
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p
    real(8), intent(out), device                      :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    integer i, j, k
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4
    real(8), dimension(4,3) :: V4
    real(8)                 :: Normal(3) = (/1.d0, 0.d0, 0.d0/)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    if (2 <= i .and. i <= nx-2) then
      rho4(:) = rho(i-1:i+2,j,k)
      p4(:) = p(i-1:i+2,j,k)
      V4(:,1) = u(i-1:i+2,j,k)
      V4(:,2) = v(i-1:i+2,j,k)
      V4(:,3) = w(i-1:i+2,j,k)
      E(i,j-offset,k-offset,:) = KEEP4(1,rho4,p4,V4,Normal)
    else
      rho2(:) = rho(i:i+1,j,k)
      p2(:) = p(i:i+1,j,k)
      V2(:,1) = u(i:i+1,j,k)
      V2(:,2) = v(i:i+1,j,k)
      V2(:,3) = w(i:i+1,j,k)
      E(i,j-offset,k-offset,:) = KEEP2(1,rho2,p2,V2,Normal)
    endif
  end subroutine calc_E_NoMUSCL

  attributes(global) subroutine calc_F_NoMUSCL(id_muscl, nx, ny, nz, rho, u, v, w, p, F)
    integer(kind=2), intent(in), value                :: id_muscl
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p
    real(8), intent(out), device                      :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    integer i, j, k
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4
    real(8), dimension(4,3) :: V4
    real(8)                 :: Normal(3) = (/0.d0, 1.d0, 0.d0/)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset
    if (2 <= j .and. j <= ny-2) then
      rho4(:) = rho(i,j-1:j+2,k)
      p4(:) = p(i,j-1:j+2,k)
      V4(:,1) = u(i,j-1:j+2,k)
      V4(:,2) = v(i,j-1:j+2,k)
      V4(:,3) = w(i,j-1:j+2,k)
      F(i-offset,j,k-offset,:) = KEEP4(2,rho4,p4,V4,Normal)
    else
      rho2(:) = rho(i,j:j+1,k)
      p2(:) = p(i,j:j+1,k)
      V2(:,1) = u(i,j:j+1,k)
      V2(:,2) = v(i,j:j+1,k)
      V2(:,3) = w(i,j:j+1,k)
      F(i-offset,j,k-offset,:) = KEEP2(2,rho2,p2,V2,Normal)
    endif
  end subroutine calc_F_NoMUSCL

  attributes(global) subroutine calc_G_NoMUSCL(id_muscl, nx, ny, nz, rho, u, v, w, p, G)
    integer(kind=2), intent(in), value                :: id_muscl
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p 
    real(8), intent(out), device                      :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    integer i, j, k
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8), dimension(4)   :: rho4, p4
    real(8), dimension(4,3) :: V4
    real(8)                 :: Normal(3) = (/0.d0, 0.d0, 1.d0/)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (2 <= k .and. k <= nz-2) then
      rho4(:) = rho(i,j,k-1:k+2)
      p4(:) = p(i,j,k-1:k+2)
      V4(:,1) = u(i,j,k-1:k+2)
      V4(:,2) = v(i,j,k-1:k+2)
      V4(:,3) = w(i,j,k-1:k+2)
      G(i-offset,j-offset,k,:) = KEEP4(3,rho4,p4,V4,Normal)
    else
      rho2(:) = rho(i,j,k:k+1)
      p2(:) = p(i,j,k:k+1)
      V2(:,1) = u(i,j,k:k+1)
      V2(:,2) = v(i,j,k:k+1)
      V2(:,3) = w(i,j,k:k+1)
      G(i-offset,j-offset,k,:) = KEEP2(3,rho2,p2,V2,Normal)
    endif
  end subroutine calc_G_NoMUSCL

!MUSCL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(global) subroutine calc_E_MUSCL(id_muscl, nx, ny, nz, rho, u, v, w, p, E)
    integer(kind=4), intent(in), value                :: id_muscl
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p
    real(8), intent(out), device                      :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    integer i, j, k
    real(8), dimension(5)   :: Q1, Q2, Q3, Q4, Ql, Qr
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8) :: Normal5(5) = (/0.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(8) :: Normal3(3) = (/1.d0, 0.d0, 0.d0/)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset

    Q2 = (/rho(i,j,k),   u(i,j,k),   v(i,j,k),   w(i,j,k),   p(i,j,k)/)
    Q3 = (/rho(i+1,j,k), u(i+1,j,k), v(i+1,j,k), w(i+1,j,k), p(i+1,j,k)/)

    if (2 <= i .and. i <= nx-2) then
      Q1 = (/rho(i-1,j,k), u(i-1,j,k), v(i-1,j,k), w(i-1,j,k), p(i-1,j,k)/)
      Q4 = (/rho(i+2,j,k), u(i+2,j,k), v(i+2,j,k), w(i+2,j,k), p(i+2,j,k)/)
      call Qlr_mid(Q1,Q2,Q3,Q4,Ql,Qr)
    elseif (i == 1) then
      Q4 = (/rho(i+2,j,k), u(i+2,j,k), v(i+2,j,k), w(i+2,j,k), p(i+2,j,k)/)
      call Qlr_left(Q2,Q3,Q4,Ql,Qr)
    else
      Q1 = (/rho(i-1,j,k), u(i-1,j,k), v(i-1,j,k), w(i-1,j,k), p(i-1,j,k)/)
      call Qlr_right(Q1,Q2,Q3,Ql,Qr)
    endif

    if (id_scheme == 2) then
      E(i,j-offset,k-offset,:) = Roe(1,Ql,Qr,Normal5)
    elseif (id_scheme == 3) then
      E(i,j-offset,k-offset,:) = SLAU(1,Ql,Qr,Normal5)
    elseif (id_scheme == 4) then
      rho2(:) = (/Ql(1), Qr(1)/)
      V2(:,1) = (/Ql(2), Qr(2)/)
      V2(:,2) = (/Ql(3), Qr(3)/)
      V2(:,3) = (/Ql(4), Qr(4)/)
      p2(:)   = (/Ql(5), Qr(5)/)
      E(i,j-offset,k-offset,:) = KEEP2(1,rho2,p2,V2,Normal3)
    endif
  end subroutine calc_E_MUSCL

  attributes(global) subroutine calc_F_MUSCL(id_muscl, nx, ny, nz, rho, u, v, w, p, F)
    integer(kind=4), intent(in), value                :: id_muscl
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p
    real(8), intent(out), device                      :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    integer i, j, k
    real(8), dimension(5)   :: Q1, Q2, Q3, Q4, Ql, Qr
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8) :: Normal5(5) = (/0.d0, 0.d0, 1.d0, 0.d0, 0.d0/)
    real(8) :: Normal3(3) = (/0.d0, 1.d0, 0.d0/)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset

    Q2 = (/rho(i,j,k),   u(i,j,k),   v(i,j,k),   w(i,j,k),   p(i,j,k)/) 
    Q3 = (/rho(i,j+1,k), u(i,j+1,k), v(i,j+1,k), w(i,j+1,k), p(i,j+1,k)/) 

    if (2 <= j .and. j <= ny-2) then
      Q1 = (/rho(i,j-1,k), u(i,j-1,k), v(i,j-1,k), w(i,j-1,k), p(i,j-1,k)/) 
      Q4 = (/rho(i,j+2,k), u(i,j+2,k), v(i,j+2,k), w(i,j+2,k), p(i,j+2,k)/) 
      call Qlr_mid(Q1,Q2,Q3,Q4,Ql,Qr)
    elseif (j == 1) then
      Q4 = (/rho(i,j+2,k), u(i,j+2,k), v(i,j+2,k), w(i,j+2,k), p(i,j+2,k)/) 
      call Qlr_left(Q2,Q3,Q4,Ql,Qr)
    else
      Q1 = (/rho(i,j-1,k), u(i,j-1,k), v(i,j-1,k), w(i,j-1,k), p(i,j-1,k)/) 
      call Qlr_right(Q1,Q2,Q3,Ql,Qr)
    endif

    if (id_scheme == 2) then
      F(i-offset,j,k-offset,:) = Roe(2,Ql,Qr,Normal5)
    elseif (id_scheme == 3) then
      F(i-offset,j,k-offset,:) = SLAU(2,Ql,Qr,Normal5)
    elseif (id_scheme == 4) then
      rho2(:) = (/Ql(1), Qr(1)/)
      V2(:,1) = (/Ql(2), Qr(2)/)
      V2(:,2) = (/Ql(3), Qr(3)/)
      V2(:,3) = (/Ql(4), Qr(4)/)
      p2(:)   = (/Ql(5), Qr(5)/)
      F(i-offset,j,k-offset,:) = KEEP2(2,rho2,p2,V2,Normal3)
    endif
  end subroutine calc_F_MUSCL
  
  attributes(global) subroutine calc_G_MUSCL(id_muscl, nx, ny, nz, rho, u, v, w, p, G)
    integer(kind=4), intent(in), value                :: id_muscl
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p
    real(8), intent(out), device                      :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    integer i, j, k
    real(8), dimension(5)   :: Q1, Q2, Q3, Q4, Ql, Qr
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8) :: Normal5(5) = (/0.d0, 0.d0, 0.d0, 1.d0, 0.d0/)
    real(8) :: Normal3(3) = (/0.d0, 0.d0, 1.d0/)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z

    Q2 = (/rho(i,j,k),   u(i,j,k),   v(i,j,k),   w(i,j,k),   p(i,j,k)/) 
    Q3 = (/rho(i,j,k+1), u(i,j,k+1), v(i,j,k+1), w(i,j,k+1), p(i,j,k+1)/) 

    if (2 <= k .and. k <= nz-2) then
      Q1 = (/rho(i,j,k-1), u(i,j,k-1), v(i,j,k-1), w(i,j,k-1), p(i,j,k-1)/) 
      Q4 = (/rho(i,j,k+2), u(i,j,k+2), v(i,j,k+2), w(i,j,k+2), p(i,j,k+2)/) 
      call Qlr_mid(Q1,Q2,Q3,Q4,Ql,Qr)
    elseif (k == 1) then
      Q4 = (/rho(i,j,k+2), u(i,j,k+2), v(i,j,k+2), w(i,j,k+2), p(i,j,k+2)/) 
      call Qlr_left(Q2,Q3,Q4,Ql,Qr)
    else
      Q1 = (/rho(i,j,k-1), u(i,j,k-1), v(i,j,k-1), w(i,j,k-1), p(i,j,k-1)/) 
      call Qlr_right(Q1,Q2,Q3,Ql,Qr)
    endif

    if (id_scheme == 2) then
      G(i-offset,j-offset,k,:) = Roe(3,Ql,Qr,Normal5)
    elseif (id_scheme == 3) then
      G(i-offset,j-offset,k,:) = SLAU(3,Ql,Qr,Normal5)
    elseif (id_scheme == 4) then
      rho2(:) = (/Ql(1), Qr(1)/)
      V2(:,1) = (/Ql(2), Qr(2)/)
      V2(:,2) = (/Ql(3), Qr(3)/)
      V2(:,3) = (/Ql(4), Qr(4)/)
      p2(:)   = (/Ql(5), Qr(5)/)
      G(i-offset,j-offset,k,:) = KEEP2(3,rho2,p2,V2,Normal3)
    endif
  end subroutine calc_G_MUSCL

!MUSCL4th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  attributes(global) subroutine calc_E_MUSCL_4th(id_muscl, nx, ny, nz, rho, u, v, w, p, E)
    integer(kind=8), intent(in), value                :: id_muscl
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p
    real(8), intent(out), device                      :: E(nx-accuracy+1,ny-accuracy,nz-accuracy,5)
    integer i, j, k
    real(8), dimension(5)   :: Q1, Q2, Q3, Q4, Q5, Q6, Ql, Qr
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8) :: Normal5(5) = (/0.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(8) :: Normal3(3) = (/1.d0, 0.d0, 0.d0/)
    real(8) :: zero(5) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset

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

    if (id_scheme == 2) then
      E(i,j-offset,k-offset,:) = Roe(1,Ql,Qr,Normal5)
    elseif (id_scheme == 3) then
      E(i,j-offset,k-offset,:) = SLAU(1,Ql,Qr,Normal5)
    elseif (id_scheme == 4) then
      rho2(:) = (/Ql(1), Qr(1)/)
      V2(:,1) = (/Ql(2), Qr(2)/)
      V2(:,2) = (/Ql(3), Qr(3)/)
      V2(:,3) = (/Ql(4), Qr(4)/)
      p2(:)   = (/Ql(5), Qr(5)/)
      E(i,j-offset,k-offset,:) = KEEP2(1,rho2,p2,V2,Normal3)
    endif
  end subroutine calc_E_MUSCL_4th

  attributes(global) subroutine calc_F_MUSCL_4th(id_muscl, nx, ny, nz, rho, u, v, w, p, F)
    integer(kind=8), intent(in), value                :: id_muscl
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p
    real(8), intent(out), device                      :: F(nx-accuracy,ny-accuracy+1,nz-accuracy,5)
    integer i, j, k
    real(8), dimension(5)   :: Q1, Q2, Q3, Q4, Q5, Q6, Ql, Qr
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8) :: Normal5(5) = (/0.d0, 0.d0, 1.d0, 0.d0, 0.d0/)
    real(8) :: Normal3(3) = (/0.d0, 1.d0, 0.d0/)
    real(8) :: zero(5) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + offset

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

    if (id_scheme == 2) then
      F(i-offset,j,k-offset,:) = Roe(2,Ql,Qr,Normal5)
    elseif (id_scheme == 3) then
      F(i-offset,j,k-offset,:) = SLAU(2,Ql,Qr,Normal5)
    elseif (id_scheme == 4) then
      rho2(:) = (/Ql(1), Qr(1)/)
      V2(:,1) = (/Ql(2), Qr(2)/)
      V2(:,2) = (/Ql(3), Qr(3)/)
      V2(:,3) = (/Ql(4), Qr(4)/)
      p2(:)   = (/Ql(5), Qr(5)/)
      F(i-offset,j,k-offset,:) = KEEP2(2,rho2,p2,V2,Normal3)
    endif
  end subroutine calc_F_MUSCL_4th
  
  attributes(global) subroutine calc_G_MUSCL_4th(id_muscl, nx, ny, nz, rho, u, v, w, p, G)
    integer(kind=8), intent(in), value                :: id_muscl
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx,ny,nz), device  :: rho, u, v, w, p
    real(8), intent(out), device                      :: G(nx-accuracy,ny-accuracy,nz-accuracy+1,5)
    integer i, j, k
    real(8), dimension(5)   :: Q1, Q2, Q3, Q4, Q5, Q6, Ql, Qr
    real(8), dimension(2)   :: rho2, p2
    real(8), dimension(2,3) :: V2
    real(8) :: Normal5(5) = (/0.d0, 0.d0, 0.d0, 1.d0, 0.d0/)
    real(8) :: Normal3(3) = (/0.d0, 0.d0, 1.d0/)
    real(8) :: zero(5) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z

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

    if (id_scheme == 2) then
      G(i-offset,j-offset,k,:) = Roe(3,Ql,Qr,Normal5)
    elseif (id_scheme == 3) then
      G(i-offset,j-offset,k,:) = SLAU(3,Ql,Qr,Normal5)
    elseif (id_scheme == 4) then
      rho2(:) = (/Ql(1), Qr(1)/)
      V2(:,1) = (/Ql(2), Qr(2)/)
      V2(:,2) = (/Ql(3), Qr(3)/)
      V2(:,3) = (/Ql(4), Qr(4)/)
      p2(:)   = (/Ql(5), Qr(5)/)
      G(i-offset,j-offset,k,:) = KEEP2(3,rho2,p2,V2,Normal3)
    endif
  end subroutine calc_G_MUSCL_4th
end module calc_flux

