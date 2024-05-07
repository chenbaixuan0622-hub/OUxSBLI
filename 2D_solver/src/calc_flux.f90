module calc_flux
  use mod_globals, only : accuracy, offset, id_accuracy, id_scheme, nx, ny
  use calc_qlr
  use calc_slau
  use calc_roe
  use calc_keep
  implicit none
  interface calc_E
    module procedure calc_E_MUSCL, calc_E_MUSCL_4th, calc_E_NoMUSCL
  end interface
  interface calc_F
    module procedure calc_F_MUSCL, calc_F_MUSCL_4th, calc_F_NoMUSCL
  end interface
contains
!NoMUSCL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(global) subroutine calc_E_NoMUSCL(id_muscl, rho, u, v, p, xix, Jacobian, E)
    integer(kind=2), intent(in), value :: id_muscl
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(in), device :: xix(nx), Jacobian(nx,ny)
    real(8), intent(out), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    integer i, j
    real(8), dimension(accuracy) :: rhos, ps
    real(8), dimension(2) :: Normal
    real(8), dimension(accuracy,2) :: Vs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    
    rhos = rho(i:i+accuracy-1,j)
    ps = p(i:i+accuracy-1,j)
    Vs(:,1) = u(i:i+accuracy-1,j)
    Vs(:,2) = v(i:i+accuracy-1,j)
    Normal(:) = (/0.5d0 * (xix(i) / Jacobian(i,j) + xix(i+1) / Jacobian(i+1,j)), 0.d0/)
    E(i,j-offset,:) = KEEP(id_accuracy,1,rhos,ps,Vs,Normal)
  end subroutine calc_E_NoMUSCL

  attributes(global) subroutine calc_F_NoMUSCL(id_muscl, rho, u, v, p, etay, Jacobian, F)
    integer(kind=2), intent(in), value :: id_muscl
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(in), device :: etay(ny), Jacobian(nx,ny)
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    integer i, j
    real(8), dimension(accuracy) :: rhos, ps
    real(8), dimension(2) :: Normal
    real(8), dimension(accuracy,2) :: Vs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    
    rhos = rho(i,j:j+accuracy-1)
    ps = p(i,j:j+accuracy-1)
    Vs(:,1) = u(i,j:j+accuracy-1)
    Vs(:,2) = v(i,j:j+accuracy-1)
    Normal(:) = (/0.d0, 0.5d0 * (etay(j) / Jacobian(i,j) + etay(j+1) / Jacobian(i,j+1))/)
    F(i-offset,j,:) = KEEP(id_accuracy,2,rhos,ps,Vs,Normal)
  end subroutine calc_F_NoMUSCL

!MUSCL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(global) subroutine calc_E_MUSCL(id_muscl, rho, u, v, p, xix, Jacobians, E)
    integer(kind=4), intent(in), value :: id_muscl
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(in), device :: xix(nx), Jacobians(nx,ny)
    real(8), intent(out), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    integer i, j
    real(8) Jacobian
    real(8), dimension(4) :: Q1, Q2, Q3, Q4, Ql, Qr, Normal
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    
    Q2 = (/rho(i,j),   u(i,j),   v(i,j),   p(i,j)/)
    Q3 = (/rho(i+1,j), u(i+1,j), v(i+1,j), p(i+1,j)/)

    if (2 <= i .and. i <= nx-2) then
      Q1 = (/rho(i-1,j), u(i-1,j), v(i-1,j), p(i-1,j)/)
      Q4 = (/rho(i+2,j), u(i+2,j), v(i+2,j), p(i+2,j)/)
      call Qlr_mid(Q1,Q2,Q3,Q4,Ql,Qr)
    elseif (i == 1) then
      Q4 = (/rho(i+2,j), u(i+2,j), v(i+2,j), p(i+2,j)/)
      call Qlr_left(Q2,Q3,Q4,Ql,Qr)
    else
      Q1 = (/rho(i-1,j), u(i-1,j), v(i-1,j), p(i-1,j)/)
      call Qlr_right(Q1,Q2,Q3,Ql,Qr)
    endif

    Jacobian = 0.5d0 * (Jacobians(i,j) + Jacobians(i+1,j))
    Normal(:) = (/0.d0, 0.5d0 * (xix(i) + xix(i+1)), 0.d0, 0.d0/)
    if (id_scheme == 2) then
      E(i,j-offset,:) = Roe(1,Ql,Qr,Normal,Jacobian)
    else
      E(i,j-offset,:) = SLAU(1,Ql,Qr,Normal,Jacobian)
    endif
  end subroutine calc_E_MUSCL

  attributes(global) subroutine calc_F_MUSCL(id_muscl, rho, u, v, p, etay, Jacobians, F)
    integer(kind=4), intent(in), value :: id_muscl
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(in), device :: etay(ny), Jacobians(nx,ny)
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    integer i, j
    real(8) Jacobian
    real(8), dimension(4) :: Q1, Q2, Q3, Q4, Ql, Qr, Normal
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y

    Q2 = (/rho(i,j),   u(i,j),   v(i,j),   p(i,j)/) 
    Q3 = (/rho(i,j+1), u(i,j+1), v(i,j+1), p(i,j+1)/) 

    if (2 <= j .and. j <= ny-2) then
      Q1 = (/rho(i,j-1), u(i,j-1), v(i,j-1), p(i,j-1)/) 
      Q4 = (/rho(i,j+2), u(i,j+2), v(i,j+2), p(i,j+2)/) 
      call Qlr_mid(Q1,Q2,Q3,Q4,Ql,Qr)
    elseif (j == 1) then
      Q4 = (/rho(i,j+2), u(i,j+2), v(i,j+2), p(i,j+2)/) 
      call Qlr_left(Q2,Q3,Q4,Ql,Qr)
    else
      Q1 = (/rho(i,j-1), u(i,j-1), v(i,j-1), p(i,j-1)/) 
      call Qlr_right(Q1,Q2,Q3,Ql,Qr)
    endif

    Jacobian = 0.5d0 * (Jacobians(i,j) + Jacobians(i,j+1))
    Normal(:) = (/0.d0, 0.d0, 0.5d0 * (etay(j) + etay(j+1)), 0.d0/)
    if (id_scheme == 2) then
      F(i-offset,j,:) = Roe(2,Ql,Qr,Normal,Jacobian)
    else
      F(i-offset,j,:) = SLAU(2,Ql,Qr,Normal,Jacobian)
    endif
  end subroutine calc_F_MUSCL

!MUSCL4th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(global) subroutine calc_E_MUSCL_4th(id_muscl, rho, u, v, p, xix, Jacobians, E)
    integer(kind=8), intent(in), value :: id_muscl
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(in), device :: xix(nx), Jacobians(nx,ny)
    real(8), intent(out), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    integer i, j
    real(8) Jacobian
    real(8), dimension(4) :: Q1, Q2, Q3, Q4, Q5, Q6, Ql, Qr, Normal
    real(8) :: zero(4) = (/0.d0, 0.d0, 0.d0, 0.d0/)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + offset
    
    Q3 = (/rho(i,j),   u(i,j),   v(i,j),   p(i,j)/)
    Q4 = (/rho(i+1,j), u(i+1,j), v(i+1,j), p(i+1,j)/)

    if (3 <= i .and. i <= nx-3) then
      ! 4th-order MUSCL
      Q1 = (/rho(i-2,j), u(i-2,j), v(i-2,j), p(i-2,j)/)
      Q2 = (/rho(i-1,j), u(i-1,j), v(i-1,j), p(i-1,j)/)
      Q5 = (/rho(i+2,j), u(i+2,j), v(i+2,j), p(i+2,j)/)
      Q6 = (/rho(i+3,j), u(i+3,j), v(i+3,j), p(i+3,j)/)
      call Qlr_mid(1.d0,1.d0,Q1,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
    elseif (i == 2) then
      ! left 3rd-order MUSCL
      ! right 4th-order MUSCL
      Q2 = (/rho(i-1,j), u(i-1,j), v(i-1,j), p(i-1,j)/)
      Q5 = (/rho(i+2,j), u(i+2,j), v(i+2,j), p(i+2,j)/)
      Q6 = (/rho(i+3,j), u(i+3,j), v(i+3,j), p(i+3,j)/)
      call Qlr_mid(0.d0,1.d0,zero,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
    elseif (i == nx-2) then
      ! left 4th-order MUSCL
      ! right 3rd-order MUSCL
      Q1 = (/rho(i-2,j), u(i-2,j), v(i-2,j), p(i-2,j)/)
      Q2 = (/rho(i-1,j), u(i-1,j), v(i-1,j), p(i-1,j)/)
      Q5 = (/rho(i+2,j), u(i+2,j), v(i+2,j), p(i+2,j)/)
      call Qlr_mid(1.d0,0.d0,Q1,Q2,Q3,Q4,Q5,zero,Ql,Qr)
    elseif (i == 1) then
      Q5 = (/rho(i+2,j), u(i+2,j), v(i+2,j), p(i+2,j)/)
      call Qlr_left(Q3,Q4,Q5,Ql,Qr)
    else
      Q2 = (/rho(i-1,j), u(i-1,j), v(i-1,j), p(i-1,j)/)
      call Qlr_right(Q2,Q3,Q4,Ql,Qr)
    endif

    Jacobian = 0.5d0 * (Jacobians(i,j) + Jacobians(i+1,j))
    Normal(:) = (/0.d0, 0.5d0 * (xix(i) + xix(i+1)), 0.d0, 0.d0/)
    if (id_scheme == 2) then
      E(i,j-offset,:) = Roe(1,Ql,Qr,Normal,Jacobian)
    else
      E(i,j-offset,:) = SLAU(1,Ql,Qr,Normal,Jacobian)
    endif
  end subroutine calc_E_MUSCL_4th

  attributes(global) subroutine calc_F_MUSCL_4th(id_muscl, rho, u, v, p, etay, Jacobians, F)
    integer(kind=8), intent(in), value :: id_muscl
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(in), device :: etay(ny), Jacobians(nx,ny)
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    integer i, j
    real(8) Jacobian
    real(8), dimension(4) :: Q1, Q2, Q3, Q4, Q5, Q6, Ql, Qr, Normal
    real(8) :: zero(4) = (/0.d0, 0.d0, 0.d0, 0.d0/)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + offset
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y

    Q3 = (/rho(i,j),   u(i,j),   v(i,j),   p(i,j)/) 
    Q4 = (/rho(i,j+1), u(i,j+1), v(i,j+1), p(i,j+1)/) 

    if (3 <= j .and. j <= ny-3) then
      ! 4th-order MUSCL
      Q1 = (/rho(i,j-2), u(i,j-2), v(i,j-2), p(i,j-2)/) 
      Q2 = (/rho(i,j-1), u(i,j-1), v(i,j-1), p(i,j-1)/) 
      Q5 = (/rho(i,j+2), u(i,j+2), v(i,j+2), p(i,j+2)/) 
      Q6 = (/rho(i,j+3), u(i,j+3), v(i,j+3), p(i,j+3)/) 
      call Qlr_mid(1.d0,1.d0,Q1,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
    elseif (j == 2) then
      ! left 3rd-order MUSCL
      ! right 4th-order MUSCL
      Q2 = (/rho(i,j-1), u(i,j-1), v(i,j-1), p(i,j-1)/) 
      Q5 = (/rho(i,j+2), u(i,j+2), v(i,j+2), p(i,j+2)/) 
      Q6 = (/rho(i,j+3), u(i,j+3), v(i,j+3), p(i,j+3)/) 
      call Qlr_mid(0.d0,1.d0,zero,Q2,Q3,Q4,Q5,Q6,Ql,Qr)
    elseif (j == ny-2) then
      ! left 4th-order MUSCL
      ! right 3rd-order MUSCL
      Q1 = (/rho(i,j-2), u(i,j-2), v(i,j-2), p(i,j-2)/) 
      Q2 = (/rho(i,j-1), u(i,j-1), v(i,j-1), p(i,j-1)/) 
      Q5 = (/rho(i,j+2), u(i,j+2), v(i,j+2), p(i,j+2)/) 
      call Qlr_mid(1.d0,0.d0,Q1,Q2,Q3,Q4,Q5,zero,Ql,Qr)
    elseif (j == 1) then
      Q5 = (/rho(i,j+2), u(i,j+2), v(i,j+2), p(i,j+2)/) 
      call Qlr_left(Q3,Q4,Q5,Ql,Qr)
    else
      Q2 = (/rho(i,j-1), u(i,j-1), v(i,j-1), p(i,j-1)/) 
      call Qlr_right(Q2,Q3,Q4,Ql,Qr)
    endif

    Jacobian = 0.5d0 * (Jacobians(i,j) + Jacobians(i,j+1))
    Normal(:) = (/0.d0, 0.d0, 0.5d0 * (etay(j) + etay(j+1)), 0.d0/)
    if (id_scheme == 2) then
      F(i-offset,j,:) = Roe(2,Ql,Qr,Normal,Jacobian)
    else
      F(i-offset,j,:) = SLAU(2,Ql,Qr,Normal,Jacobian)
    endif
  end subroutine calc_F_MUSCL_4th
end module calc_flux

