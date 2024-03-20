module calc_riemann_solver
  use mod_globals, only : accuracy, offset, nx, ny
  use calc_qlr
  use calc_slau
  use calc_roe
  implicit none
contains
  attributes(global) subroutine calc_E(id_accuracy, rho, u, v, p, E)
    integer(kind=accuracy), intent(in), value :: id_accuracy
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(out), dimension(nx-accuracy+1,ny-accuracy,4), device :: E
    integer i, j
    real(8) :: Normal(4) = (/0.d0, 1.d0, 0.d0, 0.d0/)
    real(8), dimension(4) :: Q1, Q2, Q3, Q4, Ql, Qr
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

    !E(i,j-offset,:) = SLAU(1,Ql,Qr,Normal)
    E(i,j-offset,:) = Roe(1,Ql,Qr,Normal)
  end subroutine calc_E

  attributes(global) subroutine calc_F(id_accuracy, rho, u, v, p, F)
    integer(kind=accuracy), intent(in), value :: id_accuracy
    real(8), intent(in), dimension(nx,ny), device :: rho, u, v, p
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy+1,4), device :: F
    integer i, j
    real(8) :: Normal(4) = (/0.d0, 0.d0, 1.d0, 0.d0/)
    real(8), dimension(4) :: Q1, Q2, Q3, Q4, Ql, Qr
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

    !F(i-offset,j,:) = SLAU(2,Ql,Qr,Normal)
    F(i-offset,j,:) = Roe(2,Ql,Qr,Normal)
  end subroutine calc_F
end module calc_riemann_solver

