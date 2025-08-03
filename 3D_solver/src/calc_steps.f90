module calc_steps
  use cudafor
  use mod_globals, only : accuracy, offset, dt
  implicit none
contains
  subroutine calc_R(nx, ny, nz, dx, dy, dz, E, F, G, R)
    integer, intent(in), value                                                      :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device                                    :: dx
    real(8), intent(in), dimension(ny-1), device                                    :: dy
    real(8), intent(in), dimension(nz-1), device                                    :: dz
    real(8), intent(in), dimension(5,nx-accuracy+1,ny-accuracy,nz-accuracy), device :: E
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy+1,nz-accuracy), device :: F
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy+1), device :: G
    real(8), intent(out), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy),device   :: R
    real(8) dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2*offset
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          dtdydz = dt / (dy(j) * dz(k))
          dtdzdx = dt / (dz(k) * dx(i))
          dtdxdy = dt / (dx(i) * dy(j))
          do l = 1, 5
            R(l,i,j,k) = &
            &   dtdydz * (-E(l,i,j,k) + E(l,i+1,j,k)) &
            & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
            & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
    enddo;enddo;enddo;enddo
  end subroutine calc_R

  subroutine calc_R_forcing(nx, ny, nz, dx, dy, dz, E, F, G, fx, fy, fz, R)
    integer, intent(in), value                                                      :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device                                    :: dx
    real(8), intent(in), dimension(ny-1), device                                    :: dy
    real(8), intent(in), dimension(nz-1), device                                    :: dz
    real(8), intent(in), dimension(5,nx-accuracy+1,ny-accuracy,nz-accuracy), device :: E
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy+1,nz-accuracy), device :: F
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy+1), device :: G
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy), device     :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy), device     :: fy
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy), device     :: fz
    real(8), intent(out), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy),device   :: R
    real(8) dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2*offset
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          dtdydz = dt / (dy(j) * dz(k))
          dtdzdx = dt / (dz(k) * dx(i))
          dtdxdy = dt / (dx(i) * dy(j))
          do l = 1, 5
            R(l,i,j,k) = &
            &   dtdydz * (-E(l,i,j,k) + E(l,i+1,j,k) + fx(i,j,k) / dx(i)) &
            & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k) + fy(i,j,k) / dy(j)) &
            & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1) + fz(i,j,k) / dz(k))
    enddo;enddo;enddo;enddo
  end subroutine calc_R_forcing

  subroutine calc_step1(nx, ny, nz, coef, dx, dy, dz, E, F, G, Q, Q2)
    integer, intent(in), value                                                      :: nx, ny, nz
    real(8), intent(in), value                                                      :: coef
    real(8), intent(in), dimension(nx-1), device                                    :: dx
    real(8), intent(in), dimension(ny-1), device                                    :: dy
    real(8), intent(in), dimension(nz-1), device                                    :: dz
    real(8), intent(in), dimension(5,nx-accuracy+1,ny-accuracy,nz-accuracy), device :: E
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy+1,nz-accuracy), device :: F
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy+1), device :: G
    real(8), intent(in), dimension(5,nx,ny,nz), device                              :: Q
    real(8), intent(out), dimension(5,nx,ny,nz), device                             :: Q2
    real(8) R, dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2*offset
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          dtdydz = dt / (dy(j) * dz(k))
          dtdzdx = dt / (dz(k) * dx(i))
          dtdxdy = dt / (dx(i) * dy(j))
          do l = 1, 5
            R = dtdydz * (-E(l,i,j,k) + E(l,i+1,j,k)) &
            & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
            & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
            Q2(l,i+offset,j+offset,k+offset) = Q(l,i+offset,j+offset,k+offset) - coef * R
    enddo;enddo;enddo;enddo
  end subroutine calc_step1
  
  subroutine calc_step1_forcing(nx, ny, nz, coef, dx, dy, dz, E, F, G, fx, fy, fz, Q, Q2)
    integer, intent(in), value                                                      :: nx, ny, nz
    real(8), intent(in), value                                                      :: coef
    real(8), intent(in), dimension(nx-1), device                                    :: dx
    real(8), intent(in), dimension(ny-1), device                                    :: dy
    real(8), intent(in), dimension(nz-1), device                                    :: dz
    real(8), intent(in), dimension(5,nx-accuracy+1,ny-accuracy,nz-accuracy), device :: E
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy+1,nz-accuracy), device :: F
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy+1), device :: G
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy), device     :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy), device     :: fy
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy), device     :: fz
    real(8), intent(in), dimension(5,nx,ny,nz), device                              :: Q
    real(8), intent(out), dimension(5,nx,ny,nz), device                             :: Q2
    real(8) R, dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2*offset
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          dtdydz = dt / (dy(j) * dz(k))
          dtdzdx = dt / (dz(k) * dx(i))
          dtdxdy = dt / (dx(i) * dy(j))
          do l = 1, 5
            R = dtdydz * (-E(l,i,j,k) + E(l,i+1,j,k) + fx(i,j,k) / dx(i)) &
            & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k) + fy(i,j,k) / dy(j)) &
            & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1) + fz(i,j,k) / dz(k))
            Q2(l,i+offset,j+offset,k+offset) = Q(l,i+offset,j+offset,k+offset) - coef * R
    enddo;enddo;enddo;enddo
  end subroutine calc_step1_forcing
  
  subroutine calc_step(nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, Q, Q2, Rs)
    integer, intent(in), value                                                       :: nx, ny, nz
    real(8), intent(in), value                                                       :: coef1, coef2
    real(8), intent(in), dimension(nx-1), device                                     :: dx
    real(8), intent(in), dimension(ny-1), device                                     :: dy
    real(8), intent(in), dimension(nz-1), device                                     :: dz
    real(8), intent(in), dimension(5,nx-accuracy+1,ny-accuracy,nz-accuracy), device  :: E
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy+1,nz-accuracy), device  :: F
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy+1), device  :: G
    real(8), intent(in), dimension(5,nx,ny,nz), device                               :: Q
    real(8), intent(out), dimension(5,nx,ny,nz), device                              :: Q2
    real(8), intent(inout), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy), device :: Rs
    real(8) R, dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2*offset
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          dtdydz = dt / (dy(j) * dz(k))
          dtdzdx = dt / (dz(k) * dx(i))
          dtdxdy = dt / (dx(i) * dy(j))
          do l = 1, 5
            R = dtdydz * (-E(l,i,j,k) + E(l,i+1,j,k)) &
            & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
            & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
            Q2(l,i+offset,j+offset,k+offset) = Q(l,i+offset,j+offset,k+offset) - coef1 * R
            Rs(l,i,j,k) = Rs(l,i,j,k) + coef2 * R
    enddo;enddo;enddo;enddo
  end subroutine calc_step
  
  subroutine calc_step_forcing(nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, fx, fy, fz, Q, Q2, Rs)
    integer, intent(in), value                                                       :: nx, ny, nz
    real(8), intent(in), value                                                       :: coef1, coef2
    real(8), intent(in), dimension(nx-1), device                                     :: dx
    real(8), intent(in), dimension(ny-1), device                                     :: dy
    real(8), intent(in), dimension(nz-1), device                                     :: dz
    real(8), intent(in), dimension(5,nx-accuracy+1,ny-accuracy,nz-accuracy), device  :: E
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy+1,nz-accuracy), device  :: F
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy+1), device  :: G
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy), device      :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy), device      :: fy
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy), device      :: fz
    real(8), intent(in), dimension(5,nx,ny,nz), device                               :: Q
    real(8), intent(out), dimension(5,nx,ny,nz), device                              :: Q2
    real(8), intent(inout), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy), device :: Rs
    real(8) R, dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2*offset
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          dtdydz = dt / (dy(j) * dz(k))
          dtdzdx = dt / (dz(k) * dx(i))
          dtdxdy = dt / (dx(i) * dy(j))
          do l = 1, 5
            R = dtdydz * (-E(l,i,j,k) + E(l,i+1,j,k) + fx(i,j,k) / dx(i)) &
            & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k) + fy(i,j,k) / dy(j)) &
            & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1) + fz(i,j,k) / dz(k))
            Q2(l,i+offset,j+offset,k+offset) = Q(l,i+offset,j+offset,k+offset) - coef1 * R
            Rs(l,i,j,k) = Rs(l,i,j,k) + coef2 * R
    enddo;enddo;enddo;enddo
  end subroutine calc_step_forcing
  
  subroutine calc_step2_3(nx, ny, nz, coef1, coef2, coef3, coef4, dx, dy, dz, E, F, G, Qin, Qout)
    integer, intent(in), value                                                      :: nx, ny, nz
    real(8), intent(in), value                                                      :: coef1, coef2, coef3, coef4
    real(8), intent(in), dimension(nx-1), device                                    :: dx
    real(8), intent(in), dimension(ny-1), device                                    :: dy
    real(8), intent(in), dimension(nz-1), device                                    :: dz
    real(8), intent(in), dimension(5,nx-accuracy+1,ny-accuracy,nz-accuracy), device :: E
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy+1,nz-accuracy), device :: F
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy+1), device :: G
    real(8), intent(in), dimension(5,nx,ny,nz), device                              :: Qin
    real(8), intent(inout), dimension(5,nx,ny,nz), device                           :: Qout
    real(8) R, dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2*offset
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          dtdydz = dt / (dy(j) * dz(k))
          dtdzdx = dt / (dz(k) * dx(i))
          dtdxdy = dt / (dx(i) * dy(j))
          do l = 1, 5
            R = dtdydz * (-E(l,i,j,k) + E(l,i+1,j,k)) &
            & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
            & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
            Qout(l,i+offset,j+offset,k+offset) = (coef1 * Qin(l,i+offset,j+offset,k+offset) + coef2 * Qout(l,i+offset,j+offset,k+offset) - coef3 * R) / coef4
    enddo;enddo;enddo;enddo
  end subroutine calc_step2_3
  
  subroutine calc_step2_3_forcing(nx, ny, nz, coef1, coef2, coef3, coef4, dx, dy, dz, E, F, G, fx, fy, fz, Qin, Qout)
    integer, intent(in), value                                                      :: nx, ny, nz
    real(8), intent(in), value                                                      :: coef1, coef2, coef3, coef4
    real(8), intent(in), dimension(nx-1), device                                    :: dx
    real(8), intent(in), dimension(ny-1), device                                    :: dy
    real(8), intent(in), dimension(nz-1), device                                    :: dz
    real(8), intent(in), dimension(5,nx-accuracy+1,ny-accuracy,nz-accuracy), device :: E
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy+1,nz-accuracy), device :: F
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy+1), device :: G
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy), device     :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy), device     :: fy
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy), device     :: fz
    real(8), intent(in), dimension(5,nx,ny,nz), device                              :: Qin
    real(8), intent(inout), dimension(5,nx,ny,nz), device                           :: Qout
    real(8) R, dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2*offset
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          dtdydz = dt / (dy(j) * dz(k))
          dtdzdx = dt / (dz(k) * dx(i))
          dtdxdy = dt / (dx(i) * dy(j))
          do l = 1, 5
            R = dtdydz * (-E(l,i,j,k) + E(l,i+1,j,k) + fx(i,j,k) / dx(i)) &
            & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k) + fy(i,j,k) / dy(j)) &
            & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1) + fz(i,j,k) / dz(k))
            Qout(l,i+offset,j+offset,k+offset) = (coef1 * Qin(l,i+offset,j+offset,k+offset) + coef2 * Qout(l,i+offset,j+offset,k+offset) - coef3 * R) / coef4
    enddo;enddo;enddo;enddo
  end subroutine calc_step2_3_forcing
  
  subroutine calc_step4(nx, ny, nz, dx, dy, dz, E, F, G, Rs, Q)
    integer, intent(in), value                                                       :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device                                     :: dx
    real(8), intent(in), dimension(ny-1), device                                     :: dy
    real(8), intent(in), dimension(nz-1), device                                     :: dz
    real(8), intent(in), dimension(5,nx-accuracy+1,ny-accuracy,nz-accuracy), device  :: E
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy+1,nz-accuracy), device  :: F
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy+1), device  :: G
    real(8), intent(inout), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy), device :: Rs
    real(8), intent(inout), dimension(5,nx,ny,nz), device                            :: Q
    real(8) R, dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2*offset
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          dtdydz = dt / (dy(j) * dz(k))
          dtdzdx = dt / (dz(k) * dx(i))
          dtdxdy = dt / (dx(i) * dy(j))
          do l = 1, 5
            R = dtdydz * (-E(l,i,j,k) + E(l,i+1,j,k)) &
            & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
            & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
            Rs(l,i,j,k) = Rs(l,i,j,k) + R
            Q(l,i+offset,j+offset,k+offset) = Q(l,i+offset,j+offset,k+offset) - Rs(l,i,j,k) / 6.d0
            Rs(l,i,j,k) = 0.d0
    enddo;enddo;enddo;enddo
  end subroutine calc_step4

  subroutine calc_step4_forcing(nx, ny, nz, dx, dy, dz, E, F, G, fx, fy, fz, Rs, Q)
    integer, intent(in), value                                                       :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device                                     :: dx
    real(8), intent(in), dimension(ny-1), device                                     :: dy
    real(8), intent(in), dimension(nz-1), device                                     :: dz
    real(8), intent(in), dimension(5,nx-accuracy+1,ny-accuracy,nz-accuracy), device  :: E
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy+1,nz-accuracy), device  :: F
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy+1), device  :: G
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy), device      :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy), device      :: fy
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy), device      :: fz
    real(8), intent(inout), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy), device :: Rs
    real(8), intent(inout), dimension(5,nx,ny,nz), device                            :: Q
    real(8) R, dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2*offset
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          dtdydz = dt / (dy(j) * dz(k))
          dtdzdx = dt / (dz(k) * dx(i))
          dtdxdy = dt / (dx(i) * dy(j))
          do l = 1, 5
            R = dtdydz * (-E(l,i,j,k) + E(l,i+1,j,k) + fx(i,j,k) / dx(i)) &
            & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k) + fy(i,j,k) / dy(j)) &
            & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1) + fz(i,j,k) / dz(k))
            Rs(l,i,j,k) = Rs(l,i,j,k) + R
            Q(l,i+offset,j+offset,k+offset) = Q(l,i+offset,j+offset,k+offset) - Rs(l,i,j,k) / 6.d0
            Rs(l,i,j,k) = 0.d0
    enddo;enddo;enddo;enddo
  end subroutine calc_step4_forcing

  subroutine calc_error(nx, ny, nz, R1, R2, R1_new, R2_new, err)
    integer, intent(in), value                                                    :: nx, ny, nz
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy), device :: R1, R2, R1_new, R2_new
    real(8), intent(out)                                                          :: err
    integer i, j, k, l
    err = 0.d0
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2*offset
      do j = 1, ny-2*offset
        do i = 1, nx-2*offset
          do l = 1, 5
            err = err + sqrt((R1(l,i,j,k) - R1_new(l,i,j,k)**2)) + sqrt((R2(l,i,j,k) - R2_new(l,i,j,k))**2)
    enddo;enddo;enddo;enddo
    err = err / (dble(nx - accuracy) * dble(ny * accuracy) * dble(nz * accuracy) * 5.d0)
  end subroutine calc_error

  subroutine calc_Gauss_step(nx, ny, nz, a1, a2, dx, dy, dz, R1, R2, Q, Q2)
    integer, intent(in), value                                                    :: nx, ny, nz
    real(8), intent(in), value                                                    :: a1, a2
    real(8), intent(in), dimension(nx-1), device                                  :: dx
    real(8), intent(in), dimension(ny-1), device                                  :: dy
    real(8), intent(in), dimension(nz-1), device                                  :: dz
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy), device :: R1
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy), device :: R2
    real(8), intent(in), dimension(5,nx,ny,nz), device                            :: Q
    real(8), intent(out), dimension(5,nx,ny,nz), device                           :: Q2
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          do l = 1, 5
            Q2(l,i,j,k) = Q(l,i,j,k) - (a1 * R1(l,i-offset,j-offset,k-offset) + a2 * R2(l,i-offset,j-offset,k-offset))
    enddo;enddo;enddo;enddo
  end subroutine calc_Gauss_step
  
  subroutine calc_Gauss_step_Q(nx, ny, nz, a1, a2, dx, dy, dz, R1, R2, Q)
    integer, intent(in), value                                                    :: nx, ny, nz
    real(8), intent(in), value                                                    :: a1, a2
    real(8), intent(in), dimension(nx-1), device                                  :: dx
    real(8), intent(in), dimension(ny-1), device                                  :: dy
    real(8), intent(in), dimension(nz-1), device                                  :: dz
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy), device :: R1
    real(8), intent(in), dimension(5,nx-accuracy,ny-accuracy,nz-accuracy), device :: R2
    real(8), intent(inout), dimension(5,nx,ny,nz), device                         :: Q
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          do l = 1, 5
            Q(l,i,j,k) = Q(l,i,j,k) - (a1 * R1(l,i-offset,j-offset,k-offset) + a2 * R2(l,i-offset,j-offset,k-offset))
    enddo;enddo;enddo;enddo
  end subroutine calc_Gauss_step_Q
end module calc_steps

