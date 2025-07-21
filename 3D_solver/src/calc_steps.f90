module calc_steps
  use cudafor
  use mod_globals, only : accuracy, offset, dt
  implicit none
  interface
    subroutine calc_step(nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, Q, Q2, Rs)
      integer, intent(in), value                                                                 :: nx, ny, nz
      real(8), intent(in), value                                                                 :: coef1, coef2
      real(8), intent(in), dimension(nx-1), device                                               :: dx
      real(8), intent(in), dimension(ny-1), device                                               :: dy
      real(8), intent(in), dimension(nz-1), device                                               :: dz
      real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device            :: E
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device            :: F
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device            :: G
      real(8), intent(in), dimension(nx,ny,nz,5), device                                         :: Q
      real(8), intent(out), dimension(nx,ny,nz,5), device                                        :: Q2
      real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), optional, device :: Rs
    end subroutine
  end interface
  
  interface
    subroutine calc_step_forcing(nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, fx, fy, fz, Q, Q2, Rs)
      integer, intent(in), value                                                                 :: nx, ny, nz
      real(8), intent(in), value                                                                 :: coef1, coef2
      real(8), intent(in), dimension(nx-1), device                                               :: dx
      real(8), intent(in), dimension(ny-1), device                                               :: dy
      real(8), intent(in), dimension(nz-1), device                                               :: dz
      real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device            :: E
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device            :: F
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device            :: G
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device              :: fx
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device              :: fy
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device              :: fz
      real(8), intent(in), dimension(nx,ny,nz,5), device                                         :: Q
      real(8), intent(out), dimension(nx,ny,nz,5), device                                        :: Q2
      real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), optional, device :: Rs
    end subroutine
  end interface
  
  interface
    subroutine calc_step2(nx, ny, nz, coef1, coef2, coef3, coef4, dx, dy, dz, E, F, G, Q, Q2, Rs)
      integer, intent(in), value                                                                 :: nx, ny, nz
      real(8), intent(in), value                                                                 :: coef1, coef2, coef3, coef4
      real(8), intent(in), dimension(nx-1), device                                               :: dx
      real(8), intent(in), dimension(ny-1), device                                               :: dy
      real(8), intent(in), dimension(nz-1), device                                               :: dz
      real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device            :: E
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device            :: F
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device            :: G
      real(8), intent(in), dimension(nx,ny,nz,5), device                                         :: Q
      real(8), intent(inout), dimension(nx,ny,nz,5), device                                      :: Q2
      real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), optional, device :: Rs
    end subroutine
  end interface
  
  interface
    subroutine calc_step2_forcing(nx, ny, nz, coef1, coef2, coef3, coef4, dx, dy, dz, E, F, G, fx, fy, fz, Q, Q2, Rs)
      integer, intent(in), value                                                                 :: nx, ny, nz
      real(8), intent(in), value                                                                 :: coef1, coef2, coef3, coef4
      real(8), intent(in), dimension(nx-1), device                                               :: dx
      real(8), intent(in), dimension(ny-1), device                                               :: dy
      real(8), intent(in), dimension(nz-1), device                                               :: dz
      real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device            :: E
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device            :: F
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device            :: G
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device              :: fx
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device              :: fy
      real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device              :: fz
      real(8), intent(in), dimension(nx,ny,nz,5), device                                         :: Q
      real(8), intent(inout), dimension(nx,ny,nz,5), device                                      :: Q2
      real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), optional, device :: Rs
    end subroutine
  end interface
contains
  subroutine calc_R(nx, ny, nz, dx, dy, dz, E, F, G, R)
    integer, intent(in), value                                                      :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device                                    :: dx
    real(8), intent(in), dimension(ny-1), device                                    :: dy
    real(8), intent(in), dimension(nz-1), device                                    :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5),device   :: R
    integer i, j, k, l
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1, nz-2*offset
        do j = 1, ny-2*offset
          do i = 1, nx-2*offset
            R(i,j,k,l) = &
            &   dt / (dy(j) * dz(k)) * (-E(i,j,k,l) + E(i+1,j,k,l)) &
            & + dt / (dz(k) * dx(i)) * (-F(i,j,k,l) + F(i,j+1,k,l)) &
            & + dt / (dx(i) * dy(j)) * (-G(i,j,k,l) + G(i,j,k+1,l))
    enddo;enddo;enddo;enddo
  end subroutine calc_R

  subroutine calc_R_forcing(nx, ny, nz, dx, dy, dz, E, F, G, fx, fy, fz, R)
    integer, intent(in), value                                                      :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device                                    :: dx
    real(8), intent(in), dimension(ny-1), device                                    :: dy
    real(8), intent(in), dimension(nz-1), device                                    :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device   :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device   :: fy
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device   :: fz
    real(8), intent(out), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5),device   :: R
    integer i, j, k, l
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1, nz-2*offset
        do j = 1, ny-2*offset
          do i = 1, nx-2*offset
            R(i,j,k,l) = &
            &   dt / (dy(j) * dz(k)) * (-E(i,j,k,l) + E(i+1,j,k,l) + dx(i) * fx(i,j,k,l)) &
            & + dt / (dz(k) * dx(i)) * (-F(i,j,k,l) + F(i,j+1,k,l) + dy(j) * fy(i,j,k,l)) &
            & + dt / (dx(i) * dy(j)) * (-G(i,j,k,l) + G(i,j,k+1,l) + dz(k) * fz(i,j,k,l))
    enddo;enddo;enddo;enddo
  end subroutine calc_R_forcing

  subroutine calc_step(nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, Q, Q2, Rs)
    integer, intent(in), value                                                                 :: nx, ny, nz
    real(8), intent(in), value                                                                 :: coef1, coef2
    real(8), intent(in), dimension(nx-1), device                                               :: dx
    real(8), intent(in), dimension(ny-1), device                                               :: dy
    real(8), intent(in), dimension(nz-1), device                                               :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device            :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device            :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device            :: G
    real(8), intent(in), dimension(nx,ny,nz,5), device                                         :: Q
    real(8), intent(out), dimension(nx,ny,nz,5), device                                        :: Q2
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), optional, device :: Rs
    integer i, j, k, l
    real(8) R
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1, nz-2*offset
        do j = 1, ny-2*offset
          do i = 1, nx-2*offset
            R = dt / (dy(j) * dz(k)) * (-E(i,j,k,l) + E(i+1,j,k,l)) &
            & + dt / (dz(k) * dx(i)) * (-F(i,j,k,l) + F(i,j+1,k,l)) &
            & + dt / (dx(i) * dy(j)) * (-G(i,j,k,l) + G(i,j,k+1,l))
            Q2(i+offset,j+offset,k+offset,l) = Q(i+offset,j+offset,k+offset,l) - coef1 * R
            if (present(Rs)) then
              Rs(i,j,k,l) = Rs(i,j,k,l) + coef2 * R
            endif
    enddo;enddo;enddo;enddo
  end subroutine calc_step
  
  subroutine calc_step_forcing(nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, fx, fy, fz, Q, Q2, Rs)
    integer, intent(in), value                                                                 :: nx, ny, nz
    real(8), intent(in), value                                                                 :: coef1, coef2
    real(8), intent(in), dimension(nx-1), device                                               :: dx
    real(8), intent(in), dimension(ny-1), device                                               :: dy
    real(8), intent(in), dimension(nz-1), device                                               :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device            :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device            :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device            :: G
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device              :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device              :: fy
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device              :: fz
    real(8), intent(in), dimension(nx,ny,nz,5), device                                         :: Q
    real(8), intent(out), dimension(nx,ny,nz,5), device                                        :: Q2
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), optional, device :: Rs
    integer i, j, k, l
    real(8) R
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1, nz-2*offset
        do j = 1, ny-2*offset
          do i = 1, nx-2*offset
            R = dt / (dy(j) * dz(k)) * (-E(i,j,k,l) + E(i+1,j,k,l) + dx(i) * fx(i,j,k,l)) &
            & + dt / (dz(k) * dx(i)) * (-F(i,j,k,l) + F(i,j+1,k,l) + dy(j) * fy(i,j,k,l)) &
            & + dt / (dx(i) * dy(j)) * (-G(i,j,k,l) + G(i,j,k+1,l) + dz(k) * fz(i,j,k,l))
            Q2(i+offset,j+offset,k+offset,l) = Q(i+offset,j+offset,k+offset,l) - coef1 * R
            if (present(Rs)) then
              Rs(i,j,k,l) = Rs(i,j,k,l) + coef2 * R
            endif
    enddo;enddo;enddo;enddo
  end subroutine calc_step_forcing
  
  subroutine calc_step2(nx, ny, nz, coef1, coef2, coef3, coef4, dx, dy, dz, E, F, G, Q, Q2, Rs)
    integer, intent(in), value                                                                 :: nx, ny, nz
    real(8), intent(in), value                                                                 :: coef1, coef2, coef3, coef4
    real(8), intent(in), dimension(nx-1), device                                               :: dx
    real(8), intent(in), dimension(ny-1), device                                               :: dy
    real(8), intent(in), dimension(nz-1), device                                               :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device            :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device            :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device            :: G
    real(8), intent(in), dimension(nx,ny,nz,5), device                                         :: Q
    real(8), intent(inout), dimension(nx,ny,nz,5), device                                      :: Q2
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), optional, device :: Rs
    integer i, j, k, l
    real(8) R
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1, nz-2*offset
        do j = 1, ny-2*offset
          do i = 1, nx-2*offset
            R = dt / (dy(j) * dz(k)) * (-E(i,j,k,l) + E(i+1,j,k,l)) &
            & + dt / (dz(k) * dx(i)) * (-F(i,j,k,l) + F(i,j+1,k,l)) &
            & + dt / (dx(i) * dy(j)) * (-G(i,j,k,l) + G(i,j,k+1,l))
            Q2(i+offset,j+offset,k+offset,l) = (coef1 * Q(i+offset,j+offset,k+offset,l) + coef2 * Q2(i+offset,j+offset,k+offset,l) - coef3 * R) / coef4
            if (present(Rs)) then
              Rs(i,j,k,l) = R
            endif
    enddo;enddo;enddo;enddo
  end subroutine calc_step2
  
  subroutine calc_step2_forcing(nx, ny, nz, coef1, coef2, coef3, coef4, dx, dy, dz, E, F, G, fx, fy, fz, Q, Q2, Rs)
    integer, intent(in), value                                                                 :: nx, ny, nz
    real(8), intent(in), value                                                                 :: coef1, coef2, coef3, coef4
    real(8), intent(in), dimension(nx-1), device                                               :: dx
    real(8), intent(in), dimension(ny-1), device                                               :: dy
    real(8), intent(in), dimension(nz-1), device                                               :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device            :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device            :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device            :: G
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device              :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device              :: fy
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device              :: fz
    real(8), intent(in), dimension(nx,ny,nz,5), device                                         :: Q
    real(8), intent(inout), dimension(nx,ny,nz,5), device                                      :: Q2
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), optional, device :: Rs
    integer i, j, k, l
    real(8) R
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1, nz-2*offset
        do j = 1, ny-2*offset
          do i = 1, nx-2*offset
            R = dt / (dy(j) * dz(k)) * (-E(i,j,k,l) + E(i+1,j,k,l) + dx(i) * fx(i,j,k,l)) &
            & + dt / (dz(k) * dx(i)) * (-F(i,j,k,l) + F(i,j+1,k,l) + dy(j) * fy(i,j,k,l)) &
            & + dt / (dx(i) * dy(j)) * (-G(i,j,k,l) + G(i,j,k+1,l) + dz(j) * fz(i,j,k,l))
            Q2(i+offset,j+offset,k+offset,l) = (coef1 * Q(i+offset,j+offset,k+offset,l) + coef2 * Q2(i+offset,j+offset,k+offset,l) - coef3 * R) / coef4
            if (present(Rs)) then
              Rs(i,j,k,l) = R
            endif
    enddo;enddo;enddo;enddo
  end subroutine calc_step2_forcing
  
  subroutine calc_step3(nx, ny, nz, dx, dy, dz, E, F, G, Q3, Q)
    integer, intent(in), value                                                      :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device                                    :: dx
    real(8), intent(in), dimension(ny-1), device                                    :: dy
    real(8), intent(in), dimension(nz-1), device                                    :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    real(8), intent(in), dimension(nx,ny,nz,5), device                              :: Q3
    real(8), intent(inout), dimension(nx,ny,nz,5), device                           :: Q
    integer i, j, k, l
    real(8) R
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1, nz-2*offset
        do j = 1, ny-2*offset
          do i = 1, nx-2*offset
            R = dt / (dy(j) * dz(k)) * (-E(i,j,k,l) + E(i+1,j,k,l)) &
            & + dt / (dz(k) * dx(i)) * (-F(i,j,k,l) + F(i,j+1,k,l)) &
            & + dt / (dx(i) * dy(j)) * (-G(i,j,k,l) + G(i,j,k+1,l))
            Q(i+offset,j+offset,k+offset,l) = (Q(i+offset,j+offset,k+offset,l) + 2.0d0 * Q3(i+offset,j+offset,k+offset,l) - 2.d0 * R) / 3.d0
    enddo;enddo;enddo;enddo
  end subroutine calc_step3
 
  subroutine calc_step3_forcing(nx, ny, nz, dx, dy, dz, E, F, G, fx, fy, fz, Q3, Q)
    integer, intent(in), value                                                      :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device                                    :: dx
    real(8), intent(in), dimension(ny-1), device                                    :: dy
    real(8), intent(in), dimension(nz-1), device                                    :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device :: G
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device   :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device   :: fy
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device   :: fz
    real(8), intent(in), dimension(nx,ny,nz,5), device                              :: Q3
    real(8), intent(inout), dimension(nx,ny,nz,5), device                           :: Q
    integer i, j, k, l
    real(8) R
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1, nz-2*offset
        do j = 1, ny-2*offset
          do i = 1, nx-2*offset
            R = dt / (dy(j) * dz(k)) * (-E(i,j,k,l) + E(i+1,j,k,l) + dx(i) * fx(i,j,k,l)) &
            & + dt / (dz(k) * dx(i)) * (-F(i,j,k,l) + F(i,j+1,k,l) + dy(j) * fy(i,j,k,l)) &
            & + dt / (dx(i) * dy(j)) * (-G(i,j,k,l) + G(i,j,k+1,l) + dz(k) * fz(i,j,k,l))
            Q(i+offset,j+offset,k+offset,l) = (Q(i+offset,j+offset,k+offset,l) + 2.0d0 * Q3(i+offset,j+offset,k+offset,l) - 2.d0 * R) / 3.d0
    enddo;enddo;enddo;enddo
  end subroutine calc_step3_forcing
 
  subroutine calc_step4(nx, ny, nz, dx, dy, dz, E, F, G, Rs, Q)
    integer, intent(in), value                                                       :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device                                     :: dx
    real(8), intent(in), dimension(ny-1), device                                     :: dy
    real(8), intent(in), dimension(nz-1), device                                     :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device  :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device  :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device  :: G
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device :: Rs
    real(8), intent(inout), dimension(nx,ny,nz,5), device                            :: Q
    integer i, j, k, l
    real(8) R
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1, nz-2*offset
        do j = 1, ny-2*offset
          do i = 1, nx-2*offset
            R = dt / (dy(j) * dz(k)) * (-E(i,j,k,l) + E(i+1,j,k,l)) &
            & + dt / (dz(k) * dx(i)) * (-F(i,j,k,l) + F(i,j+1,k,l)) &
            & + dt / (dx(i) * dy(j)) * (-G(i,j,k,l) + G(i,j,k+1,l))
            Rs(i,j,k,l) = Rs(i,j,k,l) + R
            Q(i+offset,j+offset,k+offset,l) = Q(i+offset,j+offset,k+offset,l) - Rs(i,j,k,l) / 6.d0
            Rs(i,j,k,l) = 0.d0
    enddo;enddo;enddo;enddo
  end subroutine calc_step4

  subroutine calc_step4_forcing(nx, ny, nz, dx, dy, dz, E, F, G, fx, fy, fz, Rs, Q)
    integer, intent(in), value                                                       :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device                                     :: dx
    real(8), intent(in), dimension(ny-1), device                                     :: dy
    real(8), intent(in), dimension(nz-1), device                                     :: dz
    real(8), intent(in), dimension(nx-accuracy+1,ny-accuracy,nz-accuracy,5), device  :: E
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy+1,nz-accuracy,5), device  :: F
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy+1,5), device  :: G
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device    :: fx
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device    :: fy
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device    :: fz
    real(8), intent(inout), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device :: Rs
    real(8), intent(inout), dimension(nx,ny,nz,5), device                            :: Q
    integer i, j, k, l
    real(8) R
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1, nz-2*offset
        do j = 1, ny-2*offset
          do i = 1, nx-2*offset
            R = dt / (dy(j) * dz(k)) * (-E(i,j,k,l) + E(i+1,j,k,l) + dx(i) * fx(i,j,k,l)) &
            & + dt / (dz(k) * dx(i)) * (-F(i,j,k,l) + F(i,j+1,k,l) + dy(j) * fy(i,j,k,l)) &
            & + dt / (dx(i) * dy(j)) * (-G(i,j,k,l) + G(i,j,k+1,l) + dz(k) * fz(i,j,k,l))
            Rs(i,j,k,l) = Rs(i,j,k,l) + R
            Q(i+offset,j+offset,k+offset,l) = Q(i+offset,j+offset,k+offset,l) - Rs(i,j,k,l) / 6.d0
            Rs(i,j,k,l) = 0.d0
    enddo;enddo;enddo;enddo
  end subroutine calc_step4_forcing

  subroutine calc_error(nx, ny, nz, R1, R2, R1_new, R2_new, err)
    integer, intent(in), value                                                    :: nx, ny, nz
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device :: R1, R2, R1_new, R2_new
    real(8), intent(out)                                                          :: err
    integer i, j, k, l
    err = 0.d0
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1, nz-2*offset
        do j = 1, ny-2*offset
          do i = 1, nx-2*offset
            err = err + sqrt((R1(i,j,k,l) - R1_new(i,j,k,l)**2)) + sqrt((R2(i,j,k,l) - R2_new(i,j,k,l))**2)
    enddo;enddo;enddo;enddo
    err = err / (dble(nx - accuracy) * dble(ny * accuracy) * dble(nz * accuracy) * 5.d0)
  end subroutine calc_error

  subroutine calc_Gauss_step(nx, ny, nz, a1, a2, dx, dy, dz, R1, R2, Q, Q2)
    integer, intent(in), value                                                    :: nx, ny, nz
    real(8), intent(in), value                                                    :: a1, a2
    real(8), intent(in), dimension(nx-1), device                                  :: dx
    real(8), intent(in), dimension(ny-1), device                                  :: dy
    real(8), intent(in), dimension(nz-1), device                                  :: dz
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device :: R1
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device :: R2
    real(8), intent(in), dimension(nx,ny,nz,5), device                            :: Q
    real(8), intent(out), dimension(nx,ny,nz,5), device                           :: Q2
    integer i, j, k, l
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1+offset, nz-offset
        do j = 1+offset, ny-offset
          do i = 1+offset, nx-offset
            Q2(i,j,k,l) = Q(i,j,k,l) - (a1 * R1(i-offset,j-offset,k-offset,l) + a2 * R2(i-offset,j-offset,k-offset,l))
    enddo;enddo;enddo;enddo
  end subroutine calc_Gauss_step
  
  subroutine calc_Gauss_step_Q(nx, ny, nz, a1, a2, dx, dy, dz, R1, R2, Q)
    integer, intent(in), value                                                    :: nx, ny, nz
    real(8), intent(in), value                                                    :: a1, a2
    real(8), intent(in), dimension(nx-1), device                                  :: dx
    real(8), intent(in), dimension(ny-1), device                                  :: dy
    real(8), intent(in), dimension(nz-1), device                                  :: dz
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device :: R1
    real(8), intent(in), dimension(nx-accuracy,ny-accuracy,nz-accuracy,5), device :: R2
    real(8), intent(inout), dimension(nx,ny,nz,5), device                         :: Q
    integer i, j, k, l
    !$cuf kernel do(4) <<<*,*>>>
    do l = 1, 5
      do k = 1+offset, nz-offset
        do j = 1+offset, ny-offset
          do i = 1+offset, nx-offset
            Q(i,j,k,l) = Q(i,j,k,l) - (a1 * R1(i-offset,j-offset,k-offset,l) + a2 * R2(i-offset,j-offset,k-offset,l))
    enddo;enddo;enddo;enddo
  end subroutine calc_Gauss_step_Q
end module calc_steps

