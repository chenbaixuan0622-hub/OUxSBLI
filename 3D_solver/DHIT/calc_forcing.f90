module calc_forcing
  use cufft
  use cudafor
  use mod_globals, only : nx, ny, nz, eps_s, kf_min, kf_max, threads, &
                          gamma, R, C_T, dt, T_ref_const => T
  use mod_constant, only : id_accuracy
  implicit none
  integer,    save :: Nf, off
  integer,    save :: fft_plan
  real(8),    device, allocatable, save :: fx_d(:,:,:), fy_d(:,:,:), fz_d(:,:,:)
  complex(8), device, allocatable, save :: ux_k(:,:,:), uy_k(:,:,:), uz_k(:,:,:)
contains

  subroutine init_forcing()
    integer :: istat
    if     (kind(id_accuracy) == 2) then; off = 1
    elseif (kind(id_accuracy) == 4) then; off = 2
    elseif (kind(id_accuracy) == 8) then; off = 3
    end if
    Nf = nx - 2*off
    allocate(fx_d(nx,ny,nz), fy_d(nx,ny,nz), fz_d(nx,ny,nz))
    allocate(ux_k(Nf,Nf,Nf), uy_k(Nf,Nf,Nf), uz_k(Nf,Nf,Nf))
    fx_d = 0.d0;  fy_d = 0.d0;  fz_d = 0.d0
    istat = cufftPlan3d(fft_plan, Nf, Nf, Nf, CUFFT_Z2Z)
    print *, '[forcing] init: Nf =', Nf, '  eps_s =', eps_s, &
             '  kf_min =', kf_min, '  kf_max =', kf_max
  end subroutine init_forcing


  ! GPU kernel: extract velocity components u,v,w from QJ into complex device arrays
  attributes(global) subroutine extract_vel_k(nx_, ny_, nz_, off_, QJ, ux, uy, uz)
    integer, intent(in), value                 :: nx_, ny_, nz_, off_
    real(8),    intent(in),  device, contiguous :: QJ(nx_,5,ny_,nz_)
    complex(8), intent(out), device, contiguous :: ux(nx_-2*off_, ny_-2*off_, nz_-2*off_)
    complex(8), intent(out), device, contiguous :: uy(nx_-2*off_, ny_-2*off_, nz_-2*off_)
    complex(8), intent(out), device, contiguous :: uz(nx_-2*off_, ny_-2*off_, nz_-2*off_)
    integer :: i, j, k
    real(8) :: rho_inv
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx_-2*off_ < i .or. ny_-2*off_ < j .or. nz_-2*off_ < k) return
    rho_inv = 1.d0 / QJ(i+off_, 1, j+off_, k+off_)
    ux(i,j,k) = cmplx(QJ(i+off_, 2, j+off_, k+off_) * rho_inv, 0.d0, 8)
    uy(i,j,k) = cmplx(QJ(i+off_, 3, j+off_, k+off_) * rho_inv, 0.d0, 8)
    uz(i,j,k) = cmplx(QJ(i+off_, 4, j+off_, k+off_) * rho_inv, 0.d0, 8)
  end subroutine extract_vel_k


  ! GPU kernel: normalize IFFT output and store real part into full-grid forcing arrays
  attributes(global) subroutine store_forcing_k(nx_, ny_, nz_, off_, norm_, ux, uy, uz, fx, fy, fz)
    integer, intent(in), value                 :: nx_, ny_, nz_, off_
    real(8), intent(in), value                 :: norm_
    complex(8), intent(in),  device, contiguous :: ux(nx_-2*off_, ny_-2*off_, nz_-2*off_)
    complex(8), intent(in),  device, contiguous :: uy(nx_-2*off_, ny_-2*off_, nz_-2*off_)
    complex(8), intent(in),  device, contiguous :: uz(nx_-2*off_, ny_-2*off_, nz_-2*off_)
    real(8),    intent(out), device, contiguous :: fx(nx_,ny_,nz_)
    real(8),    intent(out), device, contiguous :: fy(nx_,ny_,nz_)
    real(8),    intent(out), device, contiguous :: fz(nx_,ny_,nz_)
    integer :: i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx_-2*off_ < i .or. ny_-2*off_ < j .or. nz_-2*off_ < k) return
    fx(i+off_, j+off_, k+off_) = real(ux(i,j,k)) * norm_
    fy(i+off_, j+off_, k+off_) = real(uy(i,j,k)) * norm_
    fz(i+off_, j+off_, k+off_) = real(uz(i,j,k)) * norm_
  end subroutine store_forcing_k


  ! GPU kernel: relax energy toward T_ref — energy equation only, no momentum change
  attributes(global) subroutine add_cooling_k(nx_, ny_, nz_, coef_, QJ)
    integer, intent(in), value                 :: nx_, ny_, nz_
    real(8), intent(in), value                 :: coef_
    real(8), intent(inout), device, contiguous :: QJ(nx_,5,ny_,nz_)
    integer :: i, j, k
    real(8) :: rho, u, v, w, e_int, T_loc
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx_-2 < i .or. ny_-2 < j .or. nz_-2 < k) return
    rho   = QJ(i+1, 1, j+1, k+1)
    u     = QJ(i+1, 2, j+1, k+1) / rho
    v     = QJ(i+1, 3, j+1, k+1) / rho
    w     = QJ(i+1, 4, j+1, k+1) / rho
    e_int = QJ(i+1, 5, j+1, k+1) / rho - 0.5d0*(u*u + v*v + w*w)
    T_loc = e_int * (gamma - 1.d0) / R
    QJ(i+1, 5, j+1, k+1) = QJ(i+1, 5, j+1, k+1) &
                           - coef_ * dt * rho * C_T * (T_loc - T_ref_const)
  end subroutine add_cooling_k


  subroutine calc_forcing_rhs(nx_, ny_, nz_, QJ)
    integer, intent(in)  :: nx_, ny_, nz_
    real(8), device      :: QJ(nx_,5,ny_,nz_)
    integer :: istat, kxi, kyi, kzi
    integer(8) :: ikx, iky, ikz
    real(8) :: kmag2, kmag, kdu_re, kdu_im, Ks, As, norm
    complex(8) :: kdotu
    complex(8), allocatable :: ux_h(:,:,:), uy_h(:,:,:), uz_h(:,:,:)
    type(dim3) :: blk

    blk = dim3((Nf+31)/32, (Nf+3)/4, Nf)

    ! Step 1: extract velocity from QJ (device→device via GPU kernel)
    call extract_vel_k<<<blk, threads>>>(nx_, ny_, nz_, off, QJ, ux_k, uy_k, uz_k)

    ! Step 2: forward 3D FFT (in-place)
    istat = cufftExecZ2Z(fft_plan, ux_k, ux_k, CUFFT_FORWARD)
    istat = cufftExecZ2Z(fft_plan, uy_k, uy_k, CUFFT_FORWARD)
    istat = cufftExecZ2Z(fft_plan, uz_k, uz_k, CUFFT_FORWARD)

    ! Step 3: copy spectral data to host
    allocate(ux_h(Nf,Nf,Nf), uy_h(Nf,Nf,Nf), uz_h(Nf,Nf,Nf))
    ux_h = ux_k;  uy_h = uy_k;  uz_h = uz_k

    ! Step 4: solenoidal projection, forcing band mask, K_s reduction (CPU)
    Ks = 0.d0
    do kzi = 1, Nf
      ikz = kzi - 1;  if (ikz > Nf/2) ikz = ikz - Nf
      do kyi = 1, Nf
        iky = kyi - 1;  if (iky > Nf/2) iky = iky - Nf
        do kxi = 1, Nf
          ikx = kxi - 1;  if (ikx > Nf/2) ikx = ikx - Nf
          kmag2 = dble(ikx*ikx + iky*iky + ikz*ikz)
          if (kmag2 < 0.5d0) cycle
          kmag = sqrt(kmag2)
          ! Solenoidal projection: û_s = û - k(k·û)/|k|²
          kdu_re = dble(ikx)*real(ux_h(kxi,kyi,kzi)) + dble(iky)*real(uy_h(kxi,kyi,kzi)) &
                 + dble(ikz)*real(uz_h(kxi,kyi,kzi))
          kdu_im = dble(ikx)*aimag(ux_h(kxi,kyi,kzi)) + dble(iky)*aimag(uy_h(kxi,kyi,kzi)) &
                 + dble(ikz)*aimag(uz_h(kxi,kyi,kzi))
          kdotu = cmplx(kdu_re, kdu_im, 8) / kmag2
          ux_h(kxi,kyi,kzi) = ux_h(kxi,kyi,kzi) - dble(ikx)*kdotu
          uy_h(kxi,kyi,kzi) = uy_h(kxi,kyi,kzi) - dble(iky)*kdotu
          uz_h(kxi,kyi,kzi) = uz_h(kxi,kyi,kzi) - dble(ikz)*kdotu
          ! Forcing band mask and K_s accumulation
          if (kmag < dble(kf_min)-0.5d0 .or. kmag > dble(kf_max)+0.5d0) then
            ux_h(kxi,kyi,kzi) = 0.d0;  uy_h(kxi,kyi,kzi) = 0.d0;  uz_h(kxi,kyi,kzi) = 0.d0
          else
            Ks = Ks + abs(ux_h(kxi,kyi,kzi))**2 + abs(uy_h(kxi,kyi,kzi))**2 &
                    + abs(uz_h(kxi,kyi,kzi))**2
          end if
        enddo
      enddo
    enddo
    ! Parseval's theorem: K_s = 0.5 * (1/Nf^3)^2 * sum|û_s|^2
    Ks = 0.5d0 * Ks / dble(Nf)**6

    ! Step 5: compute A_s = eps_s / (2*K_s) and scale spectral modes
    As = eps_s / max(2.d0 * Ks, 1.d-30)
    ux_h = As * ux_h;  uy_h = As * uy_h;  uz_h = As * uz_h

    ! Step 6: copy scaled spectral forcing back to device
    ux_k = ux_h;  uy_k = uy_h;  uz_k = uz_h
    deallocate(ux_h, uy_h, uz_h)

    ! Step 7: inverse FFT (in-place)
    istat = cufftExecZ2Z(fft_plan, ux_k, ux_k, CUFFT_INVERSE)
    istat = cufftExecZ2Z(fft_plan, uy_k, uy_k, CUFFT_INVERSE)
    istat = cufftExecZ2Z(fft_plan, uz_k, uz_k, CUFFT_INVERSE)

    ! Step 8: normalize and store real part in full-grid forcing arrays
    norm = 1.d0 / dble(Nf)**3
    call store_forcing_k<<<blk, threads>>>(nx_, ny_, nz_, off, norm, ux_k, uy_k, uz_k, &
                                            fx_d, fy_d, fz_d)
  end subroutine calc_forcing_rhs


  subroutine apply_cooling(nx_, ny_, nz_, coef, QJ)
    integer, intent(in) :: nx_, ny_, nz_
    real(8), intent(in) :: coef
    real(8), device     :: QJ(nx_,5,ny_,nz_)
    type(dim3) :: blk
    blk = dim3((nx_-2+31)/32, (ny_-2+3)/4, nz_-2)
    call add_cooling_k<<<blk, threads>>>(nx_, ny_, nz_, coef, QJ)
  end subroutine apply_cooling


  subroutine finalize_forcing()
    integer :: istat
    istat = cufftDestroy(fft_plan)
    deallocate(fx_d, fy_d, fz_d, ux_k, uy_k, uz_k)
  end subroutine finalize_forcing

end module calc_forcing
