module set_init_dhit
  use cufft
  use mod_globals, only : id_accuracy, gamma, RHO0, p0, urms, pi, Re_lambda, &
                          pope_L, pope_eta, pope_cL, pope_p0, pope_beta, pope_ceta
  implicit none
contains
  subroutine init_spectral_velocity(nx, ny, nz, Q)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(out) :: Q(nx,5,ny,nz)
    integer   :: offset, Nf, plan, istat, icomp, seed_size
    integer, allocatable :: seed_arr(:)
    integer   :: kxi, kyi, kzi, kxc, kyc, kzc, i, j, k
    integer(8) :: ikx, iky, ikz
    real(8)   :: kmag2, kmag, amp, r(2), ph, kdu_re, kdu_im
    real(8)   :: kL, keta, fL, feta, kmax_dealias
    real(8)   :: norm, urms_sq, uscale, u, v, w
    complex(8), allocatable        :: uk(:,:,:,:)    ! (Nf,Nf,Nf,3) host spectral field
    complex(8), allocatable        :: vel_c(:,:,:,:) ! (Nf,Nf,Nf,3) host IFFT output
    complex(8), device, allocatable :: uk_d(:,:,:)   ! one component on device
    complex(8), device, allocatable :: vk_d(:,:,:)   ! IFFT output on device

    if     (kind(id_accuracy) == 2) then; offset = 1
    elseif (kind(id_accuracy) == 4) then; offset = 2
    elseif (kind(id_accuracy) == 8) then; offset = 3
    end if
    Nf = nx - 2*offset
    allocate(uk(Nf,Nf,Nf,3))
    allocate(vel_c(Nf,Nf,Nf,3))
    allocate(uk_d(Nf,Nf,Nf))
    allocate(vk_d(Nf,Nf,Nf))
    ! Step 1: generate random Fourier modes on CPU
    call random_seed(size=seed_size)
    allocate(seed_arr(seed_size))
    seed_arr = 42
    call random_seed(put=seed_arr)
    uk = cmplx(0.d0, 0.d0, kind=8)
    kmax_dealias = dble(Nf) / 3.d0   ! (2/3) * k_Nyquist = (2/3) * (Nf/2) = Nf/3

    do kzi = 1, Nf
      ikz = kzi - 1; if (ikz > Nf/2) ikz = ikz - Nf
      do kyi = 1, Nf
        iky = kyi - 1; if (iky > Nf/2) iky = iky - Nf
        do kxi = 1, Nf
          ikx = kxi - 1; if (ikx > Nf/2) ikx = ikx - Nf
          kmag2 = dble(ikx*ikx + iky*iky + ikz*ikz)
          if (kmag2 < 0.5d0) cycle
          kmag = sqrt(kmag2)
          if (kmag > kmax_dealias) cycle
          ! Pope (2000): E(k) = k^(-5/3) * f_L(kL) * f_eta(k*eta)
          ! amp = sqrt(E(k)/k^2) = k^(-11/6) * sqrt(f_L * f_eta)
          kL   = kmag * pope_L
          keta = kmag * pope_eta
          fL   = (kL / sqrt(kL**2 + pope_cL))**(5.d0/3.d0 + pope_p0)
          feta = exp(-pope_beta * (((keta)**4 + pope_ceta**4)**0.25d0 - pope_ceta))
          amp  = kmag**(-11.d0/6.d0) * sqrt(fL * feta)
          do icomp = 1, 3
            call random_number(r)
            ph = 2.d0 * pi * r(2)
            uk(kxi,kyi,kzi,icomp) = amp * cmplx(r(1)*cos(ph), r(1)*sin(ph), kind=8)
    enddo;enddo;enddo;enddo

    ! Step 2: divergence-free projection in Fourier space
    ! û_i^⊥ = û_i - k_i*(k·û)/|k|^2
    do kzi = 1, Nf
      ikz = kzi - 1; if (ikz > Nf/2) ikz = ikz - Nf
      do kyi = 1, Nf
        iky = kyi - 1; if (iky > Nf/2) iky = iky - Nf
        do kxi = 1, Nf
          ikx = kxi - 1; if (ikx > Nf/2) ikx = ikx - Nf
          kmag2 = dble(ikx*ikx + iky*iky + ikz*ikz)
          if (kmag2 < 0.5d0) cycle
          if (sqrt(kmag2) > kmax_dealias) cycle
          kdu_re = dble(ikx)*real(uk(kxi,kyi,kzi,1)) &
                 + dble(iky)*real(uk(kxi,kyi,kzi,2)) &
                 + dble(ikz)*real(uk(kxi,kyi,kzi,3))
          kdu_im = dble(ikx)*aimag(uk(kxi,kyi,kzi,1)) &
                 + dble(iky)*aimag(uk(kxi,kyi,kzi,2)) &
                 + dble(ikz)*aimag(uk(kxi,kyi,kzi,3))
          uk(kxi,kyi,kzi,1) = uk(kxi,kyi,kzi,1) &
                             - dble(ikx)*cmplx(kdu_re,kdu_im,8)/kmag2
          uk(kxi,kyi,kzi,2) = uk(kxi,kyi,kzi,2) &
                             - dble(iky)*cmplx(kdu_re,kdu_im,8)/kmag2
          uk(kxi,kyi,kzi,3) = uk(kxi,kyi,kzi,3) &
                             - dble(ikz)*cmplx(kdu_re,kdu_im,8)/kmag2
    enddo;enddo;enddo

    ! Step 3: enforce Hermitian symmetry so IFFT produces a real field
    ! uk_sym(k) = 0.5*(uk(k) + conj(uk(-k)))
    ! conjugate index: kxc = mod(Nf-kxi+1, Nf)+1  (1-indexed, size Nf)
    do kzi = 1, Nf
      kzc = mod(Nf-kzi+1, Nf) + 1
      do kyi = 1, Nf
        kyc = mod(Nf-kyi+1, Nf) + 1
        do kxi = 1, Nf
          kxc = mod(Nf-kxi+1, Nf) + 1
          do icomp = 1, 3
            uk(kxi,kyi,kzi,icomp) = 0.5d0 * &
              (uk(kxi,kyi,kzi,icomp) + conjg(uk(kxc,kyc,kzc,icomp)))
    enddo;enddo;enddo;enddo

    ! Step 4: inverse cuFFT per velocity component (GPU)
    istat = cufftPlan3d(plan, Nf, Nf, Nf, CUFFT_Z2Z)
    norm  = 1.d0 / dble(Nf)**3
    do icomp = 1, 3
      uk_d = uk(:,:,:,icomp)
      istat = cufftExecZ2Z(plan, uk_d, vk_d, CUFFT_INVERSE)
      vel_c(:,:,:,icomp) = vk_d
    enddo
    istat = cufftDestroy(plan)

    ! Step 5: compute actual u_rms and rescale to urms, then pack conservative Q
    urms_sq = 0.d0
    do icomp = 1, 3
      do k = 1, Nf
        do j = 1, Nf
          do i = 1, Nf
            urms_sq = urms_sq + (real(vel_c(i,j,k,icomp)) * norm)**2
    enddo;enddo;enddo;enddo
    urms_sq = urms_sq / dble(3 * Nf**3)
    uscale  = urms / sqrt(urms_sq)
    print *, '[DHIT init] Re_lambda =', Re_lambda, &
             '  pope_eta =', pope_eta, '  1/pope_eta =', 1.d0/pope_eta, &
             '  kmax*eta =', dble(Nf)/3.d0 * pope_eta

    do k = 1, Nf
      do j = 1, Nf
        do i = 1, Nf
          u = real(vel_c(i,j,k,1)) * norm * uscale
          v = real(vel_c(i,j,k,2)) * norm * uscale
          w = real(vel_c(i,j,k,3)) * norm * uscale
          Q(i+offset, 1, j+offset, k+offset) = RHO0
          Q(i+offset, 2, j+offset, k+offset) = RHO0 * u
          Q(i+offset, 3, j+offset, k+offset) = RHO0 * v
          Q(i+offset, 4, j+offset, k+offset) = RHO0 * w
          Q(i+offset, 5, j+offset, k+offset) = p0/(gamma-1.d0) &
              + 0.5d0 * RHO0 * (u*u + v*v + w*w)
    enddo;enddo;enddo

    deallocate(uk, vel_c, uk_d, vk_d, seed_arr)
  end subroutine init_spectral_velocity
end module set_init_dhit
