# Plan: DHIT Test Case with cuFFT Spectral Initialization

## Context

Add a new test case `3D_solver/DHIT/` (Decaying Homogeneous Isotropic Turbulence) to the OUxSBLI solver. The case is weakly compressible (Ma=0.1), triply periodic, with no external forcing. Unlike the deterministic Taylor-Green vortex IC in NSTGV, DHIT requires a random velocity field drawn from a prescribed energy spectrum E(k) = A k^4 exp[-2(k/kp)^2] with kp=4. The initialization must be done spectrally via cuFFT (Z2Z double-precision complex-to-complex transform) with a divergence-free projection enforced in Fourier space and Hermitian symmetry enforced before the IFFT.

---

## Files to Create

```
3D_solver/DHIT/
  mod_globals.f90      — parameters (NSTGV copy, M0=0.1, kp=4 added)
  set_init_dhit.f90    — new module: cuFFT spectral IC
  set.f90              — grid/IC/BC (NSTGV copy, set_init replaced)
  Makefile             — NSTGV copy + set_init_dhit.o + -lcufft
  calc.sh              — same as NSTGV
```

---

## File Details

### 1. `3D_solver/DHIT/mod_globals.f90`

Copy `3D_solver/NSTGV/mod_globals.f90`. Changes:

| Line | Original (NSTGV) | DHIT value | Notes |
|------|-----------------|------------|-------|
| `id_visc` | `integer(4)` | `integer(4)` | NS, unchanged |
| `id_scheme` | `real(2)` | `integer(2)` | Switch to KEEP (better for smooth vortical flow) |
| `L0` | `1.524d-3` | `1.0d0` | Unit reference length; domain is exactly [0, 2π] |
| `M0` | `1.25d0` | `0.1d0` | Weakly compressible |
| `dt` | `CFL*(Lx/(nx-1))/V0` | `CFL*(Lx/dble(nx-1))*M0/V0` | Acoustic-limited; c_s = V0/M0 dominates |
| Add | — | `real(8), parameter :: kp = 4.d0` | Peak wavenumber for E(k) |

Keep all `id_bc_*=.false.`, `nx=ny=nz=66`, `Lx=Ly=Lz=2π*L0`, `Re=1600`, `T=294.4K`, GPU thread configs, and time-integration flags unchanged. With `integer(kind=8) :: id_accuracy`, offset=3 and interior Nf=60 per direction.

---

### 2. `3D_solver/DHIT/set_init_dhit.f90`

New module encapsulating the spectral initialization. Called once from `set_init` in `set.f90`.

```fortran
module set_init_dhit
  use cufft
  implicit none
contains

  subroutine init_spectral_velocity(nx, ny, nz, Q)
    use mod_globals, only : id_accuracy, gamma, R, RHO0, p0, V0, kp, pi
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(out) :: Q(nx,5,ny,nz)

    integer  :: offset, Nf, plan, istat, icomp, seed_arr(1)
    integer  :: kxi, kyi, kzi
    integer(8) :: ikx, iky, ikz
    real(8)  :: kmag2, kmag, amp, r(2), ph, kdu_re, kdu_im, norm, urms_sq, uscale
    integer  :: kxc, kyc, kzc, i, j, k

    complex(8), allocatable        :: uk(:,:,:,:)    ! (Nf,Nf,Nf,3) host spectral field
    complex(8), allocatable        :: vel_c(:,:,:,:) ! (Nf,Nf,Nf,3) host physical (complex)
    complex(8), device, allocatable :: uk_d(:,:,:)   ! one component, device
    complex(8), device, allocatable :: vk_d(:,:,:)   ! IFFT output, device

    ! Determine ghost-cell offset from id_accuracy kind
    if     (kind(id_accuracy) == 2) then; offset = 1
    elseif (kind(id_accuracy) == 4) then; offset = 2
    elseif (kind(id_accuracy) == 8) then; offset = 3
    end if
    Nf = nx - 2*offset   ! interior (FFT) size per direction

    allocate(uk(Nf,Nf,Nf,3), vel_c(Nf,Nf,Nf,3))
    allocate(uk_d(Nf,Nf,Nf), vk_d(Nf,Nf,Nf))

    ! ── Step 1: generate random spectral field (CPU) ─────────────────────────
    seed_arr(1) = 42
    call random_seed(put=seed_arr)
    uk = cmplx(0.d0, 0.d0, kind=8)

    do kzi = 1, Nf
      ikz = kzi - 1; if (ikz > Nf/2) ikz = ikz - Nf
      do kyi = 1, Nf
        iky = kyi - 1; if (iky > Nf/2) iky = iky - Nf
        do kxi = 1, Nf
          ikx = kxi - 1; if (ikx > Nf/2) ikx = ikx - Nf
          kmag2 = dble(ikx*ikx + iky*iky + ikz*ikz)
          if (kmag2 < 0.5d0) cycle
          kmag = sqrt(kmag2)
          ! E(k) = k^4 * exp(-2*(k/kp)^2); amplitude ∝ sqrt(E(k)/k^2)
          amp = kmag**2 * exp(-dble(kmag/kp)**2)   ! = sqrt(k^4 * exp(-2(k/kp)^2) / k^2)
          do icomp = 1, 3
            call random_number(r)              ! r(1): magnitude, r(2): phase
            ph = 2.d0 * pi * r(2)
            uk(kxi,kyi,kzi,icomp) = amp * cmplx(r(1)*cos(ph), r(1)*sin(ph), kind=8)
          end do
        end do
      end do
    end do

    ! ── Step 2: divergence-free projection (CPU, Fourier space) ──────────────
    do kzi = 1, Nf
      ikz = kzi - 1; if (ikz > Nf/2) ikz = ikz - Nf
      do kyi = 1, Nf
        iky = kyi - 1; if (iky > Nf/2) iky = iky - Nf
        do kxi = 1, Nf
          ikx = kxi - 1; if (ikx > Nf/2) ikx = ikx - Nf
          kmag2 = dble(ikx*ikx + iky*iky + ikz*ikz)
          if (kmag2 < 0.5d0) cycle
          ! k·û (complex dot product: k is real, û is complex)
          kdu_re = dble(ikx)*real(uk(kxi,kyi,kzi,1)) &
                 + dble(iky)*real(uk(kxi,kyi,kzi,2)) &
                 + dble(ikz)*real(uk(kxi,kyi,kzi,3))
          kdu_im = dble(ikx)*aimag(uk(kxi,kyi,kzi,1)) &
                 + dble(iky)*aimag(uk(kxi,kyi,kzi,2)) &
                 + dble(ikz)*aimag(uk(kxi,kyi,kzi,3))
          ! û_i^⊥ = û_i - k_i*(k·û)/|k|^2
          uk(kxi,kyi,kzi,1) = uk(kxi,kyi,kzi,1) - dble(ikx)*cmplx(kdu_re,kdu_im,8)/kmag2
          uk(kxi,kyi,kzi,2) = uk(kxi,kyi,kzi,2) - dble(iky)*cmplx(kdu_re,kdu_im,8)/kmag2
          uk(kxi,kyi,kzi,3) = uk(kxi,kyi,kzi,3) - dble(ikz)*cmplx(kdu_re,kdu_im,8)/kmag2
        end do
      end do
    end do

    ! ── Step 3: enforce Hermitian symmetry (CPU) ──────────────────────────────
    ! Symmetrize: uk_sym(k) = 0.5*(uk(k) + conj(uk(-k)))
    ! For 1-indexed size-Nf array: conjugate index = mod(Nf-idx+1, Nf)+1
    do kzi = 1, Nf
      kzc = mod(Nf-kzi+1, Nf) + 1
      do kyi = 1, Nf
        kyc = mod(Nf-kyi+1, Nf) + 1
        do kxi = 1, Nf
          kxc = mod(Nf-kxi+1, Nf) + 1
          do icomp = 1, 3
            uk(kxi,kyi,kzi,icomp) = 0.5d0 * &
              (uk(kxi,kyi,kzi,icomp) + conjg(uk(kxc,kyc,kzc,icomp)))
          end do
        end do
      end do
    end do

    ! ── Step 4: inverse cuFFT per component (GPU) ────────────────────────────
    ! For cubic Nf^3 array, dimension order is unambiguous
    istat = cufftPlan3d(plan, Nf, Nf, Nf, CUFFT_Z2Z)
    norm  = 1.d0 / dble(Nf)**3

    do icomp = 1, 3
      uk_d = uk(:,:,:,icomp)                           ! host → device
      istat = cufftExecZ2Z(plan, uk_d, vk_d, CUFFT_INVERSE)
      vel_c(:,:,:,icomp) = vk_d                         ! device → host
    end do
    istat = cufftDestroy(plan)

    ! ── Step 5: normalize, rescale to target u_rms = V0, pack Q (CPU) ────────
    urms_sq = 0.d0
    do icomp = 1, 3
      urms_sq = urms_sq + sum(real(vel_c(:,:,:,icomp))**2)
    end do
    urms_sq = urms_sq * norm**2 / dble(3*Nf**3)   ! mean squared velocity per component
    uscale  = V0 / sqrt(urms_sq)                   ! rescale so u_rms = V0

    do k = 1, Nf
      do j = 1, Nf
        do i = 1, Nf
          Q(i+offset,1,j+offset,k+offset) = RHO0
          Q(i+offset,2,j+offset,k+offset) = RHO0 * real(vel_c(i,j,k,1)) * norm * uscale
          Q(i+offset,3,j+offset,k+offset) = RHO0 * real(vel_c(i,j,k,2)) * norm * uscale
          Q(i+offset,4,j+offset,k+offset) = RHO0 * real(vel_c(i,j,k,3)) * norm * uscale
          Q(i+offset,5,j+offset,k+offset) = p0/(gamma-1.d0) + 0.5d0/RHO0 * ( &
              Q(i+offset,2,j+offset,k+offset)**2 &
            + Q(i+offset,3,j+offset,k+offset)**2 &
            + Q(i+offset,4,j+offset,k+offset)**2 )
        end do
      end do
    end do

    deallocate(uk, vel_c, uk_d, vk_d)
  end subroutine init_spectral_velocity

end module set_init_dhit
```

---

### 3. `3D_solver/DHIT/set.f90`

Copy `3D_solver/NSTGV/set.f90`. Changes:
- Add `use set_init_dhit` to the module's `use` list (line 1–5 of the module)
- Replace the body of `set_init` (lines 21–46) with:

```fortran
  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    use mod_globals, only : id_accuracy
    use set_init_dhit
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,5,ny,nz)
    call init_spectral_velocity(nx, ny, nz, Q)
    call set_bc_cyclic(id_accuracy, nx, ny, nz, Q)
  end subroutine set_init
```

`set_grid`, `set_bc`, and `set_bc_mut` are unchanged.

---

### 4. `3D_solver/DHIT/Makefile`

Copy `3D_solver/NSTGV/Makefile`. Changes:

1. **Add `set_init_dhit.o`** to `OBJ` list, immediately before `set.o`
2. **Link line** — append `-lcufft`:
   ```makefile
   a.out: $(OBJ)
       $(FC) -cuda -acc -fast -Minfo -gpu=rdc,lto -lnvhpcwrapnvtx -lcufft -I./ $^
   ```
3. **Add dependency rules**:
   ```makefile
   set_init_dhit.o: mod_globals.mod
   set.o: mod_globals.mod mod_constant.mod calc_para.mod set_bc_common.mod set_coordinate.mod set_init_dhit.mod
   ```
   (The `set.o` line replaces the existing NSTGV `set.o` rule.)

---

### 5. `3D_solver/DHIT/calc.sh`

Identical to `3D_solver/NSTGV/calc.sh`:
```bash
nohup mpirun -n 2 a.out &
cp ./set.f90 ./data
cp ./mod_globals.f90 ./data
```

Create `data/` subdirectory before first run: `mkdir -p data`.

---

## MPI Initialization Semantics

With 2 MPI ranks and 1D x-decomposition, each rank allocates the full `Q(nx,5,ny,nz)` array and both ranks independently call `init_spectral_velocity` with the **same fixed seed (42)**. Both ranks produce identical `Q` arrays. The time-stepping MPI decomposition (`calc_para.f90`) then handles which x-strip each rank actively updates. This mirrors the NSTGV pattern where both ranks evaluate the analytical IC at all grid points.

---

## Critical File References

| File | Purpose |
|------|---------|
| [3D_solver/NSTGV/mod_globals.f90](3D_solver/NSTGV/mod_globals.f90) | Template for mod_globals |
| [3D_solver/NSTGV/set.f90](3D_solver/NSTGV/set.f90) | Template for set.f90; set_init body replaced |
| [3D_solver/NSTGV/Makefile](3D_solver/NSTGV/Makefile) | Template for Makefile |
| [3D_solver/src/set_bc_common.f90](3D_solver/src/set_bc_common.f90) | `set_bc_cyclic` (unchanged; reused) |
| [src/set_coordinate.f90](src/set_coordinate.f90) | `set_grid_cyclic` (unchanged; reused) |

---

## Verification

```bash
cd 3D_solver/DHIT
mkdir -p data
make clean && make          # should compile cleanly; cufft module resolves from NVHPC SDK
bash calc.sh                # launches 2 MPI ranks

# In ParaView: open the VTK output in data/; check velocity field looks isotropic
# and spatially random (not sinusoidal like TGV)

# Quick sanity check at t=0: mean velocity should be ~0, u_rms ≈ V0 per component,
# no large-scale structure in any single direction
```
