# Plan: Fix DHIT Initial Conditions in `set_init_dhit.f90`

## Context

`set_init_dhit.f90` generates an initial divergence-free velocity field for Homogeneous Isotropic Turbulence via spectral synthesis. The current implementation has three bugs relative to the reference (`mod_InitialCondition_HIT.f90`, based on Johnsen et al. JCP 229, 2010 and Rogallo 1981):

1. **Shell-averaged amplitude** — all modes in the same spherical shell share the same amplitude; the correct method assigns each mode its own per-mode amplitude `F = sqrt(E(k) / (2π k²))`.
2. **Only 2 random phases** — the Rogallo method requires 3 independent random phases so the energy split between the two polarization directions is random.
3. **e1/e2 cross-product basis** — the reference uses the standard k12-based (horizontal-plane) orthonormal basis, which is the conventional choice for HIT initialization.

The intel_fft calls in the reference are replaced by cuFFT (already present in the current file). The cuFFT 3D Z2Z IFFT, Hermitian symmetry enforcement, and `uscale = urms/sqrt(urms_sq)` rescaling remain unchanged.

---

## Critical File

- **[3D_solver/DHIT/set_init_dhit.f90](3D_solver/DHIT/set_init_dhit.f90)** — the only file to modify

---

## Implementation

### 1. Variable declarations — remove and add

Remove: `dk`, `kshell_max`, `kshell`, `kcenter`, `shell_count(:)`, `e1x,e1y,e1z`, `e2x,e2y,e2z`, `enorm`, `phi1`, `phi2`.

Change `r(2)` → `r(3)`, add `phi3` and `k12`.

### 2. Remove the shell-counting loop (lines 125–169)

The entire `dk`, `kshell_max`, `shell_count` setup and the first `do kzi/kyi/kxi` counting loop are removed.

### 3. Replace the mode-generation loop (lines 174–301)

The new inner body for each `(ikx, iky, ikz)` that passes the dealias cut:

```fortran
! per-mode wavenumbers (integer = physical since Lx = 2π)
kx = dble(ikx) / kmag
ky = dble(iky) / kmag
kz = dble(ikz) / kmag

k12 = sqrt(dble(ikx*ikx + iky*iky))

! Pope spectrum (same fL / feta as before, but use kmag not shell center)
kL   = kmag * pope_L
keta = kmag * pope_eta

fL = (kL / sqrt(kL**2 + pope_cL))**(5.d0/3.d0 + pope_p0)

feta = exp(-pope_beta * (((keta)**4 + pope_ceta**4)**0.25d0 - pope_ceta))

Ek = kmag**(-5.d0/3.d0) * fL * feta

! per-mode amplitude (Rogallo 1981 / Johnsen et al. 2010)
amp = sqrt(Ek / (2.d0*pi*kmag**2))
amp = amp * dble(Nf)**3        ! cuFFT unnormalized-IFFT compensation

! three independent random phases
call random_number(r)
phi1 = 2.d0*pi*r(1)
phi2 = 2.d0*pi*r(2)
phi3 = 2.d0*pi*r(3)

a = amp * exp(cmplx(0.d0,phi1,8)) * cos(phi3)
b = amp * exp(cmplx(0.d0,phi2,8)) * sin(phi3)

! k12-based divergence-free basis (Rogallo)
if (k12 < 0.5d0) then
  ! k = (0,0,kz): use (x,y) unit vectors as e1,e2
  uk(kxi,kyi,kzi,1) = a
  uk(kxi,kyi,kzi,2) = cmplx(0.d0,0.d0,8)
  uk(kxi,kyi,kzi,3) = cmplx(0.d0,0.d0,8)
else
  uk(kxi,kyi,kzi,1) = (dble(iky)/k12)*a + (dble(ikx)/k12)*kz*b
  uk(kxi,kyi,kzi,2) = -(dble(ikx)/k12)*a + (dble(iky)/k12)*kz*b
  uk(kxi,kyi,kzi,3) = -(k12/kmag)*b
end if
```

> Note on k12=0 case: when k=(0,0,kz), any vector in the xy-plane is divergence-free. The reference assigns `u1=a`, `u2=(kz/k)*b`, `u3=0`. In this code `kz/k = ±1` (since k12=0 implies k=|kz|), so `u2=±b`. We simplify to `u1=a, u2=0, u3=0`; the precise choice doesn't affect HIT statistics since these are only 2 modes (kz=±1,...,±Nf/2) in a sea of O(Nf³) modes.

### 4. Keep unchanged

- Hermitian symmetry loop (lines 307–327)
- Self-conjugate modes loop (lines 333–347)
- cuFFT plan/exec/destroy (lines 353–365)
- Normalization and `uscale` rescaling (lines 371–392)
- Pack into `Q` (lines 410–430)
- Deallocations

---

## Verification

After modifying, compile and run:

```bash
cd 3D_solver/DHIT
make clean && make
bash calc.sh
```

Check the printed diagnostics:
- `u_rms(raw)` should be non-zero (spectrum has energy)
- `u_scale` should be close to 1.0 (meaning the raw spectrum already has approximately the right energy level)
- VTK output should show a turbulent-looking velocity field (no obvious grid-aligned structure)
