# Plan: Transpose E, F, G Flux Arrays to AoS Layout

## Context

E, F, G are the x-, y-, z-direction flux arrays written by convective and viscous kernels and read by the time-integration (RK) kernels. Their current layout is `E(nx-1, 5, ny-2, nz-2)` — spatial index first, component index second (SoA-like). In AoS the 5 conserved-variable components become the first (fastest-varying) dimension: `E(5, nx-1, ny-2, nz-2)`. This makes all 5 writes per spatial point contiguous in memory, which is the dominant access pattern in the convective kernels (KEEP, SLAU, Hybrid each write all 5 components for one point per thread).

## Layout Change Summary

| | Current (SoA) | Target (AoS) |
|---|---|---|
| E | `E(nx-1, 5, ny-2, nz-2)` | `E(5, nx-1, ny-2, nz-2)` |
| F | `F(nx-2, 5, ny-1, nz-2)` | `F(5, nx-2, ny-1, nz-2)` |
| G | `G(nx-2, 5, ny-2, nz-1)` | `G(5, nx-2, ny-2, nz-1)` |

Index transformation: `E(i, m, j, k)` → `E(m, i, j, k)` (swap first two indices everywhere).

Array section: `E(i, :, j, k)` → `E(:, i, j, k)`.

## Files to Modify

14 files in `3D_solver/src/`:

### 1. `preprocess.f90`
- Line 48: allocation
  - `E(nx-1,5,ny-2,nz-2)` → `E(5,nx-1,ny-2,nz-2)`
  - `F(nx-2,5,ny-1,nz-2)` → `F(5,nx-2,ny-1,nz-2)`
  - `G(nx-2,5,ny-2,nz-1)` → `G(5,nx-2,ny-2,nz-1)`

### 2. `calc_flux_base.f90`
- ~9 subroutine signatures, each declaring E, F, G — shape change only.
- No direct element indexing in this file.

### 3. `calc_steps.f90`
- 4 subroutine signatures: shape change.
- 4 flux-divergence loops (`do l = 1, 5`):
  - `E(i,l,j,k)` → `E(l,i,j,k)`, `E(i+1,l,j,k)` → `E(l,i+1,j,k)`
  - `F(i,l,j,k)` → `F(l,i,j,k)`, `F(i,l,j+1,k)` → `F(l,i,j+1,k)`
  - `G(i,l,j,k)` → `G(l,i,j,k)`, `G(i,l,j,k+1)` → `G(l,i,j,k+1)`

### 4. `calc_keep_kernel.f90`
- Multiple subroutine signatures: shape change.
- Array section writes: `E(i,:,j-1,k-1)` → `E(:,i,j-1,k-1)`, similarly for F and G.

### 5. `calc_keep_kernel_internal.f90`
- Subroutine signatures: shape change.
- Array section writes: same pattern as above.

### 6. `calc_slau_kernel.f90`
- Multiple subroutine signatures: shape change.
- Individual component writes:
  - `E(i,1,j-1,k-1), E(i,2,j-1,k-1), E(i,3,j-1,k-1), E(i,4,j-1,k-1), E(i,5,j-1,k-1)` → `E(1,i,j-1,k-1), ...`
  - Similar for F, G.

### 7. `calc_slau_kernel_internal.f90`
- Subroutine signatures + same index swap as calc_slau_kernel.f90.

### 8. `calc_hybrid_kernel.f90`
- Multiple subroutine signatures: shape change.
- Mix of array sections and individual component writes: same swap pattern.

### 9. `calc_hybrid_kernel_internal.f90`
- Subroutine signatures + same swap.

### 10. `calc_roe_kernel.f90`
- Multiple subroutine signatures + individual component writes: same swap (keep consistent even though Roe is unused).

### 11. `calc_visc2.f90`
- 6 subroutine signatures: shape change.
- Accumulation writes `E(i,2,j-1,k-1) = E(i,2,j-1,k-1) - ...` → `E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - ...`
- Same for components 3, 4, 5; and for F and G.

### 12. `calc_visc4.f90`
- 6 subroutine signatures: shape change.
- Same accumulation pattern as visc2.

### 13. `calc_visc4_internal.f90`
- Subroutine signatures + same accumulation swap.

### 14. `calc_visc4_les_internal.f90`
- Subroutine signatures + same accumulation swap.

## Index Transformation Rules (for implementation)

Apply these text substitutions consistently across all 14 files:

**Declarations:**
- `E(nx-1,5,ny-2,nz-2)` → `E(5,nx-1,ny-2,nz-2)`
- `F(nx-2,5,ny-1,nz-2)` → `F(5,nx-2,ny-1,nz-2)`
- `G(nx-2,5,ny-2,nz-1)` → `G(5,nx-2,ny-2,nz-1)`

**Array section writes (KEEP/Hybrid):**
- `E(i,:,j-1,k-1)` → `E(:,i,j-1,k-1)`
- `E(i-1,:,j-1,k-1)` → `E(:,i-1,j-1,k-1)` (if any)
- `F(i-1,:,j,k-1)` → `F(:,i-1,j,k-1)`
- `G(i-1,:,j-1,k)` → `G(:,i-1,j-1,k)`

**Element access (SLAU/Roe/Hybrid/visc/steps):**
- `E(i,N,j,k)` → `E(N,i,j,k)` for N in {1,2,3,4,5,l}
- `E(i+1,N,j,k)` → `E(N,i+1,j,k)`
- `F(i,N,j,k)` → `F(N,i,j,k)`
- `F(i,N,j+1,k)` → `F(N,i,j+1,k)`
- `F(i-1,N,j,k)` → `F(N,i-1,j,k)`
- `G(i,N,j,k)` → `G(N,i,j,k)`
- `G(i,N,j,k+1)` → `G(N,i,j,k+1)`
- `G(i-1,N,j-1,k)` → `G(N,i-1,j-1,k)`

## Out of Scope

- `Rs(nx-2,5,ny-2,nz-2)` in `calc_time_dev.f90` — not part of E/F/G; leave unchanged.
- `Q(nx,5,ny,nz)` and `ruvwp(5,nx,ny,nz)` — unrelated arrays; leave unchanged.

## Verification

1. `cd 3D_solver/NSTGV && make clean && make` — must compile without errors.
2. `bash calc.sh` — run simulation; output VTK files must be produced.
3. Compare physical quantities (e.g., kinetic energy, divergence) with a pre-change reference to confirm numerical identity.
