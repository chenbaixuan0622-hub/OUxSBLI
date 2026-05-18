# Bug & Performance Survey: OUxSBLI optim2 Branch

## Context

Full survey of the `optim2` branch after the original NSTGV result bug was fixed. Goal: find remaining bugs and performance bottlenecks in `3D_solver/src/`. Per `CLAUDE.md`, the main targets are KEEP, SLAU, and Hybrid schemes in `calc_flux_base.f90`, `calc_steps.f90`, and the `calc_*_kernel.f90` files.

---

## Already-Fixed Bugs (this conversation)

| Bug | File | Fix |
|-----|------|-----|
| `coef3` dropped in TVD RK3 Stage 3 | `calc_steps.f90` | Added `coef3 *` to three assignments in `calc_step2_3` |
| `eps = 1.0_sp-12` evaluated as −11.0 | `calc_hybrid.f90:23` | Changed to `1.0e-12_sp` |
| Stale imports `accuracy, offset` | `calc_hybrid.f90:5` | Removed from `use mod_globals, only:` |
| `fdz` uninitialized in `calc_roe_z4` | `calc_roe_kernel.f90` | Added `fdz = 0.5_sp * (sensor(i,j,k) + sensor(i,j,k+1))` before `if` block |
| `threshold` type mismatch with `sp` | all 6 `*/mod_globals.f90` | Changed to `real(sp), parameter :: threshold = 0.4_sp` |

---

## Remaining Confirmed Bugs

### Bug 1: `pure attributes(global)` on `calc_Ducros`

**File:** `3D_solver/src/calc_hybrid.f90:12`

```fortran
! Current (wrong):
pure attributes(global) subroutine calc_Ducros(nx, ny, nz, dx, dy, dz, Q, fd)

! Fix (remove `pure`):
attributes(global) subroutine calc_Ducros(nx, ny, nz, dx, dy, dz, Q, fd)
```

A CUDA `global` kernel writes to device memory (`fd(i,j,k) = ...`), which is a side effect. Fortran's `pure` guarantees no side effects — the two attributes are contradictory. The NVHPC compiler may silently accept it, but `pure` must be removed to be correct and to allow the compiler to apply normal GPU kernel optimizations.

### Bug 2: Memory leak in `calc_Gaussian_filter_x`

**File:** `3D_solver/src/set_init_common.f90:36`

```fortran
! Current (leaks wg):
deallocate(phi_tmp)

! Fix:
deallocate(phi_tmp, wg)
```

`wg` is allocated at line 18 but not freed. The y-variant (line 69) and z-variant (line 102) correctly deallocate both. This is a CPU-side host memory leak.

---

## Performance Bottlenecks

### Bottleneck 1: Warp divergence from stencil-size branching (highest impact)

**Files:** all x6/y6/z6 kernels in `calc_keep_kernel.f90`, `calc_slau_kernel.f90`, `calc_hybrid_kernel.f90`

Each 6-point kernel has a 3-way branch per thread:

```fortran
if (3 <= i .and. i <= nx-3) then      ! 6-point stencil (interior)
  ...
elseif (2 <= i .and. i <= nx-2) then  ! 4-point stencil (near-boundary)
  ...
else                                    ! 2-point stencil (boundary)
  ...
endif
```

Threads in the same warp hit different branches depending on their global index. This serializes execution within a warp. Since all x/y/z flux kernels share this pattern, it affects every flux computation.

**Mitigation:** Split the domain into an interior launch (guaranteed 6-point, no branches) and a thin-shell boundary launch. Interior blocks can run with zero divergence. However, this is a significant restructuring — evaluate ROI with profiling before implementing.

### Bottleneck 2: Warp divergence from KEEP/SLAU scheme switching (hybrid kernels)

**Files:** `calc_hybrid_kernel.f90` x6/y6/z6

Inside the stencil branch, there is a second diverging branch:

```fortran
if (fdx <= threshold) then   ! KEEP path
  ...
else                          ! SLAU path (via delta6/delta4)
  ...
endif
```

This stacks on top of Bottleneck 1, causing up to 6 distinct execution paths in a single warp. In shock-free regions (NSTGV, KHI) the sensor is uniformly low, so all threads take the KEEP path — no divergence. In mixed regions (SBLI, TBL) divergence is significant.

### Bottleneck 3: SLAU x/y/z kernel sensor dual-use (`fdx` reused as wiggle_detector)

**File:** `calc_slau_kernel.f90` x6/y6/z6, lines ~73-98

The Ducros shock sensor is loaded into `fdx`, passed to `delta6`/`delta4` for MUSCL reconstruction, then overwritten with `wiggle_detector(p(...))` before being passed to `SLAU`. This is intentional (documented by comment `!< 1st use: shock sensor, 2nd use: wiggle detector`), but adds extra computation (wiggle_detector reads pressure from shared memory again). No fix needed unless profiling shows it's a bottleneck.

---

## False Positives (not bugs)

| Claim | Reality |
|-------|---------|
| RK3 Stage 3 coef wrong (coef1=2, coef2=1) | Correct — Qin=Q^(2), Qout=Q^n; formula gives `(2/3)Q^(2)+(1/3)Q^n - (2/3)R` ✓ |
| Block scope violation in hybrid kernels | Correct — rhol/ul/etc. copied into shared memory before block ends (lines 133-137), SLAU uses shared arrays after block |
| Ducros corner cells unset | Not accessed — x-flux uses j∈[2,ny-1], y-flux uses i∈[2,nx-1], so corners never reached simultaneously |
| id_accuracy4/id_accuracy2 uninitialized | Typed dispatch tokens — KIND matters, VALUE is irrelevant; x4 uses `id_accuracy` directly (which has kind=4 for x4 variant) |

---

## Fix Plan

### Step 1: Remove `pure` from `calc_Ducros` (5 min)

```fortran
! calc_hybrid.f90:12
! Change:
pure attributes(global) subroutine calc_Ducros(...)
! To:
attributes(global) subroutine calc_Ducros(...)
```

### Step 2: Fix memory leak in `calc_Gaussian_filter_x` (2 min)

```fortran
! set_init_common.f90:36
! Change:
deallocate(phi_tmp)
! To:
deallocate(phi_tmp, wg)
```

### Step 3 (optional, high effort): Kernel splitting for warp divergence

Only worth implementing after profiling with `nsys`/`ncu` to confirm divergence overhead exceeds the cost of extra kernel launch overhead.

---

## Verification

```bash
cd 3D_solver/NSTGV
make clean && make
bash calc.sh
# Compare kinetic energy output against optim branch reference
```

For TBL/SBLI (shock cases), compare wall pressure or skin friction against reference to verify sensor/hybrid behavior is unchanged.
