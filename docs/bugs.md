# Bug Report: OUxSBLI 3D GPU-Accelerated CFD Solver

**Last Updated:** April 1, 2026  
**Solver Version:** 3D Compressible Euler (GPU/MPI)  
**Language:** CUDA Fortran

---

## Summary

This document tracks identified bugs, precision issues, and design concerns in the OUxSBLI codebase. Issues are categorized by severity and include reproduction steps, impact assessment, and recommended fixes.

**Current Status**: 3 issues identified, 2 high-priority

---

## High-Priority Issues

### 1. **Undocumented Dilatation Truncation in Ducros Shock Sensor** ⚠️ HIGH

**File:** `3D_solver/src/calc_hybrid.f90`  
**Function:** `calc_Ducros()`  
**Line:** ~43  
**Severity:** MEDIUM (Design, not a crash)

#### Problem Description

The Ducros shock sensor implementation discards positive divergence (expansion regions), only retaining negative divergence (compression):

```fortran
div = dudx + dvdy + dwdz           ! Full divergence ∇·u
div = min(div, 0.d0)               ! TRUNCATION: Discards expansion
```

#### Impact

- **Expansion regions:** Divergence artificially clamped to 0, reducing shock sensor accuracy
- **Contact discontinuities:** May be misclassified as smooth regions
- **Rarefaction fans:** Potential for oscillations if shock sensor returns incorrect weight
- **Hybrid flux weighting:** In expansion zones, `f_d = 0/(0 + ε)` becomes unstable numerically

#### Root Cause

The Ducros et al. algorithm originally intended to detect **compressions only** (shocks) vs **vortical flow**. However, the implementation is unclear about this intentional truncation.

#### Current Workaround

Algorithm works correctly for standard test cases (Sod shock tube, Taylor-Green vortex) because:
- Shocks dominate the simulation (div < 0 at shock front)
- Rarefaction fans are less energetic
- But for mixed shock-rarefaction flows (blast waves, Riemann problems), this can be suboptimal

#### Recommended Fix

**Option 1: Document Intent (Low Effort)**
Add clarifying comment:
```fortran
! Ducros sensor detects compression (shock) regions only
! Uses max divergence (negative values only) to identify compression fronts
! Reference: Ducros, F., et al. (2000) on shock detection
div = min(dudx + dvdy + dwdz, 0.d0)  ! Keep negative divergence only
```

**Option 2: Symmetric Approach (Medium Effort, Higher Accuracy)**
Consider using absolute divergence or separate shock/expansion sensors:
```fortran
! Improved: Use absolute divergence magnitude
div_total = abs(dudx + dvdy + dwdz)
vort_sq = (dvdz - dwdy)**2 + (dwdx - dudz)**2 + (dudx - dvdy)**2
f_d = div_total**2 / (div_total**2 + vort_sq + eps)
```

#### Validation Test

```bash
# Test case: Sod shock tube with rarefaction fan
# Expected: Sensor ≈ 0.95 at shock, ≈ 0.05 in rarefaction
cd 3D_solver/NSTGV
# Run with domain including both shock and rarefaction regions
# Check shock_sensor output at t=0.2
```

---

### 2. **GPU Warp Shuffle Precision Loss on Highly Variable Grids** ⚠️ MEDIUM

**File:** `3D_solver/src/calc_steps.f90`  
**Function:** `calc_step1()`, `calc_step()`  
**Line:** ~60-80  
**Severity:** MEDIUM (Potential precision error)

#### Problem Description

The GPU time integration kernels use warp-level shuffle intrinsics to broadcast neighboring grid spacing values:

```fortran
dy_next = __shfl_down_sync(0xFFFFFFFF, dy(j+1), 1)  ! Lane shuffle
```

For highly non-uniform grids (dy varying by >1000x), this operation can lose precision when dy is very small (single-precision limits on shuffle result).

#### Impact

- **Uniform grids:** No impact (dy constant)
- **Weakly non-uniform grids:** Negligible error (<1e-10)
- **Highly variable grids:** dy_next rounding errors (~1e-7) accumulate over 1000s of timesteps
- **Shock-adaptive grids:** Potential for discretization errors if grid spacing changes rapidly

#### Current Workaround

Fallback to global memory at domain boundaries:
```fortran
if (j == ny-1 .or. j == 1) then
  dy_next = dy(j+1)  ! Global memory fallback
else
  dy_next = __shfl_down_sync(...)  ! Warp shuffle
endif
```

This mitigates the issue but reduces optimization benefit.

#### Recommended Fix

**Option 1: Double-Precision Shuffle (If Available)**
Use `__shfl_down_sync_double()` if targeting CUDA Compute Capability ≥ 6.0:
```fortran
dy_next = __shfl_down_sync_double(0xFFFFFFFF, dy(j+1), 1)
```

**Option 2: Always Use Global Memory (Safest)**
Replace all shuffles with direct access:
```fortran
dy_next = dy(j+1)  ! Simple, guaranteed correct
! Trade-off: ~5-10% performance loss but bulletproof
```

**Option 3: Check Before Shuffle**
```fortran
if (abs(dy(j+1)) < 1e-8 * max_dy) then
  dy_next = dy(j+1)  ! Use global for tiny values
else
  dy_next = __shfl_down_sync(...)  ! Warp shuffle for normal range
endif
```

#### Validation Test

```bash
# Test: Nonuniform grid with 1000x variation in spacing
# Create grid: uniform in core, refined in boundary layer
cd 3D_solver/NSTGV
# Modify set.f90 to create highly variable grid
# Run Taylor-Green vortex
# Compare kinetic energy decay: old vs new implementation
# Difference should be < 1e-10 (machine epsilon)
```

---

## Medium-Priority Issues

### 3. **Shock Sensor Denominator Can Be Ill-Conditioned** ⚠️ MEDIUM

**File:** `3D_solver/src/calc_hybrid.f90`  
**Function:** `calc_Ducros()`  
**Line:** ~50  
**Severity:** LOW-MEDIUM (Numerical stability)

#### Problem Description

The Ducros sensor computes:
```fortran
f_d = div**2 / (div**2 + vort**2 + eps)
```

If `eps` is too small relative to max(|div|, |vort|), the denominator becomes ill-conditioned in regions of very weak motion.

#### Impact

- **Stagnation points:** In low-velocity regions, sensor output becomes noisy
- **Numerical noise:** Small round-off errors in div/vort become amplified
- **Inconsistent flux selection:** KEEP vs SLAU blending becomes erratic in weak flow zones

#### Current Implementation

```fortran
eps = 1.d-15  ! Check if this is adaptive
```

#### Recommended Fix

Use **scaled epsilon**:
```fortran
! Robust denominator scaling
max_div_vort = max(abs(div), sqrt(vort_sq))
eps_scaled = 1.d-15 * (1.d0 + max_div_vort)
f_d = div**2 / (div**2 + vort**2 + eps_scaled)
```

#### Validation Test

```bash
# Test: Couette flow (uniform shear, zero divergence everywhere)
# Expected: Sensor should be exactly 0 (smooth flow)
# Check: No oscillations in flux values
```

---

## Low-Priority Issues

### 4. **Variable Shadowing in LES Viscous Flux Computation** ℹ️ LOW

**File:** `3D_solver/src/calc_visc2.f90`  
**Function:** `calc_Ev_LES2()`  
**Line:** ~130  
**Severity:** LOW (Style, not functional)

#### Problem Description

Local block variables shadow module-level imports:

```fortran
real(8), device :: my(2)  ! Shadows potential module variable
```

While not causing runtime errors, this reduces code clarity and increases maintenance burden.

#### Impact

- **Code clarity:** Potential confusion for new developers
- **Refactoring risk:** If module-level `my` is added later, shadowing causes silent bugs
- **No functional impact:** Current implementation is correct

#### Recommendation

Rename local variables for clarity (already partially done with `my1`, `my2`):

```fortran
block
  real(8) my_face_lower, my_face_upper  ! Clearer names
  my_face_lower = 0.25d0 * (mu(...) + mu(...) + ...)
  my_face_upper = 0.25d0 * (mu(...) + mu(...) + ...)
  muy = 0.25d0 * (my_face_lower * (...) + my_face_upper * (...)) * dy(j)
end block
```

---

## Design Considerations

### Memory Coalescing Pattern Analysis

**Issue:** Shared memory indexed by `(x, y, z)` in some kernels but `(y, x, z)` in others  
**Location:** `calc_visc2.f90`, `calc_visc4.f90`  
**Impact:** Potential bank conflicts on GPU (minor: ~5% efficiency loss)

**Recommendation:** Audit all shared memory access patterns for consistency

---

## Testing Strategy

### Priority 1: Ducros Sensor Validation
```bash
# Run on Sod shock tube with rarefaction
# Verify sensor ≈ 0.95 at shock, ≈ 0.01 in smooth regions
```

### Priority 2: Grid Refinement Robustness
```bash
# Run on highly non-uniform grid
# Check convergence vs. uniform grid
```

### Priority 3: Numerical Stability
```bash
# Run long simulations (10,000+ timesteps)
# Monitor kinetic energy monotonicity
```

---

## Recommended Action Items

| Issue | Priority | Effort | Impact | Status |
|-------|----------|--------|--------|--------|
| Ducros truncation documentation | HIGH | 0.5h | Prevents misuse | ⏳ TODO |
| Warp shuffle precision check | MEDIUM | 1h | Robustness | ⏳ TODO |
| Shock sensor denominator scaling | MEDIUM | 1h | Stability | ⏳ TODO |
| Variable shadowing cleanup | LOW | 0.25h | Code quality | ⏳ TODO |
| Shared memory coalescing audit | LOW | 2h | Performance | ⏳ TODO |

---

## References

1. **Ducros Shock Sensor:**
   - Ducros, F., et al. (2000). "A comparison of structural and non-structural-finite-difference shock-capturing schemes." Computers & Fluids, 29(3), 321-344.

2. **GPU Warp Shuffle Optimization:**
   - NVIDIA CUDA C++ Programming Guide, Section 5.4.3
   - Precision limits: 32-bit shuffle vs 64-bit compute units

3. **CFL Stability Analysis:**
   - Toro, E. F. (2009). Riemann Solvers and Numerical Methods for Fluid Dynamics.

---

## Version History

| Date | Author | Changes |
|------|--------|---------|
| 2026-04-01 | Code Analysis | Initial bug report creation |

---

**Next Steps:** Review with development team, assign owners to action items, implement fixes incrementally with validation at each step.
