# Plan: Optimize GPU bottleneck kernels via code restructuring

## Context

Two kernels dominate runtime on the 513³ NSTGV case (H100 GPU, TVD RK3):

| Kernel | Per-call | Calls/step | Total/step |
|---|---|---|---|
| `calc_step2_3` | 12.2 ms | 2 | **24.4 ms** |
| `calc_Gv4_in`  | 9.74 ms | 3 | **29.2 ms** |

Both are register-pressure-limited. The goal is to reduce register usage via code restructuring (not `maxregcount`) so the compiler can schedule more blocks per SM.

Register thresholds (H100: 65,536 regs/SM, 4 warps/block):
- `calc_step2_3`: 42 regs → 1,536 regs/warp → 10 blocks → 62.5% occupancy.
  At ≤40 regs → 1,280 regs/warp → **12 blocks → 75% occupancy**.
- `calc_Gv4_in`: 60 regs → 2,048 regs/warp → 8 blocks → 50% occupancy.
  At ≤48 regs → 1,536 regs/warp → **10 blocks → 62.5% occupancy** (shared mem caps at 10 anyway).

---

## Root Cause: What is consuming the registers

### `calc_step2_3` (`3D_solver/src/calc_steps.f90`, line 111)

The inlined `calc_R` call fills `R(5)` (a 5-element double array) by doing all-x then all-y then all-z accumulations across all components simultaneously. This forces **10 registers** for `R(1:5)` to remain live until the update loop consumes them.

```fortran
! Current pattern (10 regs for R alive from calc_R return until loop end):
real(8) R(5), ...
call calc_R(..., E, F, G, R)       ! fills all R(1:5) at once
do l = 1, 5
  Qout(...) = (... - R(l)) * coef4_inv
enddo
```

### `calc_Gv4_in` (`3D_solver/src/calc_visc4_internal.f90`, line 132)

Seven stress/flux scalars (`tzx, tzy, tzz, utzx, vtzy, wtzz, kTz`) are declared at subroutine scope and assigned sequentially inside a `block`, but all are used together only in the final G updates **outside** the block. This forces **14 registers** for all 7 scalars to be live simultaneously from the first assignment through the final G write.

```fortran
! Current pattern (all 7 tau outputs live simultaneously = 14 regs):
real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz
...
block
  kTz = flux4(...)          ! kTz assigned
  call calc_tau_straight(..., tzz, wtzz)   ! tzz, wtzz assigned
  call calc_tau_cross(..., tzx, utzx)      ! tzx, utzx assigned
  call calc_tau_cross(..., tzy, vtzy)      ! tzy, vtzy assigned
end block
G(2,...) -= tzx
G(3,...) -= tzy
G(4,...) -= tzz
G(5,...) -= (utzx + vtzy + wtzz + kTz)   ! all 7 needed here
```

---

## Changes

### Change 1 — `calc_step2_3`: replace `R(5)` array with scalar `Rl` in loop

**File:** `3D_solver/src/calc_steps.f90` — lines 127–144

Replace `real(8) R(5)` with `real(8) Rl`, remove the `call calc_R(...)`, and inline the flux divergence directly as a scalar expression inside the loop:

```fortran
! Before:
real(8) R(5), coef3_dtdxdy, coef3_dtdydz, coef3_dtdzdx
integer i, j, k, l
...
coef3_dtdxdy = coef3 * dtdxdy(i,j)
coef3_dtdydz = coef3 * dtdydz(j,k)
coef3_dtdzdx = coef3 * dtdzdx(i,k)
call calc_R(nx, ny, nz, i, j, k, coef3_dtdxdy, coef3_dtdydz, coef3_dtdzdx, E, F, G, R)
do l = 1, 5
  Qout(i+1,l,j+1,k+1) = (coef1 * Qin(i+1,l,j+1,k+1) + coef2 * Qout(i+1,l,j+1,k+1) - R(l)) * coef4_inv
enddo

! After:
real(8) Rl, coef3_dtdxdy, coef3_dtdydz, coef3_dtdzdx
integer i, j, k, l
...
coef3_dtdxdy = coef3 * dtdxdy(i,j)
coef3_dtdydz = coef3 * dtdydz(j,k)
coef3_dtdzdx = coef3 * dtdzdx(i,k)
do l = 1, 5
  Rl = coef3_dtdydz * (E(l,i+1,j,k) - E(l,i,j,k)) &
     + coef3_dtdzdx * (F(l,i,j+1,k) - F(l,i,j,k)) &
     + coef3_dtdxdy * (G(l,i,j,k+1) - G(l,i,j,k))
  Qout(i+1,l,j+1,k+1) = (coef1 * Qin(i+1,l,j+1,k+1) + coef2 * Qout(i+1,l,j+1,k+1) - Rl) * coef4_inv
enddo
```

**Register saving:** `R(5)` = 10 regs → `Rl` = 2 regs → **saves ~8 registers** (42 → ~34).
This is well below the ≤40 threshold. The `calc_R` device subroutine is kept unchanged (still used by `calc_step1`, `calc_step`, `calc_step4`).

Also **remove `maxregcount=40`** from the subroutine attribute (it was added in the previous session and is no longer needed).

**Sign note:** `calc_R` uses `(-E(l,i,...) + E(l,i+1,...))` for the x-direction, which equals `(E(l,i+1,...) - E(l,i,...))`. The inline version must match this convention exactly.

---

### Change 2 — `calc_Ev4_in`, `calc_Fv4_in`, `calc_Gv4_in`: move tau outputs into nested `block` scopes

**File:** `3D_solver/src/calc_visc4_internal.f90`

All three kernels share the same anti-pattern: 7 stress/flux scalars declared at subroutine scope are set inside a `block` but all consumed together in one expression outside the block, keeping 14 registers live simultaneously. The fix is identical for each: move each pair of outputs into its own `block`, apply the flux array update immediately, then let the block end free those registers.

Also remove unused loop-index variables (`ii` in `calc_Ev4_in`, `jj` in `calc_Fv4_in`, `kk` in `calc_Gv4_in` — all declared but never assigned in the kernel body).
Also **remove `maxregcount=48`** from `calc_Gv4_in`'s attribute (no longer needed).

#### `calc_Ev4_in` (lines 40–69)

```fortran
! Before:
integer i, j, k, it, jt, kt, ii, idx, offset_yz
real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx
...
if (3 <= i ... ) then
  block
    real(8) :: mu3(3)
    mu3(:) = ...
    block ! dQdx
      real(8) :: kTx3(3)
      kTx3(:) = ... * dx(i)
      kTx     = flux4(kTx3(:))
    end block
    call calc_tau_straight(mu3, u(idx-2:idx+3), vy(idx-2:idx+3), wz(idx-2:idx+3), dx(i), txx, utxx)
    call calc_tau_cross(mu3, v(idx-2:idx+3), uy(idx-2:idx+3), dx(i), txy, vtxy)
    call calc_tau_cross(mu3, w(idx-2:idx+3), uz(idx-2:idx+3), dx(i), txz, wtxz)
  end block
  E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
  E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
  E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
  E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx)
endif

! After:
integer i, j, k, it, jt, kt, idx, offset_yz      ! ii removed, no tau vars at subroutine scope
...
if (3 <= i ... ) then
  block
    real(8) :: mu3(3)
    mu3(:) = ...
    block
      real(8) :: kTx3(3), kTx
      kTx3(:) = ... * dx(i)
      kTx     = flux4(kTx3(:))
      E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - kTx
    end block
    block
      real(8) :: txx, utxx
      call calc_tau_straight(mu3, u(idx-2:idx+3), vy(idx-2:idx+3), wz(idx-2:idx+3), dx(i), txx, utxx)
      E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
      E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - utxx
    end block
    block
      real(8) :: txy, vtxy
      call calc_tau_cross(mu3, v(idx-2:idx+3), uy(idx-2:idx+3), dx(i), txy, vtxy)
      E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
      E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - vtxy
    end block
    block
      real(8) :: txz, wtxz
      call calc_tau_cross(mu3, w(idx-2:idx+3), uz(idx-2:idx+3), dx(i), txz, wtxz)
      E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
      E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - wtxz
    end block
  end block
endif
```

#### `calc_Fv4_in` (lines 97–126) — identical pattern for F flux

```fortran
! Before:
integer i, j, k, it, jt, kt, jj, idx, offset_xz
real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy
...
  block
    ...
    block ! dQdy
      real(8) :: kTy3(3)
      kTy3(:) = ... * dy(j)
      kTy     = flux4(kTy3(:))
    end block
    call calc_tau_straight(mu3, v(idx-2:idx+3), wz(idx-2:idx+3), ux(idx-2:idx+3), dy(j), tyy, vtyy)
    call calc_tau_cross(mu3, u(idx-2:idx+3), vx(idx-2:idx+3), dy(j), tyx, utyx)
    call calc_tau_cross(mu3, w(idx-2:idx+3), vz(idx-2:idx+3), dy(j), tyz, wtyz)
  end block
  F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
  F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
  F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
  F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy)

! After:
integer i, j, k, it, jt, kt, idx, offset_xz      ! jj removed, no tau vars at subroutine scope
...
  block
    ...
    block
      real(8) :: kTy3(3), kTy
      kTy3(:) = ... * dy(j)
      kTy     = flux4(kTy3(:))
      F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - kTy
    end block
    block
      real(8) :: tyy, vtyy
      call calc_tau_straight(mu3, v(idx-2:idx+3), wz(idx-2:idx+3), ux(idx-2:idx+3), dy(j), tyy, vtyy)
      F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
      F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - vtyy
    end block
    block
      real(8) :: tyx, utyx
      call calc_tau_cross(mu3, u(idx-2:idx+3), vx(idx-2:idx+3), dy(j), tyx, utyx)
      F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
      F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - utyx
    end block
    block
      real(8) :: tyz, wtyz
      call calc_tau_cross(mu3, w(idx-2:idx+3), vz(idx-2:idx+3), dy(j), tyz, wtyz)
      F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
      F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - wtyz
    end block
  end block
```

#### `calc_Gv4_in` (lines 154–183) — same pattern for G flux

```fortran
! Before:
integer i, j, k, it, jt, kt, kk, idx, offset_xy
real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz
...
  block
    ...
    block ! dQdz
      real(8) :: kTz3(3)
      kTz3(:) = ... * dz(k)
      kTz     = flux4(kTz3(:))
    end block
    call calc_tau_straight(mu3, w(idx-2:idx+3), ux(idx-2:idx+3), vy(idx-2:idx+3), dz(k), tzz, wtzz)
    call calc_tau_cross(mu3, u(idx-2:idx+3), wx(idx-2:idx+3), dz(k), tzx, utzx)
    call calc_tau_cross(mu3, v(idx-2:idx+3), wy(idx-2:idx+3), dz(k), tzy, vtzy)
  end block
  G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
  G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
  G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
  G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz)

! After:
integer i, j, k, it, jt, kt, idx, offset_xy      ! kk removed, no tau vars at subroutine scope
...
  block
    ...
    block
      real(8) :: kTz3(3), kTz
      kTz3(:) = ... * dz(k)
      kTz     = flux4(kTz3(:))
      G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - kTz
    end block
    block
      real(8) :: tzz, wtzz
      call calc_tau_straight(mu3, w(idx-2:idx+3), ux(idx-2:idx+3), vy(idx-2:idx+3), dz(k), tzz, wtzz)
      G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
      G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - wtzz
    end block
    block
      real(8) :: tzx, utzx
      call calc_tau_cross(mu3, u(idx-2:idx+3), wx(idx-2:idx+3), dz(k), tzx, utzx)
      G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
      G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - utzx
    end block
    block
      real(8) :: tzy, vtzy
      call calc_tau_cross(mu3, v(idx-2:idx+3), wy(idx-2:idx+3), dz(k), tzy, vtzy)
      G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
      G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - vtzy
    end block
  end block
```

**Register saving (all three kernels):** 7 tau outputs simultaneously live (14 regs) → max 2 outputs at a time (4 regs) + mu3(3) in outer block (6 regs). Estimated **saves ~10–12 registers**. The direction-specific spacing scalar (`dx(i)`, `dy(j)`, `dz(k)`) is re-loaded inside each inner block from L1 rather than held persistently, saving an additional 2 registers.

**Trade-off:** The energy flux component (E/F/G index 5) is written 4 separate times instead of 1 combined. Each is an L1-cached RMW to the same address by the same thread — 3 extra cache hits, negligible cost.

---

## Expected Impact

| Kernel | Metric | Current | After |
|---|---|---|---|
| `calc_step2_3` | Registers/thread | 42 | ~34 |
| `calc_step2_3` | Block limit (regs) | 10 | 12 |
| `calc_step2_3` | Theoretical occupancy | 62.5% | 75% |
| `calc_Ev4_in`, `calc_Fv4_in` | Registers/thread | similar to Gv | ~10–12 lower |
| `calc_Gv4_in` | Registers/thread | 60 | ~48–50 |
| `calc_Gv4_in` | Block limit (regs) | 8 | 10 (shm-capped) |
| `calc_Gv4_in` | Theoretical occupancy | 50% | 62.5% |

---

## Risk: What if register estimates are wrong?

The register savings are estimated from first principles; the actual count depends on nvfortran's register allocator. After compiling, run `ptxinfo` output (already in the Makefile via `-gpu=ptxinfo`) to check actual register counts.

- `calc_step2_3`: needs ≤40. We estimate ~34. Margin = 6 registers. Low risk.
- `calc_Gv4_in`: needs ≤48. We estimate ~48–50. **If the result is 50–52, add `maxregcount=48` back as a minimal nudge** — the code restructuring still reduces spilling cost vs. applying maxregcount=48 to the original 60-register version.

---

## Verification

1. `cd 3D_solver/NSTGV && make clean && make` — check ptxinfo output for register counts
2. `bash calc.sh` — run for a few time steps, confirm no divergence
3. Re-profile with `bash profile.sh` / `ncu` and check:
   - `calc_step2_3`: Block Limit Registers → 12, Achieved Occupancy ~73–75%
   - `calc_Gv4_in`: Block Limit Registers → 10, Achieved Occupancy ~60–63%
