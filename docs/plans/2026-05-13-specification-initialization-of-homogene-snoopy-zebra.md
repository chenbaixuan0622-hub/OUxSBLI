# Plan: Mean-Temperature Relaxation Cooling for Weakly Compressible HIT

## Context

The DHIT case now has Petersen-type solenoidal forcing (implemented; `calc_forcing.f90`, `#ifdef USE_FORCING` guards in `calc_time_dev.f90`). Continuous energy injection causes long-time temperature drift. This plan adds a simple thermal relaxation term to the energy equation only:

```
S_E = -ρ · C_T · (T - T_ref)
```

`T_ref` is the initial uniform temperature — already the compile-time parameter `T` in `mod_globals.f90` (`T = urms²/(γ·R·Mt²)`). The cooling is applied at every RK substep, co-located with the forcing.

---

## Scope — Files to Modify

| File | Change |
|------|--------|
| `3D_solver/DHIT/mod_globals.f90` | Add `C_T` parameter |
| `3D_solver/DHIT/calc_forcing.f90` | Extend `use` imports; add `add_cooling_k` kernel; add `apply_cooling` subroutine |
| `3D_solver/src/calc_time_dev.f90` | Add `call apply_cooling(...)` after each `call add_forcing` in `#ifdef USE_FORCING` blocks |

No new files. Stub (`3D_solver/src/calc_forcing.f90`) and Makefile are unchanged — `apply_cooling` is only called inside `#ifdef USE_FORCING` blocks, so non-DHIT cases never see it.

---

## 1. `mod_globals.f90` — New Parameter

After the `kf_max` line (~line 117):

```fortran
real(8), parameter :: C_T    = 0.02d0  ! temperature relaxation coefficient
```

---

## 2. `calc_forcing.f90` — Cooling Kernel + Wrapper

### 2a. Update `use mod_globals` (top of module)

Current:
```fortran
use mod_globals, only : id_accuracy, nx, ny, nz, eps_s, kf_min, kf_max, threads
```
Change to:
```fortran
use mod_globals, only : id_accuracy, nx, ny, nz, eps_s, kf_min, kf_max, threads, &
                        gamma, R, C_T, dt, T_ref_const => T
```
(`T` is renamed `T_ref_const` to distinguish it from any local temperature variable.)

### 2b. New GPU kernel `add_cooling_k` (add after `store_forcing_k`, before `calc_forcing_rhs`)

```fortran
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
```

### 2c. New subroutine `apply_cooling` (add after `calc_forcing_rhs`, before `finalize_forcing`)

```fortran
  subroutine apply_cooling(nx_, ny_, nz_, coef, QJ)
    integer, intent(in) :: nx_, ny_, nz_
    real(8), intent(in) :: coef
    real(8), device     :: QJ(nx_,5,ny_,nz_)
    type(dim3) :: blk
    blk = dim3((nx_-2+31)/32, (ny_-2+3)/4, nz_-2)
    call add_cooling_k<<<blk, threads>>>(nx_, ny_, nz_, coef, QJ)
  end subroutine apply_cooling
```

Block grid covers `(nx-2)×(ny-2)×(nz-2)` interior cells with `threads=dim3(32,4,1)`. For 65³: blk=(2,16,63), threads=(32,4,1).

---

## 3. `calc_time_dev.f90` — Add Cooling at Each RK Substep

Inside the three `#ifdef USE_FORCING` blocks in `RungeKutta_3rd`, add `call apply_cooling(...)` immediately after each `call add_forcing`:

**Stage 1** (currently lines 99–101):
```fortran
#ifdef USE_FORCING
          call add_forcing<<<blocks,threads>>>(nx, ny, nz, 1.d0, QJ, QJ2, fx_d, fy_d, fz_d)
          call apply_cooling(nx, ny, nz, 1.d0, QJ2)
#endif
```

**Stage 2** (currently lines 111–113):
```fortran
#ifdef USE_FORCING
          call add_forcing<<<blocks,threads>>>(nx, ny, nz, 0.25d0, QJ2, QJ2, fx_d, fy_d, fz_d)
          call apply_cooling(nx, ny, nz, 0.25d0, QJ2)
#endif
```

**Stage 3** (currently lines 122–124):
```fortran
#ifdef USE_FORCING
          call add_forcing<<<blocks,threads>>>(nx, ny, nz, 2.d0*one_third, QJ, QJ, fx_d, fy_d, fz_d)
          call apply_cooling(nx, ny, nz, 2.d0*one_third, QJ)
#endif
```

RK coefficients (1.0, 0.25, 2/3) match those of the forcing — both are source terms `dQ/dt = f`.

---

## Physics Notes

- `T_ref_const` = `T` from mod_globals = `urms²/(γ·R·Mt²)` ≈ 2489 K (initial uniform temperature)
- Cooling time scale: τ_T = 1/C_T = 50 [1/C_T units]. Spec requires C_T·τ_eddy ≪ 1; with τ_eddy ~ L0/urms = 2π/100 ≈ 0.063 s and C_T = 0.02, C_T·τ_eddy ≈ 0.0013 ≪ 1 ✓
- Cooling modifies only QJ(:,5,:,:) — momentum QJ(:,2:4,:,:) is untouched, preserving turbulence structure
- T_loc is computed locally from conservative variables: e_int = E - ½|u|², T = e_int·(γ-1)/R

---

## Verification

1. **Compile:**
   ```bash
   cd 3D_solver/DHIT && make clean && make
   ```
   Expected: no new errors; `calc_forcing.mod` contains `apply_cooling`.

2. **Run:**
   ```bash
   bash calc.sh
   ```
   Monitor stdout. With forcing+cooling, temperature should stabilize near T_ref rather than drifting.

3. **Check stationarity:**
   ```bash
   python3 check_hit.py
   ```
   Expected: TKE stationary; mean temperature stays near T_ref.

---

## Critical Files

- [`3D_solver/DHIT/mod_globals.f90`](3D_solver/DHIT/mod_globals.f90) — add `C_T` after `kf_max` (line ~117)
- [`3D_solver/DHIT/calc_forcing.f90`](3D_solver/DHIT/calc_forcing.f90) — extend use statement; add `add_cooling_k` and `apply_cooling`
- [`3D_solver/src/calc_time_dev.f90`](3D_solver/src/calc_time_dev.f90) — add `call apply_cooling(...)` at lines ~100, ~112, ~123
