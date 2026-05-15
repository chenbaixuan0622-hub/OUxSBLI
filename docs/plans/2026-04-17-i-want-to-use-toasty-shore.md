# Plan: Flatten smem arrays in visc4 files (TMA Step 1)

## Context

Preparation for Tensor Memory Accelerator (TMA) usage requires shared memory arrays to be 1D contiguous buffers. Currently all smem arrays in the visc4 kernels are declared as 3D. This plan collapses them to 1D following the same convention used in `calc_keep_kernel_internal.f90`.

## KEEP convention (reference: `calc_keep_kernel_internal.f90`)

```fortran
integer, parameter :: io = kind(id_accuracy) / 3   ! stencil order param: 0,1,2 for 2nd,4th,6th
integer, parameter :: sx = threadsE%x + 2*io + 1   ! padded tile size (includes halo)
integer, parameter :: sy = threadsE%y
integer, parameter :: sz = threadsE%z
real(8), dimension(-(io-1):sx*sy*sz-io), shared :: rho, u, ...
offset_yz = (jt-1)*sx + (kt-1)*sx*sy   ! linearizes perpendicular dims
idx       = ii + offset_yz              ! smem index (load)
idx       = it + offset_yz              ! smem index (compute)
rho(idx)  = Q(...)                      ! load: direct index
rho(idx-io:idx+io+1)                    ! consume: stencil slice
```
Key: `idx` IS the smem address directly — no `+3` offset.

## Visc4 uses the same stencil pattern as KEEP io=2 (6th order)

The visc4 stencil `u(it-2:it+3, ...)` = `u(idx-2:idx+3)` matches KEEP's `rho(idx-io:idx+io+1)` with `io=2`. All load loops `do ii = it-2, threads%x+3` also match `it-io` to `threads%x+io+1` with `io=2`.

Define `integer, parameter :: io_v = 2` in all affected subroutines/kernels.

## Thread block constants (from `mod_globals.f90`)
```
threadsEv = dim3(32, 1, 1)   threadsFv = dim3(32, 4, 1)   threadsGv = dim3(32, 1, 4)
```

## Tile size parameters per direction

| Direction | sx | sy | sz | Total |
|-----------|----|----|-----|-------|
| Ev (x-kernel) | `threadsEv%x + 2*io_v + 1` = 37 | `threadsEv%y` = 1 | `threadsEv%z` = 1 | 37 |
| Fv (y-kernel) | `threadsFv%x` = 32 | `threadsFv%y + 2*io_v + 1` = 9 | `threadsFv%z` = 1 | 288 |
| Gv (z-kernel) | `threadsGv%x` = 32 | `threadsGv%y` = 1 | `threadsGv%z + 2*io_v + 1` = 9 | 288 |

## Offset and index formulas

| Direction | offset | idx (load) | idx (compute) |
|-----------|--------|------------|---------------|
| Ev | `offset_yz = (jt-1)*sx + (kt-1)*sx*sy` = 0 | `ii + offset_yz` | `it + offset_yz` |
| Fv | `offset_xz = (it-1)*sy + (kt-1)*sy*sx` = `(it-1)*9` | `jj + offset_xz` | `jt + offset_xz` |
| Gv | `offset_xy = (jt-1)*sz + (it-1)*sz*sy` = `(it-1)*9` | `kk + offset_xy` | `kt + offset_xy` |

## Smem declaration (1D, identical for all arrays in a kernel)
```fortran
real(8), dimension(-(io_v-1):sx*sy*sz-io_v), shared :: u, v, w, ...
! Ev: dimension(-1:35), Fv: dimension(-1:286), Gv: dimension(-1:286)
```

## Index conversion table (consuming code)

| Old multi-dim access | New 1D access |
|---------------------|---------------|
| `A(it+δ, jt, kt)` | `A(idx+δ)` |
| `A(it-2:it+3, jt, kt)` | `A(idx-2:idx+3)` |
| `A(it-1:it+2, jt, kt)` | `A(idx-1:idx+2)` |
| `A(it:it+1, jt, kt)` | `A(idx:idx+1)` |
| `A(jt+δ, it, kt)` | `A(idx+δ)` |
| `A(jt-2:jt+3, it, kt)` | `A(idx-2:idx+3)` |
| `A(jt-1:jt+2, it, kt)` | `A(idx-1:idx+2)` |
| `A(jt:jt+1, it, kt)` | `A(idx:idx+1)` |
| `A(kt+δ, jt, it)` | `A(idx+δ)` |
| `A(kt-2:kt+3, jt, it)` | `A(idx-2:idx+3)` |
| `A(kt-1:kt+2, jt, it)` | `A(idx-1:idx+2)` |
| `A(kt:kt+1, jt, it)` | `A(idx:idx+1)` |

---

## Files and changes

### 1. `3D_solver/src/load_smem_visc4.f90`

Each subroutine gets local parameters `io_v`, `sx`, `sy`, `sz` and offset variable.
`uz`/`wz` in the x-subroutine had swapped 2nd/3rd dims in old code — after 1D flattening they use the same declaration as all other arrays.

#### `load_smem_visc4_x` (lines 8-46)
Add local parameters (after implicit none in the subroutine):
```fortran
integer, parameter :: io_v = 2
integer, parameter :: sx = threadsEv%x + 2*io_v + 1
integer, parameter :: sy = threadsEv%y
integer, parameter :: sz = threadsEv%z
```
Dummy arg declarations (lines 21-27) — change all 7 from `(-2:threadsEv%x+3, threadsEv%y, threadsEv%z)` to `(-(io_v-1):sx*sy*sz-io_v)`.
Add `integer :: idx, offset_yz` to locals (line 28).
Add before the `do ii` loop (line 30): `offset_yz = (jt-1)*sx + (kt-1)*sx*sy`
In the loop body: `idx = ii + offset_yz`, then replace all `X(ii, jt, kt)` with `X(idx)` (lines 33-44).

#### `load_smem_visc4_y` (lines 49-87)
```fortran
integer, parameter :: io_v = 2
integer, parameter :: sx = threadsFv%x
integer, parameter :: sy = threadsFv%y + 2*io_v + 1
integer, parameter :: sz = threadsFv%z
```
Dummy args (lines 62-68): change to `(-(io_v-1):sx*sy*sz-io_v)`.
Add `integer :: idx, offset_xz`, set `offset_xz = (it-1)*sy + (kt-1)*sy*sx` before loop.
In loop: `idx = jj + offset_xz`, replace `X(jj, it, kt)` → `X(idx)`.

#### `load_smem_visc4_z` (lines 90-128)
```fortran
integer, parameter :: io_v = 2
integer, parameter :: sx = threadsGv%x
integer, parameter :: sy = threadsGv%y
integer, parameter :: sz = threadsGv%z + 2*io_v + 1
```
Dummy args (lines 103-109): change to `(-(io_v-1):sx*sy*sz-io_v)`.
Add `integer :: idx, offset_xy`, set `offset_xy = (jt-1)*sz + (it-1)*sz*sy` before loop.
In loop: `idx = kk + offset_xy`, replace `X(kk, jt, it)` → `X(idx)`.

---

### 2. `3D_solver/src/calc_visc4.f90`

Each kernel follows the same KEEP pattern. Add `io_v`, `sx`, `sy`, `sz` as local parameters, add `integer :: idx, offset_*`, set offset after thread index assignments, set `idx` after syncthreads.

#### `calc_Ev4` and `calc_Ev_LES4`
Local params:
```fortran
integer, parameter :: io_v = 2
integer, parameter :: sx = threadsEv%x + 2*io_v + 1
integer, parameter :: sy = threadsEv%y
integer, parameter :: sz = threadsEv%z
```
Smem declarations (7 arrays): `(-2:threadsEv%x+3, threadsEv%y, threadsEv%z)` → `(-(io_v-1):sx*sy*sz-io_v)`. Eliminate the old swapped-dim variant for `uz`/`wz` (all 7 use same declaration).
Add `integer :: idx, offset_yz` to locals.
After thread index lines: `offset_yz = (jt-1)*sx + (kt-1)*sx*sy`
After syncthreads: `idx = it + offset_yz`
Replace all smem accesses using the table above.

#### `calc_Fv4` and `calc_Fv_LES4`
Local params:
```fortran
integer, parameter :: io_v = 2
integer, parameter :: sx = threadsFv%x
integer, parameter :: sy = threadsFv%y + 2*io_v + 1
integer, parameter :: sz = threadsFv%z
```
Smem declarations (7 arrays): `(-2:threadsFv%y+3, threadsFv%x, threadsFv%z)` → `(-(io_v-1):sx*sy*sz-io_v)`.
Add `integer :: idx, offset_xz`.
After thread index lines: `offset_xz = (it-1)*sy + (kt-1)*sy*sx`
After syncthreads: `idx = jt + offset_xz`
Replace all smem accesses.

#### `calc_Gv4` and `calc_Gv_LES4`
Local params:
```fortran
integer, parameter :: io_v = 2
integer, parameter :: sx = threadsGv%x
integer, parameter :: sy = threadsGv%y
integer, parameter :: sz = threadsGv%z + 2*io_v + 1
```
Smem declarations (7 arrays): `(-2:threadsGv%z+3, threadsGv%y, threadsGv%x)` → `(-(io_v-1):sx*sy*sz-io_v)`.
Add `integer :: idx, offset_xy`.
After thread index lines: `offset_xy = (jt-1)*sz + (it-1)*sz*sy`
After syncthreads: `idx = kt + offset_xy`
Replace all smem accesses.

---

### 3. `3D_solver/src/calc_visc4_internal.f90`

Apply identical changes to `calc_Ev4_in`, `calc_Fv4_in`, `calc_Gv4_in`. These have no boundary path — only 4th-order `calc_tau_*` call sites to update.

### 4. `3D_solver/src/calc_visc4_les_internal.f90`

Apply identical changes to `calc_Ev_LES4_in`, `calc_Fv_LES4_in`, `calc_Gv_LES4_in`. These have the LES Hsgs 4-elem slices in addition to the `calc_tau_*` calls.

---

## No changes needed
- `calc_visc_me4_base.f90` (`calc_tau_straight`, `calc_tau_cross`, `flux4`): receive assumed-shape 1D slices — unchanged since slice length (6 elements) is unchanged.

## Verification
```bash
cd 3D_solver/NSTGV
make clean && make   # must compile without errors
bash calc.sh         # run and compare VTK output with baseline
```
