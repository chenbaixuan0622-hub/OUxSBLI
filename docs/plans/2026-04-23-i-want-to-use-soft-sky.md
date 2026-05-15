# Plan: TMA in `calc_slau_kernel*` — pipeline Q load with MUSCL interpolation

## Context

The SLAU x-kernels currently load all 5 Q components into shared memory via a scalar loop + `syncthreads()`, then begin MUSCL. This serializes memory transfer and compute.

The user wants to pipeline: issue TMA load for variable `n+1` while computing MUSCL on variable `n`. Five separate TMA barriers — one per Q component — allow each to be consumed as soon as it is ready.

TMA is already used in `load_smem_visc4_x` ([3D_solver/src/load_smem_visc4.f90:8-79](3D_solver/src/load_smem_visc4.f90)); that pattern (single barrier) is extended here to 5 barriers.

**TMA applicability:** x-direction only. `Q(i,comp,j,k)` for consecutive `i` is contiguous. y/z have stride `nx` or `nx*5*ny` → TMA bulk load does not apply.

## Files to Modify

| File | Kernels |
|------|---------|
| [3D_solver/src/calc_slau_kernel.f90](3D_solver/src/calc_slau_kernel.f90) | `calc_slau_x6`, `calc_slau_x4` |
| [3D_solver/src/calc_slau_kernel_internal.f90](3D_solver/src/calc_slau_kernel_internal.f90) | `calc_slau_x_in` |

`calc_slau_x2` is excluded — it hardcodes `fdx = 1.d0` and does no MUSCL.

## Execution Timeline

```
Thread 0:  tma_bulk_load(bar_rho, Q(:,1,...), rho)  ─┐
           tma_bulk_load(bar_u,   Q(:,2,...),   u)   │ all 5 TMA
           tma_bulk_load(bar_v,   Q(:,3,...),   v)   │ transfers
           tma_bulk_load(bar_w,   Q(:,4,...),   w)   │ in flight
           tma_bulk_load(bar_p,   Q(:,5,...),   p)  ─┘

All threads (overlapped with TMA):
           read sensor → fdx (register)

           arrive(bar_rho) → wait → delta6 rho → rhol, rhor(idx_r)
           arrive(bar_u)   → wait → delta6 u   → ul, ur(idx_r)    ← bar_v,w,p still completing
           arrive(bar_v)   → wait → delta6 v   → vl, vr(idx_r)
           arrive(bar_w)   → wait → delta6 w   → wl, wr(idx_r)
           arrive(bar_p)   → wait → delta6 p   → pl, pr(idx_r) + wiggle_detector
           SLAU flux (unchanged)
```

No `syncthreads()` is needed between variable stages because:
- Each stage reads from `rho/u/v/w/p` (separate smem arrays, read-only after TMA)
- Each stage writes to `rhor/ur/vr/wr/pr` at unique index `idx_r` (per-thread, no conflicts)

## TMA Index Calculation (per kernel)

For `threadsE = dim3(32,1,1)`: `offset_yz = 0` always (jt=1, kt=1).

| Kernel | `start_i` | `end_i` | `start_idx` | smem lower bound |
|--------|-----------|---------|-------------|-----------------|
| `calc_slau_x6` | `max(1, i_base - 1)` | `min(nx, i_base + threadsE%x + 3)` | `start_i - i_base` | `-1` |
| `calc_slau_x4` | `max(1, i_base)` | `min(nx, i_base + threadsE%x + 2)` | `start_i - i_base` | `0` |
| `calc_slau_x_in` | `max(1, i_base + 1 - io)` | `min(nx, i_base + threadsE%x + io + 1)` | `start_i - i_base` | `-(io-1)` |

## Variable Declarations to Add (inside `#if _CUDA_ARCH_ >= 900` block only)

```fortran
integer(8), shared :: bar_rho, bar_u, bar_v, bar_w, bar_p
integer(8)         :: token
integer            :: start_i, end_i, count, start_idx, tid_linear
```

## Code Structure (x6 shown; x4 is identical with `io`→1, stencil offsets adjusted)

Replace the existing load loop + syncthreads + early-return + fdx read with:

```fortran
#if _CUDA_ARCH_ >= 900
    tid_linear = threadIdx%x + (threadIdx%y-1)*blockDim%x &
                + (threadIdx%z-1)*blockDim%x*blockDim%y
    if (tid_linear == 1) then
      call barrier_init(bar_rho, blockDim%x * blockDim%y * blockDim%z)
      call barrier_init(bar_u,   blockDim%x * blockDim%y * blockDim%z)
      call barrier_init(bar_v,   blockDim%x * blockDim%y * blockDim%z)
      call barrier_init(bar_w,   blockDim%x * blockDim%y * blockDim%z)
      call barrier_init(bar_p,   blockDim%x * blockDim%y * blockDim%z)
    endif
    call syncthreads()
    start_i   = max(1,  i_base - 1)
    end_i     = min(nx, i_base + threadsE%x + 3)
    count     = max(0, end_i - start_i + 1)
    start_idx = start_i - i_base
    if (tid_linear == 1) then
      if (count > 0 .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        call tma_bulk_load(bar_rho, Q(start_i,1,j,k), rho(start_idx), count)
        call tma_bulk_load(bar_u,   Q(start_i,2,j,k),   u(start_idx), count)
        call tma_bulk_load(bar_v,   Q(start_i,3,j,k),   v(start_idx), count)
        call tma_bulk_load(bar_w,   Q(start_i,4,j,k),   w(start_idx), count)
        call tma_bulk_load(bar_p,   Q(start_i,5,j,k),   p(start_idx), count)
      endif
    endif
    ! Overlapped with TMA: compute i, idx, fdx
    i   = (blockIdx%x-1)*blockDim%x + it
    idx   = it + offset_yz
    idx_r = it + offset_yzr
    fdx   = 0.5d0 * (sensor(i,j,k) + sensor(i+1,j,k))
    call syncthreads()
    ! Threads outside compute range must still participate in all 5 barriers
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) then
      token = barrier_arrive(bar_rho)
      do; if (barrier_try_wait_sleep(bar_rho, token, 1000000) .ne. 0) exit; enddo
      token = barrier_arrive(bar_u)
      do; if (barrier_try_wait_sleep(bar_u, token, 1000000) .ne. 0) exit; enddo
      token = barrier_arrive(bar_v)
      do; if (barrier_try_wait_sleep(bar_v, token, 1000000) .ne. 0) exit; enddo
      token = barrier_arrive(bar_w)
      do; if (barrier_try_wait_sleep(bar_w, token, 1000000) .ne. 0) exit; enddo
      token = barrier_arrive(bar_p)
      do; if (barrier_try_wait_sleep(bar_p, token, 1000000) .ne. 0) exit; enddo
      return
    endif
    ! --- rho ---
    token = barrier_arrive(bar_rho)
    do; if (barrier_try_wait_sleep(bar_rho, token, 1000000) .ne. 0) exit; enddo
    if (3 <= i .and. i <= nx-3) then
      call delta6(fdx, rho(idx-2:idx+3), rhol, rhor(idx_r))
    elseif (2 <= i .and. i <= nx-2) then
      call delta4(fdx, rho(idx-1:idx+2), rhol, rhor(idx_r))
    else
      rhol = rho(idx); rhor(idx_r) = rho(idx+1)
    endif
    ! --- u ---
    token = barrier_arrive(bar_u)
    do; if (barrier_try_wait_sleep(bar_u, token, 1000000) .ne. 0) exit; enddo
    if (3 <= i .and. i <= nx-3) then
      call delta6(fdx, u(idx-2:idx+3), ul, ur(idx_r))
    elseif (2 <= i .and. i <= nx-2) then
      call delta4(fdx, u(idx-1:idx+2), ul, ur(idx_r))
    else
      ul = u(idx); ur(idx_r) = u(idx+1)
    endif
    ! --- v ---
    token = barrier_arrive(bar_v)
    do; if (barrier_try_wait_sleep(bar_v, token, 1000000) .ne. 0) exit; enddo
    if (3 <= i .and. i <= nx-3) then
      call delta6(fdx, v(idx-2:idx+3), vl, vr(idx_r))
    elseif (2 <= i .and. i <= nx-2) then
      call delta4(fdx, v(idx-1:idx+2), vl, vr(idx_r))
    else
      vl = v(idx); vr(idx_r) = v(idx+1)
    endif
    ! --- w ---
    token = barrier_arrive(bar_w)
    do; if (barrier_try_wait_sleep(bar_w, token, 1000000) .ne. 0) exit; enddo
    if (3 <= i .and. i <= nx-3) then
      call delta6(fdx, w(idx-2:idx+3), wl, wr(idx_r))
    elseif (2 <= i .and. i <= nx-2) then
      call delta4(fdx, w(idx-1:idx+2), wl, wr(idx_r))
    else
      wl = w(idx); wr(idx_r) = w(idx+1)
    endif
    ! --- p (also updates fdx for wiggle detector) ---
    token = barrier_arrive(bar_p)
    do; if (barrier_try_wait_sleep(bar_p, token, 1000000) .ne. 0) exit; enddo
    if (3 <= i .and. i <= nx-3) then
      call delta6(fdx, p(idx-2:idx+3), pl, pr(idx_r))
      fdx = wiggle_detector(p(idx-1:idx+2))
    elseif (2 <= i .and. i <= nx-2) then
      call delta4(fdx, p(idx-1:idx+2), pl, pr(idx_r))
      fdx = wiggle_detector(p(idx-1:idx+2))
    else
      pl = p(idx); pr(idx_r) = p(idx+1); fdx = 1.d0
    endif
    ! SLAU flux (unchanged)
    associate(un1 => ul, un2 => ur(idx_r))
      call SLAU(...)
    end associate
#else
    ! original loop + syncthreads unchanged (verbatim)
    do ii = it-2, threadsE%x+3, blockDim%x
      ...
    enddo
    call syncthreads()
    i = ...; if (...) return
    idx = ...; idx_r = ...; fdx = ...
    if (3 <= i ...) then
      call delta6(fdx, rho(...), ...); ...
    endif
    ...
    call SLAU(...)
#endif
```

## `calc_slau_x_in` — per-variable inlined deltas

In `x_in`, `io = kind(id_accuracy)/3` is a compile-time constant. Use `select case (io)` (dead branches eliminated at compile time) for each variable:

```fortran
! After barrier_arrive(bar_rho) + wait:
select case (io)
case (0); rhol = rho(i1); rhor(idx_r) = rho(i2)
case (1); call delta4(fdx, rho(i1:i2), rhol, rhor(idx_r))
case (2); call delta6(fdx, rho(i1:i2), rhol, rhor(idx_r))
end select
! ... same for u, v, w ...
! After barrier_arrive(bar_p) + wait:
select case (io)
case (0); pl = p(i1); pr(idx_r) = p(i2); fdx = 1.d0
case (1); call delta4(fdx, p(i1:i2), pl, pr(idx_r)); fdx = wiggle_detector(p(idx-1:idx+2))
case (2); call delta6(fdx, p(i1:i2), pl, pr(idx_r)); fdx = wiggle_detector(p(idx-1:idx+2))
end select
```

The `call interp(id_accuracy, ...)` call is removed and replaced by the 5 per-variable sections above. The guard condition (`io+1 <= i .and. i <= nx-(io+1)`) is preserved in the early-exit path.

The second `call syncthreads()` (line 125 in x_in, before the smem reuse) and the subsequent smem-reuse block (lines 126-133) remain **after** the `#endif`, since they are used in both TMA and non-TMA paths. For the TMA path: `rhol, ul, vl, wl, pl` are already in registers from the per-variable delta calls, so the smem reuse (`rho(idx) = rhol` etc.) still works correctly.

## Verification

1. Build in `3D_solver/NSTGV` on H100/GH200 (arch ≥ 900): `make clean && make`
2. Run `bash calc.sh` — simulation completes without NaN/error
3. Diff VTK output against non-TMA baseline for numerical identity
4. `ncu` profile: per-variable TMA loads visible, MUSCL on rho overlaps u/v/w/p transfers in timeline
