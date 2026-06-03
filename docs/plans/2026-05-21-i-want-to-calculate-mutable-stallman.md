# Plan: Overlapping Computation and MPI Communication (z-decomposition)

## Context

The current solver runs with **no inter-GPU decomposition during time integration**. Each even-ranked MPI process owns the full 3D domain on a single GPU; odd ranks are idle except for I/O. The `exchange_cyclic` subroutine in `calc_para.f90` exists for x-direction decomposition but is never called from the time loop.

The goal is to introduce proper **z-direction domain decomposition** with **communication–computation overlap**: compute G-fluxes (and x/y fluxes) for the interior z-slab while the boundary ghost-cell data is in flight over MPI.

---

## Key Facts (from code reading)

| Item | Finding |
|------|---------|
| Array layout | `Q(nx, 5, ny, nz)` — z is 4th (least contiguous in Fortran) |
| Ghost cell width | `overlap` = 1 (2nd), 2 (4th), 3 (6th order) |
| Current z-BCs | `id_bc_z = .false.` (periodic) — z-decomposition fits naturally |
| Existing x-exchange | `flatten/reconstruct` in `calc_para.f90` pack/unpack x-boundaries |
| MPI rank layout | Even ranks compute (1 GPU each), odd ranks idle — **must change** |
| Pinned buffers | Host buffers allocated inside `exchange_cyclic` each call — should be pre-allocated |

### What "interior" and "halo" mean here

With z-decomposition, each rank holds `Q(nx, 5, ny, nz_local + 2*overlap)`:
- **Ghost cells**: z-indices `[1 .. overlap]` (from rank-1) and `[nz_local+overlap+1 .. nz_local+2*overlap]` (from rank+1)
- **Physical cells**: z-indices `[overlap+1 .. nz_local+overlap]`
- **Halo send data**: physical cells at z = `[overlap+1 .. 2*overlap]` (left edge) and `[nz_local+1 .. nz_local+overlap]` (right edge)

Flux computation at a z-face of index `k` in G needs Q at `k` and `k+1` (2nd order) up to `k±3` (6th order). All E- and F-flux kernels also read Q at `k±1` in z for viscous terms. Thus the **boundary slab** of fluxes that depends on ghost cells has depth equal to `overlap` in z.

---

## Four Design Proposals

---

### Design A — Non-blocking MPI with z-range kernel split *(recommended starting point)*

**Philosophy:** Minimal change to kernel code; pure MPI non-blocking; split each kernel launch into a halo-slab launch and an interior-slab launch using a z-offset argument.

**Sequence (one RK stage):**
```
1. set_bc_local (x, y BCs only — no change)
2. flatten_z_lo / flatten_z_hi  → D2H of boundary Q into pinned host buffers
3. MPI_Isend(send_lo → rank-1),  MPI_Isend(send_hi → rank+1)
   MPI_Irecv(recv_lo ← rank-1),  MPI_Irecv(recv_hi ← rank+1)   [4 non-blocking calls]
4. calc_EFG_interior_z(...)       [GPU: z-faces [overlap+1 .. nz_local-overlap], no ghost needed]
5. MPI_Waitall(4 requests)
6. H2D recv_lo → Qr1d_lo,  H2D recv_hi → Qr1d_hi
7. reconstruct_z(...)             [GPU: fill ghost cells]
8. calc_EFG_halo_z(...)           [GPU: z-faces [1..overlap] and [nz_local+1..nz_local+overlap]]
9. calc_step(...)                 [GPU: full domain, uses combined E,F,G]
```

**Kernel split mechanism:** Add an integer `klo` argument and launch with a restricted grid in z:
```fortran
! Interior G-flux (avoids ghost-dependent faces)
blocksG_int%z = (nz_local - 2*overlap + threadsG%z - 1) / threadsG%z
call calc_keep_z_koffset<<<blocksG_int, threadsG>>>(id_accuracy, nx, ny, nz, Q, T, G, klo=overlap+1)

! Halo G-flux (bottom slab)
blocksG_halo%z = (overlap + threadsG%z - 1) / threadsG%z
call calc_keep_z_koffset<<<blocksG_halo, threadsG>>>(id_accuracy, nx, ny, nz, Q, T, G, klo=1)
! (similar for top slab)
```

In the kernel, global z-index becomes `k = (blockIdx%z-1)*blockDim%z + kt + klo - 1`.

**Files to modify:**
- `src/calc_para.f90` — new `flatten_z_lo`, `flatten_z_hi`, `reconstruct_z`; new `exchange_z_nonblocking` with pre-allocated pinned buffers
- `3D_solver/src/calc_flux_base.f90` — split `calc_EFG` into `calc_EFG_interior_z` / `calc_EFG_halo_z` with klo/khi pass-through
- `3D_solver/src/calc_keep_3d.f90`, `calc_slau_3d.f90`, `calc_visc2.f90`, `calc_visc4.f90` — add `klo` argument to z-direction kernels
- `3D_solver/src/calc_time_dev.f90` — restructure RK loop
- `3D_solver/<CASE>/mod_globals.f90` — add `nz_local`, eliminate even/odd rank distinction
- `3D_solver/<CASE>/set.f90` — subdomain IC/BC

**Pros:** Portable, no special hardware, straightforward to reason about.  
**Cons:** D2H and H2D memory copies are serialized (step 2 and 6 are separate); MPI call blocks the CPU thread while the GPU runs interior kernels.

---

### Design B — Two-CUDA-stream pipeline (GPU compute + async H2D overlap)

**Philosophy:** Launch interior GPU computation on stream 0 while performing H2D/D2H memory copies for halo data on stream 1, then issue MPI calls. This additionally overlaps the H2D transfer with interior computation.

**Sequence:**
```
1. set_bc_local
2. [stream 1] flatten_z (D2H async) → cudaEventRecord(ev_d2h)
3. cudaStreamSynchronize(stream 1)  [CPU must have host buffer for MPI_Isend]
4. MPI_Isend / MPI_Irecv  (4 non-blocking)
5. [stream 0] calc_EFG_interior_z  (interior, overlaps with MPI latency on CPU side)
6. MPI_Waitall
7. [stream 1] H2D recv buffers → Qr1d  (async)
8. [stream 0] barrier: cudaStreamWaitEvent(stream 0, ev_h2d)
9. [stream 0] reconstruct_z + calc_EFG_halo_z
10. [stream 0] calc_step
```

**Key detail:** `MPI_Isend`/`MPI_Irecv` are CPU non-blocking calls; the CPU issues them and immediately launches the interior GPU kernel. The GPU runs `calc_EFG_interior_z` while the network/host is handling MPI. True overlap occurs if the GPU kernel runtime > MPI latency.

**Additional change:** Pre-allocate two pinned host buffers (`send_lo`, `send_hi`, `recv_lo`, `recv_hi`) in `pre_calc` to avoid per-call allocation overhead.

**Pros:** Best latency hiding on large grids where GPU kernel >> MPI latency; no extra GPU memory.  
**Cons:** Adds CUDA stream management complexity; Fortran CUDA stream API is limited; requires careful event synchronization.

---

### Design C — GPU-aware MPI (GPUDirect / `id_gpumpi`)

**Philosophy:** Use `id_gpumpi` (already defined in `mod_globals.f90`) to pass device pointers directly to MPI, bypassing H2D and D2H entirely. Communication happens via NVLink or RDMA over InfiniBand.

**Sequence (one RK stage):**
```
1. set_bc_local
2. [GPU] flatten_z → Qs1d_lo, Qs1d_hi  (device arrays, no H2D)
3. MPI_Isend(Qs1d_lo[device], rank-1),  MPI_Isend(Qs1d_hi[device], rank+1)
   MPI_Irecv(Qr1d_lo[device], rank-1),  MPI_Irecv(Qr1d_hi[device], rank+1)
4. [GPU stream 0] calc_EFG_interior_z
5. MPI_Waitall
6. [GPU] reconstruct_z (Qr1d → ghost cells of QJ)
7. [GPU] calc_EFG_halo_z
8. [GPU] calc_step
```

**Activation:** Conditioned on `kind(id_gpumpi) == 4` (use existing kind-dispatch pattern):
```fortran
if (kind(id_gpumpi) == 4) then
  call MPI_Isend(Qs1d_lo, ..., rank-1, ...)   ! device pointer
else
  stat = cudaMemcpy(Qs_lo_pinned, Qs1d_lo, ..., cudaMemcpyDeviceToHost)
  call MPI_Isend(Qs_lo_pinned, ..., rank-1, ...)
endif
```

**Pros:** Eliminates D2H/H2D bottleneck entirely; maximum bandwidth utilization; cleanest code path.  
**Cons:** Requires NVLink (multi-GPU node) or InfiniBand + GPUDirect RDMA; MPI library must support `CUDA_AWARE_MPI`; additional `#ifdef` or runtime check needed.

---

### Design D — Send-ahead lookahead (post-step packing)

**Philosophy:** Pack and send halo data at the **end** of each RK stage (after `calc_step`), so that MPI is in flight during the **next stage's** interior computation — shifting communication by one stage.

**Modified RK stage loop:**
```fortran
! Stage n (assuming recv from stage n-1 already done):
call calc_EFG_interior_z(...)        ! interior flux (no ghost needed)
call MPI_Waitall(recv_req)           ! wait for ghost cells from stage n-1 comm
call reconstruct_z(...)              ! unpack ghost cells
call calc_EFG_halo_z(...)            ! boundary flux (now ghost cells valid)
call calc_step(...)                  ! update Q

! Pack and fire for next stage:
call flatten_z(...)                  ! D2H
call MPI_Isend / MPI_Irecv(send_req, recv_req)  ! non-blocking
! → loop back to next stage with MPI in flight
```

**Bootstrap issue:** Before stage 1 of a new time step, an initial blocking exchange is needed to seed the pipeline. This can be a single `MPI_Sendrecv` call.

**Pros:** Maximizes overlap: interior flux computation (longest GPU kernel) fully overlaps MPI.  
**Cons:** Requires managing separate request arrays across stage boundaries; more complex bookkeeping; one extra blocking exchange per time step for bootstrap.

---

## Common Infrastructure (all designs)

### New z-direction pack/unpack kernels (`calc_para.f90`)

```fortran
! Pack: extract physical boundary slabs for sending
subroutine flatten_z(nx, ny, nz, nz_local, overlap, Q, Q1d_lo, Q1d_hi)
  ! Q1d_lo ← Q(:,:,:, overlap+1 : 2*overlap)      ! send to rank-1
  ! Q1d_hi ← Q(:,:,:, nz_local+1 : nz_local+overlap)  ! send to rank+1

! Unpack: fill ghost cell slabs from received data
subroutine reconstruct_z(nx, ny, nz, overlap, Q1d_lo, Q1d_hi, Q)
  ! Q(:,:,:, 1:overlap) ← Q1d_lo  (from rank-1)
  ! Q(:,:,:, nz-overlap+1:nz) ← Q1d_hi (from rank+1)
```

### Rank topology change

Current: even ranks compute, odd ranks idle (x-decomposition artifact).  
New: every rank computes 1 z-slab. Total compute ranks = `nranks`. Periodic z wraps rank 0 ↔ rank nranks-1.

```fortran
rank_lo = mod(myrank - 1 + nranks, nranks)  ! rank-1 in z
rank_hi = mod(myrank + 1, nranks)            ! rank+1 in z
```

### Grid size per rank

`nz_local = nz_global / nranks` (must divide evenly, enforced in `set.f90`).  
Local array: `QJ(nx, 5, ny, nz_local + 2*overlap)`.

Update block grid:
```fortran
blocksG%z = (nz_local + threadsG%z - 1) / threadsG%z
```

---

## Files to Modify

| File | Change |
|------|--------|
| `src/calc_para.f90` | Add `flatten_z`, `reconstruct_z`, `exchange_z_nb` (non-blocking); pre-allocated pinned buffers |
| `3D_solver/src/calc_time_dev.f90` | New RK loop variant with overlap for each of Design A–D; add `RungeKutta_4th_zdec` |
| `3D_solver/src/calc_flux_base.f90` | Add `calc_EFG_interior_z` / `calc_EFG_halo_z` variants; update `calc_EFG_z_interior_visc` / `calc_EFG_z_halo_visc` for 4th-order |
| `3D_solver/src/calc_keep_3d.f90` | Add `klo` offset argument to z-direction kernel |
| `3D_solver/src/calc_slau_3d.f90` | Same |
| `3D_solver/src/calc_visc2.f90` | Same |
| `3D_solver/src/calc_visc4.f90` | Add `calc_Ev4_koff`, `calc_Fv4_koff`, `calc_Gv4_koff` |
| `3D_solver/src/load_smem_visc4.f90` | Add `load_smem_visc4_z_koff` |
| `3D_solver/<CASE>/mod_globals.f90` | Add `nz_local`, `nranks_z`; remove even/odd rank logic |
| `3D_solver/<CASE>/set.f90` | Partition z-grid per rank; subdomain IC/BC |
| `src/preprocess.f90` | Allocate pinned halo buffers; compute z-subdomain offsets |

---

## Additional Subroutines: z-decomposition 4th-order support

These five new subroutines extend the z-decomposition communication-overlap path to 4th-order viscous accuracy. They are added alongside the existing 2nd-order `_koff` variants and the `RungeKutta_3rd_zdec`.

---

### 1. `load_smem_visc4_z_koff` — [3D_solver/src/load_smem_visc4.f90](3D_solver/src/load_smem_visc4.f90)

**Purpose:** Load shared memory for 4th-order viscous z-flux kernels when launched with a z-offset (koff) to cover only a sub-range of z-slabs.

**Pattern:** Mirrors `load_smem_visc2_z_koff` ([3D_solver/src/load_smem_visc2.f90:126](3D_solver/src/load_smem_visc2.f90#L126)). The only change from `load_smem_visc4_z` ([3D_solver/src/load_smem_visc4.f90:191](3D_solver/src/load_smem_visc4.f90#L191)) is:

```fortran
! In load_smem_visc4_z:
k_base = (blockIdx%z-1)*blockDim%z

! In load_smem_visc4_z_koff (add k_lo parameter):
k_base = (blockIdx%z-1)*blockDim%z + k_lo - 1
```

**New parameter:** `integer, intent(in), value :: k_lo` — appended after the shared-array arguments (same position as in `load_smem_visc2_z_koff`).

**Public list update:** Add `load_smem_visc4_z_koff` to the `public` statement at line 7.

---

### 2. `calc_Ev4_koff` — [3D_solver/src/calc_visc4.f90](3D_solver/src/calc_visc4.f90)

**Purpose:** 4th-order viscous E-flux for a restricted z-range `[k_lo, k_hi]`, enabling the interior/halo split of `calc_EFG_z_interior_visc` / `calc_EFG_z_halo_visc`.

**Pattern:** Mirrors `calc_Ev2_koff` ([3D_solver/src/calc_visc2.f90:481](3D_solver/src/calc_visc2.f90#L481)). Differences from `calc_Ev4` ([3D_solver/src/calc_visc4.f90:17](3D_solver/src/calc_visc4.f90#L17)):

| Line/item | `calc_Ev4` | `calc_Ev4_koff` |
|-----------|-----------|----------------|
| Signature | `(nx, ny, nz, dx, dy, dz, Q, T, mu, E)` | `(nx, ny, nz, dx, dy, dz, Q, T, mu, E, k_lo, k_hi)` |
| k-index | `k = (blockIdx%z-1)*blockDim%z + kt + 1` | `k = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt` |
| Bounds check | `if (... .or. nz-1 < k) return` | `if (... .or. k_hi < k) return` |
| Loader call | `load_smem_visc4_x(it, jt, kt, j, k, ...)` — k already shifted via formula above | unchanged (k is passed in) |

The high-order interior guard `if (3 <= i ... .and. 3 <= k .and. k <= nz-3)` remains as-is; it naturally clips to the correct range once k is shifted.

**Public list update:** Add `calc_Ev4_koff` to the `public` statement at line 10.

---

### 3. `calc_Fv4_koff` — [3D_solver/src/calc_visc4.f90](3D_solver/src/calc_visc4.f90)

**Purpose:** 4th-order viscous F-flux for a restricted z-range `[k_lo, k_hi]`.

**Pattern:** Same as `calc_Ev4_koff` above, applied to `calc_Fv4` ([3D_solver/src/calc_visc4.f90:231](3D_solver/src/calc_visc4.f90#L231)). Mirrors `calc_Fv2_koff` ([3D_solver/src/calc_visc2.f90:542](3D_solver/src/calc_visc2.f90#L542)):

| Item | Change |
|------|--------|
| Signature | Add `k_lo, k_hi` integer parameters |
| k-index | `k = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt` (was `+ kt + 1`) |
| Bounds check | `k_hi < k` (was `nz-1 < k`) |
| Loader `load_smem_visc4_y` | Called with `i, k` — k is already shifted, no loader change needed |

**Public list update:** Add `calc_Fv4_koff` to the `public` statement at line 10.

---

### 4. `calc_Gv4_koff` — [3D_solver/src/calc_visc4.f90](3D_solver/src/calc_visc4.f90)

**Purpose:** 4th-order viscous G-flux for a restricted z-range `[k_lo, k_hi]`.

**Pattern:** Mirrors `calc_Gv2_koff` ([3D_solver/src/calc_visc2.f90:603](3D_solver/src/calc_visc2.f90#L603)). This is the most involved because `calc_Gv4` calls `load_smem_visc4_z` with an internal k_base — the koff variant must call `load_smem_visc4_z_koff` instead:

| Item | `calc_Gv4` ([line 447](3D_solver/src/calc_visc4.f90#L447)) | `calc_Gv4_koff` |
|------|---------|----------------|
| Signature | `(nx, ny, nz, dx, dy, dz, Q, T, mu, G)` | `(nx, ny, nz, dx, dy, dz, Q, T, mu, G, k_lo, k_hi)` |
| Loader call | `load_smem_visc4_z(it, jt, kt, i, j, ...)` | `load_smem_visc4_z_koff(it, jt, kt, i, j, ..., k_lo)` |
| k-index (after load) | `k = (blockIdx%z-1)*blockDim%z + kt` | `k = (blockIdx%z-1)*blockDim%z + k_lo - 1 + kt` |
| Bounds check | `if (... .or. nz-1 < k) return` | `if (... .or. k_hi < k) return` |

All other computation (stress/heat-flux body) is identical to `calc_Gv4`.

**Public list update:** Add `calc_Gv4_koff` to the `public` statement at line 10.

---

### 5. `RungeKutta_4th_zdec` — [3D_solver/src/calc_time_dev.f90](3D_solver/src/calc_time_dev.f90)

**Purpose:** 4-stage classical RK4 time-stepping with z-direction domain decomposition and communication-computation overlap. All MPI ranks compute. Mirrors `RungeKutta_3rd_zdec` ([line 444](3D_solver/src/calc_time_dev.f90#L444)) but uses the 4-stage RK4 formula and a residual-accumulation array `Rs`.

**Interface disambiguation:** Uses the same `integer(8)` kind for `id_RungeKutta` as `RungeKutta_3rd_zdec`, but `integer(4)` for `id_rescale` (vs `integer(2)` in `RungeKutta_3rd_zdec`). The Fortran generic interface resolves by argument type combination.

```fortran
subroutine RungeKutta_4th_zdec(id_RungeKutta, id_rescale, myrank, mygpu, nx, ny, nz, &
                                 x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
  integer(8), intent(in) :: id_RungeKutta   ! kind=8 → zdec mode
  integer(4), intent(in) :: id_rescale      ! kind=4 → 4th-order (distinguishes from RK3_zdec)
```

**Memory:** Same as `RungeKutta_3rd_zdec` but adds `Rs(nx-2,5,ny-2,nz-2)` for residual accumulation (like `RungeKutta_4th` at [line 279](3D_solver/src/calc_time_dev.f90#L279)).

**4-stage loop structure** (one timestep, analogous to `RungeKutta_4th` [lines 297-316](3D_solver/src/calc_time_dev.f90#L297-L316) but with MPI overlap):

```
Stage 1: QJ → QJs (weight=1 to Rs, step=0.5*dt)
  start_exchange_z(QJ)
  calc_EFG_z_interior(QJ)
  finish_exchange_z → QJ ghost cells
  calc_EFG_z_halo(QJ)
  calc_step<<<>>>(nx,ny,nz, 0.5d0, 1.d0, ..., E,F,G, QJ, QJs, Rs)   ! QJs = Q^n + 0.5*dt*k1
  set_bc(QJs)

Stage 2: QJs → QJs (weight=2 to Rs, step=0.5*dt)
  start_exchange_z(QJs)
  calc_EFG_z_interior(QJs)
  finish_exchange_z → QJs ghost cells
  calc_EFG_z_halo(QJs)
  calc_step<<<>>>(nx,ny,nz, 0.5d0, 2.d0, ..., E,F,G, QJ, QJs, Rs)   ! QJs = Q^n + 0.5*dt*k2
  set_bc(QJs)

Stage 3: QJs → QJs (weight=2 to Rs, step=1.0*dt)
  start_exchange_z(QJs)
  calc_EFG_z_interior(QJs)
  finish_exchange_z → QJs ghost cells
  calc_EFG_z_halo(QJs)
  calc_step<<<>>>(nx,ny,nz, 1.0d0, 2.d0, ..., E,F,G, QJ, QJs, Rs)   ! QJs = Q^n + dt*k3
  set_bc(QJs)

Stage 4: QJs → QJ (final assembly, weight=1, uses calc_step4)
  start_exchange_z(QJs)
  calc_EFG_z_interior(QJs)
  finish_exchange_z → QJs ghost cells
  calc_EFG_z_halo(QJs)
  calc_step4<<<>>>(nx,ny,nz, ..., E,F,G, Rs, QJ)                     ! QJ = Q^n + (dt/6)*(k1+2k2+2k3+k4)
  set_bc(QJ)
  Rs = 0.d0   ! reset for next timestep
```

The `calc_step` and `calc_step4` kernel calls are identical to those in `RungeKutta_4th` ([lines 299-314](3D_solver/src/calc_time_dev.f90#L299-L314)).

**Interface update:** Add `RungeKutta_4th_zdec` to the `RungeKutta` generic interface at [line 23](3D_solver/src/calc_time_dev.f90#L23).

---

### 6. Update `calc_EFG_z_interior_visc` and `calc_EFG_z_halo_visc` for 4th-order — [3D_solver/src/calc_flux_base.f90](3D_solver/src/calc_flux_base.f90)

The existing `calc_EFG_z_interior_visc` ([line 404](3D_solver/src/calc_flux_base.f90#L404)) and `calc_EFG_z_halo_visc` ([line 438](3D_solver/src/calc_flux_base.f90#L438)) currently only call `calc_Ev2_koff` / `calc_Fv2_koff` / `calc_Gv2_koff`. To support 4th-order accuracy in zdec mode, add a dispatch branch using `kind(id_accuracy)`:

```fortran
! In calc_EFG_z_interior_visc (after calc_conv):
if (kind(id_accuracy) >= 4) then
  call calc_Ev4_koff<<<blocksEv_int,threadsEv>>>(nx,ny,nz,inv_dx,inv_dy,inv_dz,Q,T,mu,E, &
       overlap_fb+2, nz-overlap_fb-1)
  call calc_Fv4_koff<<<blocksFv_int,threadsFv>>>(nx,ny,nz,inv_dy,inv_dx,inv_dz,Q,T,mu,F, &
       overlap_fb+2, nz-overlap_fb-1)
  call calc_Gv4_koff<<<blocksGv_int,threadsGv>>>(nx,ny,nz,inv_dx,inv_dy,inv_dz,Q,T,mu,G, &
       overlap_fb+1, nz-overlap_fb-1)
else
  call calc_Ev2_koff<<<...>>>(...) ! existing 2nd-order path
  ...
end if
```

The same pattern applies to `calc_EFG_z_halo_visc` for the lo/hi halo slab launches.

---

## Verification

1. Run `bash test_cicd.sh` from repo root — all 9 Cartesian targets must pass.
2. Compare single-rank z-decomposition output (nranks=1) with current baseline VTK output for NSTGV (TGV) and SBLI cases.
3. Multi-rank: run 4-rank z-decomposition and verify kinetic energy decay curve matches 1-GPU reference.
4. Overlap correctness: temporarily disable interior computation (zero-out interior fluxes) and verify solution diverges — confirms halo+interior split covers all cells.
5. Timing: measure wallclock per RK stage with and without overlap using `nvtx` markers already present.
6. Compile check for new subroutines: enable 4th-order accuracy (`integer(4), parameter :: id_accuracy = 0`) and `RungeKutta_4th_zdec` (`integer(8), parameter :: id_RungeKutta = 0; integer(4), parameter :: id_rescale = 0`) in one test case's `mod_globals.f90`, then run `bash test_cicd.sh`.
